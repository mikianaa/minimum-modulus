"""Exact residual-partition construction search; no sampling.

Each resource is placed at a maximum uncovered cofactor of a reduced
multi-valued decision diagram. Integer-valued marginal diagrams eliminate
unconstrained digits by summation, so every gain and overlap is exact.
Only the independent construction_kernel verifier can certify success.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import json
from math import prod
from pathlib import Path
import time

from construction_kernel import resources, capacity, labels, verify
from owens_branch_search import primes_through
from search_large_minimum import depths_for

ROOT = Path(__file__).resolve().parents[1]
ZERO, ONE = -1, -2


class DiagramLimit(Exception):
    pass


class Diagram:
    def __init__(self, domains, limit=200000, deadline=float('inf')):
        self.domains = tuple(domains)
        self.nodes = []
        self.unique = {}
        self.limit, self.deadline = limit, deadline
        self.suffix = [1]*(len(domains)+1)
        for i in range(len(domains)-1,-1,-1):
            self.suffix[i] = self.suffix[i+1]*domains[i]

    def node(self, v, children):
        children = tuple(children)
        assert len(children) == self.domains[v]
        if len(set(children)) == 1:
            return children[0]
        key = (v,children)
        if key in self.unique:
            return self.unique[key]
        if len(self.nodes) >= self.limit:
            raise DiagramLimit('diagram node limit')
        if len(self.nodes) % 1024 == 0 and time.monotonic() >= self.deadline:
            raise DiagramLimit('time limit')
        idx = len(self.nodes)
        self.nodes.append(key)
        self.unique[key] = idx
        return idx

    @lru_cache(None)
    def count(self, node, start=0):
        if node < 0:
            return (-node-1)*self.suffix[start]
        v,children = self.nodes[node]
        assert v >= start
        return (self.suffix[start]//self.suffix[v]) * sum(self.count(c,v+1) for c in children)

    def erase(self, root, condition):
        """Return the exact Boolean function root AND NOT(condition)."""
        condition = tuple(sorted(condition.items()))
        @lru_cache(None)
        def walk(node,k):
            if node == ZERO or k == len(condition):
                return ZERO
            wanted,digit = condition[k]
            v = self.nodes[node][0] if node >= 0 else len(self.domains)
            if wanted <= v:
                cs = list(self.nodes[node][1]) if wanted == v else [node]*self.domains[wanted]
                cs[digit] = walk(cs[digit],k+1)
                return self.node(wanted,cs)
            return self.node(v,[walk(c,k) for c in self.nodes[node][1]])
        return walk(root,0)

    def maximum_cofactor(self, root, selected):
        """Maximize exact residual cardinality for the selected digit labels."""
        selected = frozenset(selected)
        q = Diagram(self.domains,self.limit,self.deadline)
        # Products over eliminated variables only.
        drop = [1]*(len(self.domains)+1)
        for i in range(len(self.domains)-1,-1,-1):
            drop[i] = drop[i+1]*(1 if i in selected else self.domains[i])

        @lru_cache(None)
        def scale(node,k):
            if k == 1 or node == ZERO:
                return node
            if node < 0:
                return -((-node-1)*k)-1
            v,cs = q.nodes[node]
            return q.node(v,[scale(c,k) for c in cs])

        @lru_cache(None)
        def add(a,b):
            if a == ZERO:
                return b
            if b == ZERO:
                return a
            if a == b:
                return scale(a,2)
            if a < 0 and b < 0:
                return a+b+1
            va = q.nodes[a][0] if a >= 0 else len(self.domains)
            vb = q.nodes[b][0] if b >= 0 else len(self.domains)
            v = min(va,vb)
            ca = q.nodes[a][1] if va == v else [a]*self.domains[v]
            cb = q.nodes[b][1] if vb == v else [b]*self.domains[v]
            return q.node(v,[add(x,y) for x,y in zip(ca,cb)])

        @lru_cache(None)
        def marginal(node,start):
            if node < 0:
                return -((-node-1)*drop[start])-1
            v,cs = self.nodes[node]
            children = [marginal(c,v+1) for c in cs]
            if v in selected:
                result = q.node(v,children)
            else:
                result = ZERO
                for child,multiplicity in Counter(children).items():
                    result = add(result,scale(child,multiplicity))
            return scale(result,drop[start]//drop[v])

        out = marginal(root,0)
        @lru_cache(None)
        def best(node):
            if node < 0:
                return -node-1
            return max(best(c) for c in q.nodes[node][1])
        value = best(out)
        assignment = {v:0 for v in selected}
        while out >= 0:
            v,cs = q.nodes[out]
            digit = max(range(len(cs)),key=lambda j:best(cs[j]))
            assignment[v] = digit
            out = cs[digit]
        return value,assignment,len(q.nodes)

    def compact(self,root):
        new = Diagram(self.domains,self.limit,self.deadline)
        @lru_cache(None)
        def clone(node):
            if node < 0:
                return node
            v,cs = self.nodes[node]
            return new.node(v,[clone(c) for c in cs])
        return new,clone(root)

    def whole_cofactor(self,root,selected):
        """Find a whole cylinder inside the residual, or prove none exists.

        Universally quantify unselected digits with Boolean AND. This
        deliberately forbids overlap in the kernel; finite tail completion
        may still overlap. It is a restricted construction policy.
        """
        selected = frozenset(selected)
        q = Diagram(self.domains,self.limit,self.deadline)
        @lru_cache(None)
        def meet(a,b):
            if a == ZERO or b == ZERO:
                return ZERO
            if a == ONE or a == b:
                return b
            if b == ONE:
                return a
            va,ca = q.nodes[a]
            vb,cb = q.nodes[b]
            v = min(va,vb)
            if va != v:
                ca = [a]*q.domains[v]
            if vb != v:
                cb = [b]*q.domains[v]
            return q.node(v,[meet(*sorted((x,y))) for x,y in zip(ca,cb)])
        @lru_cache(None)
        def project(n):
            if n < 0:
                assert n in (ZERO,ONE)
                return n
            v,cs = self.nodes[n]
            children = [project(c) for c in cs]
            if v in selected:
                return q.node(v,children)
            out = ONE
            for c in set(children):
                out = meet(*sorted((out,c)))
                if out == ZERO:
                    break
            return out
        projected = project(root)
        assignment = {v:0 for v in selected}
        if projected == ZERO:
            return 0,assignment,len(q.nodes)
        while projected >= 0:
            v,cs = q.nodes[projected]
            digit = next(i for i,c in enumerate(cs) if c != ZERO)
            assignment[v] = digit
            projected = cs[digit]
        assert projected == ONE
        volume = self.suffix[0]//prod(self.domains[v] for v in selected)
        return volume,assignment,len(q.nodes)

    def digitwise_cofactor(self,root,selected):
        """Structural coordinate choice with exact gains, not a global max.

        Optimize one selected digit at a time, condition on it, and repeat.
        Only one-axis marginals are formed, avoiding large integer-valued
        diagrams. This is deterministic and uses no sampled states.
        """
        q,current = self.compact(root)
        assignment = {}
        def cofactor(node,v,digit):
            @lru_cache(None)
            def walk(n):
                if n < 0 or q.nodes[n][0] > v:
                    return n
                k,cs = q.nodes[n]
                if k == v:
                    return cs[digit]
                return q.node(k,[walk(c) for c in cs])
            return walk(node)
        for v in sorted(selected):
            @lru_cache(None)
            def histogram(node,start):
                if node < 0 or q.nodes[node][0] > v:
                    c = q.count(node,start)//q.domains[v]
                    return (c,)*q.domains[v]
                k,cs = q.nodes[node]
                factor = q.suffix[start]//q.suffix[k]
                if k == v:
                    return tuple(factor*q.count(c,k+1) for c in cs)
                hist = [histogram(c,k+1) for c in cs]
                return tuple(factor*sum(h[j] for h in hist) for j in range(q.domains[v]))
            values = histogram(current,0)
            digit = max(range(q.domains[v]),key=lambda d:values[d])
            assignment[v] = digit
            current = cofactor(current,v,digit)
        gain = q.count(current)//prod(q.domains[v] for v in selected)
        nodes = len(q.nodes)
        q.count.cache_clear()
        return gain,assignment,nodes


def coordinates(ps,ds,layout='prime-first'):
    items = [(i,k,p if k < E-1 else p-1) for i,(p,E) in enumerate(zip(ps,ds)) for k in range(E)]
    if layout == 'large-prime-first':
        items.sort(key=lambda x:(-x[0],x[1]))
    elif layout == 'digit-level':
        items.sort(key=lambda x:(x[1],x[0]))
    domains = [d for i,k,d in items]
    index = {(i,k):v for v,(i,k,d) in enumerate(items)}
    blocks = [[index[i,k] for k in range(E)] for i,E in enumerate(ds)]
    return domains,blocks


def condition_from_row(row,ps,ds,blocks):
    out = {}
    for p,E,block,e,r in zip(ps,ds,blocks,row['exponents'],row['residues']):
        if e == E:
            r,d0 = divmod(r,p-1)
            out[block[-1]] = d0
            fixed = E-1
        else:
            fixed = e
        for k in range(fixed):
            out[block[k]] = r % p
            r //= p
    return out


def row_from_condition(es,assignment,ps,ds,blocks):
    rs = []
    for p,E,block,e in zip(ps,ds,blocks,es):
        fixed = E-1 if e == E else e
        a = sum(assignment[block[k]]*p**k for k in range(fixed))
        rs.append(a*(p-1)+assignment[block[-1]] if e == E else a)
    return {'exponents':list(es),'residues':rs}


def search(minimum,largest,cap,seconds,node_limit,order,placement='maximum',layout='prime-first',depths=None):
    start = time.monotonic()
    ps = primes_through(largest)
    ds = depths_for(ps,minimum) if depths is None else list(depths)
    if len(ds) != len(ps) or any(type(e) is not int or e < 1 for e in ds):
        raise ValueError('One positive integer depth is required for every prime')
    es = list(resources(ps,ds,minimum,cap))
    mass = capacity(ps,ds,es)
    base = lambda ex:prod(p**e for p,e in zip(ps,ex))
    volume = lambda ex:prod(labels(p,E,e) for p,E,e in zip(ps,ds,ex))
    if order == 'base':
        es.sort(key=lambda ex:(base(ex),volume(ex)))
    elif order == 'mass':
        es.sort(key=lambda ex:(volume(ex),base(ex)))
    else:
        es.sort(key=lambda ex:(max(p for p,e in zip(ps,ex) if e),base(ex)))
    result = {'minimum':minimum,'primes':ps,'depths':ds,
              'maximum_resource_base':cap,'resource_count':len(es),'capacity':str(mass),
              'order':order,'placement':placement,'layout':layout,'node_limit':node_limit,'seconds_budget':seconds,
              'method':'exact_residual_partition_no_sampling','rows':[],'trace':[]}
    if mass < 1:
        result['status'] = 'EXCLUDED_BY_EXACT_CAPACITY'
        return result
    domains,blocks = coordinates(ps,ds,layout)
    diagram = Diagram(domains,node_limit,start+seconds)
    root = ONE
    total = prod(domains)
    before = total
    outcome = 'RESOURCES_EXHAUSTED'
    last_resource = None
    processed = 0
    try:
        for j,ex in enumerate(es):
            last_resource = list(ex)
            if time.monotonic()-start >= seconds:
                raise DiagramLimit('time limit')
            selected = [v for block,e in zip(blocks,ex) for v in block[:e]]
            chooser = {'maximum':diagram.maximum_cofactor,
                       'digits':diagram.digitwise_cofactor,
                       'disjoint':diagram.whole_cofactor}[placement]
            gain,assignment,query_nodes = chooser(root,selected)
            if gain:
                row = row_from_condition(ex,assignment,ps,ds,blocks)
                assert condition_from_row(row,ps,ds,blocks) == assignment
                new_root = diagram.erase(root,assignment)
                after = diagram.count(new_root)
                assert before-after == gain
                result['rows'].append(row)
                root = new_root
                before = after
            processed = j+1
            if j % 10 == 0 or root == ZERO:
                info = {'resources_processed':processed,'rows_placed':len(result['rows']),
                        'uncovered_fraction':str(F(before,total)),
                        'uncovered_decimal':float(F(before,total)),
                        'diagram_nodes':len(diagram.nodes),'query_nodes':query_nodes,
                        'elapsed_seconds':round(time.monotonic()-start,3)}
                result['trace'].append(info)
                print(json.dumps(info),flush=True)
            if root == ZERO:
                outcome = 'EXACT_PARTITION_COVER_FOUND'
                break
            if len(diagram.nodes) > node_limit//4:
                old = diagram
                diagram,root = old.compact(root)
                old.count.cache_clear()
    except DiagramLimit as exc:
        outcome = str(exc).upper().replace(' ','_')
    result['search_outcome'] = outcome
    result['resources_processed'] = processed
    result['last_resource_considered'] = last_resource
    result['exact_uncovered_fraction'] = str(F(before,total))
    result['exact_uncovered_states'] = before
    result['kernel_states'] = total
    # This separate verifier operates on coordinate cylinders, not the
    # marginal diagrams or overlap objective used by this construction.
    result['verification'] = verify(result,node_limit=300000)
    result['status'] = ('VERIFIED' if result['verification']['status'] == 'VERIFIED'
                        else 'PARTIAL_CONSTRUCTION' if result['verification']['status'] == 'UNCOVERED'
                        else 'UNKNOWN_VERIFICATION_LIMIT')
    if outcome == 'EXACT_PARTITION_COVER_FOUND':
        assert result['status'] == 'VERIFIED'
    result['elapsed_seconds'] = round(time.monotonic()-start,3)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum',type=int,default=100)
    parser.add_argument('--largest-prime',type=int,default=31)
    parser.add_argument('--cap',type=int,default=1000)
    parser.add_argument('--seconds',type=int,default=120)
    parser.add_argument('--nodes',type=int,default=100000)
    parser.add_argument('--order',choices=['base','mass','prime-stage'],default='mass')
    parser.add_argument('--placement',choices=['maximum','digits','disjoint'],default='maximum')
    parser.add_argument('--layout',choices=['prime-first','large-prime-first','digit-level'],default='prime-first')
    parser.add_argument('--output',required=True)
    parser.add_argument('--depths',help='Optional comma-separated finite-kernel depths')
    args = parser.parse_args()
    result = search(args.minimum,args.largest_prime,args.cap,args.seconds,args.nodes,args.order,args.placement,args.layout,
                    None if args.depths is None else list(map(int,args.depths.split(','))))
    (ROOT/args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','trace')},indent=2))
