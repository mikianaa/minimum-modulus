"""Construct arrows using certified coverage in earlier prime inputs.

Own leaves use disjoint exponent families. Borrowed inputs refer only to
earlier largest-prime stages and are not duplicated as new congruences.
The guaranteed progression union is checked exhaustively. Its complement
is not asserted to miss incidental coverage of the expanded construction.
"""
import argparse
from copy import deepcopy
from functools import lru_cache
from itertools import product
import json
from math import gcd,prod
from pathlib import Path
import time

from graded_arrow_blocks import np,check as check_base
from graded_arrow_blocks import supply as absolute_supply
from owens_branch_search import primes_through
from construction_kernel import prime


def make_supply(order,fixed,borrowed,materialize=True):
    nodes=[dict(kind='unit',grade=1)]
    current=[(1,0)]
    for p in order:
        e=fixed.get(p,0)
        free=set(borrowed.get(p,()))
        assert not free or e==0
        digits=[d for d in range(1,p) if d not in free]
        if not digits:
            return None,None
        ordered=sorted(current,key=lambda z:(-z[0],z[1]))
        following=current[:]
        for j in range(1,e+1):
            for grade,child in ordered:
                idx=len(nodes) if materialize else len(following)
                if materialize:
                    nodes.append(dict(kind='fixed',prime=p,exponent=j,
                                      children=[child],grade=grade*p**j))
                following.append((grade*p**j,idx))
        k=len(digits)
        for start in range(0,len(ordered)-k+1,k):
            group=ordered[start:start+k]
            grade=min(g for g,_ in group)*p**(e+1)
            idx=len(nodes) if materialize else len(following)
            if materialize:
                nodes.append(dict(kind='tail',prime=p,exponent=e+1,
                                  children=[n for _,n in group],digits=digits,
                                  borrowed=sorted(free),grade=grade))
            following.append((grade,idx))
        current=following
    return sorted(current,key=lambda z:(-z[0],z[1])),nodes if materialize else None


class Catalog:
    def __init__(self,sources=None,independent=False):
        self.sources=[]
        self.by_prime={}
        self.independent=independent
        for source in sources or []:
            self.append(source)

    def append(self,source):
        p=source['root_prime']
        assert p not in self.by_prime
        self.by_prime[p]=len(self.sources)
        self.sources.append(source)

    @lru_cache(maxsize=250000)
    def contains(self,source_id,node_id,a,m):
        source=self.sources[source_id]
        node=source['nodes'][node_id]
        if node['kind']=='unit':
            return True
        p,e=node['prime'],node['exponent']
        old_a=source['residue']
        if node['kind']=='fixed':
            n=p**e
            return m%n==0 and (a-old_a)%n==0 and self.contains(source_id,node['children'][0],a,m)
        prefix=old_a%p**(e-1)
        if m%p**(e-1) or (a-prefix)%p**(e-1):
            return False
        f=0
        r=m
        while r%p==0:
            f+=1
            r//=p
        if self.independent:
            # Independent CRT-intersection calculation: one level beyond
            # v_p(m) represents every later possible first-nonzero level.
            possible=[d for d in range(1,p) if any(
                (a-(prefix+d*p**(k-1)))%gcd(m,p**k)==0
                for k in range(e,max(e,f+1)+1))]
        else:
            delta=(a-prefix)%p**f
            if delta==0:
                possible=range(1,p)
            else:
                while delta%p==0:
                    delta//=p
                possible=[delta%p]
        own=dict(zip(node['digits'],node['children']))
        for d in possible:
            if d in own:
                if not self.contains(source_id,own[d],a,m):
                    return False
            else:
                if d not in node['borrowed'] or p not in self.by_prime:
                    return False
                earlier=self.by_prime[p]
                assert earlier<source_id
                child=self.sources[earlier]['root_inputs'][d-1]
                if not self.contains(earlier,child,a,m):
                    return False
        return True

    def borrowed(self,a,m):
        result={}
        for s,source in enumerate(self.sources):
            p=source['root_prime']
            ds=[d for d,child in enumerate(source['root_inputs'],1)
                if self.contains(s,child,a,m)]
            if ds:
                result[p]=ds
        return result


def core_base(base):
    base=deepcopy(base)
    base['blocks']=[b for b in base['blocks'] if b is not None and b['root_prime']<=7]
    covered=np.zeros(base['period'],dtype=bool)
    for block in base['blocks']:
        for a in block['residues']:
            covered[a::block['modulus']]=True
    base['guaranteed_uncovered']=int(np.count_nonzero(~covered))
    for key in ['trace','capacity_audit','verification','elapsed_seconds']:
        base.pop(key,None)
    base['method']='core_stages_from_graded_arrow_blocks'
    base['verification']=check_base(base)
    base['status']=base['verification']['status']
    return base,covered


def certify_frozen_probe(base,probe,index=0):
    """Export a saved count witness as a fully audited relative schema."""
    original=base
    base,covered=core_base(original)
    sources=[]
    for block in original['blocks']:
        if block is None or block['root_prime']<=7:
            continue
        assert block['tail'] and block['exponent']==1 and len(block['residues'])==1
        fixed={int(p):e for p,e in block['fixed_exponents'].items()}
        values,nodes=absolute_supply(block['order'],fixed,True)
        for node in nodes:
            if node['kind']=='tail':
                node['digits']=list(range(1,node['prime']))
                node['borrowed']=[]
        q=block['root_prime']
        selected=[i for g,i in values if g>=block['input_threshold']][:q-1]
        a,m=block['residues'][0],block['modulus']
        sources.append(dict(root_prime=q,residue=a,modulus=m,
                            fixed_exponents=block['fixed_exponents'],order=block['order'],
                            borrowed={},nodes=nodes,root_inputs=selected))
        covered[a::m]=True
    witness=probe['best_targets'][index]
    q=probe['root_prime']
    assert witness['sufficient'] and not witness['already_covered_by_prior_block']
    fixed={p:e for p,e in zip(base['core_primes'],witness['fixed_exponents']) if e}
    borrowed={int(p):ds for p,ds in witness['borrowed_digits'].items()}
    values,nodes=make_supply(witness['order'],fixed,borrowed)
    assert len(values)==witness['available_sets'] and len(values)>=q-1
    a,m=witness['residue'],witness['modulus']
    sources.append(dict(root_prime=q,residue=a,modulus=m,
                        fixed_exponents={str(p):e for p,e in fixed.items()},
                        order=witness['order'],borrowed=witness['borrowed_digits'],
                        nodes=nodes,root_inputs=[i for g,i in sorted(values)[:q-1]]))
    covered[a::m]=True
    result=dict(minimum=base['minimum'],base=base,sources=sources,
                derived_progressions=[],guaranteed_uncovered=int(np.count_nonzero(~covered)),
                method='certified_frozen_base_relative_replacement',largest_prime=q,
                new_progression=dict(residue=a,modulus=m),
                recipe='Recorded own nodes and root inputs define all new congruence families. Borrowed inputs reference earlier root stages. Finitize own tail occurrences with distinct auxiliary primes above largest_prime, in source order and preorder traversal, after the base stages.')
    result['verification']=validate(result)
    result['status']=result['verification']['status']
    return result


def best_order(q,fixed,borrowed,minimum):
    ps=primes_through(q-1)
    special=[p for p in ps if p in fixed]
    ordinary=sorted((p for p in ps if p not in fixed),
                    key=lambda p:(p-1-len(borrowed.get(p,())),p))
    orders=[ps,special+ordinary,list(reversed(special))+ordinary]
    best=None
    for order in orders:
        values,_=make_supply(order,fixed,borrowed,False)
        if values is None:
            return None,order
        threshold=(minimum+q-1)//q
        good=sum(g>=threshold for g,_ in values)
        if best is None or good>best[0]:
            best=good,order
    return best


def validate(result):
    """Rebuild supply operations and audit all borrowed inputs by CRT."""
    base=result['base']
    assert result['minimum']==base['minimum']
    assert all(block['root_prime']<=7 for block in base['blocks'])
    check_base(base)
    L=base['period']
    covered=bytearray(L)
    for block in base['blocks']:
        for a in block['residues']:
            covered[a::block['modulus']]=b'\x01'*len(range(a,L,block['modulus']))
    catalog=Catalog(independent=True)
    dependency_checks=0
    for source in result['sources']:
        q=source['root_prime']
        a,m=source['residue'],source['modulus']
        assert prime(q) and q>7 and q not in catalog.by_prime
        assert q<=result['largest_prime']
        assert not catalog.sources or q>catalog.sources[-1]['root_prime']
        assert 0<=a<m and L%m==0
        fixed={int(p):e for p,e in source['fixed_exponents'].items()}
        assert all(prime(p) and type(e) is int and e>0 for p,e in fixed.items())
        assert prod(p**e for p,e in fixed.items())==m
        assert set(fixed)<=set(base['core_primes'])
        order=source['order']
        assert len(order)==len(set(order)) and all(prime(p) and p<q for p in order)
        assert set(fixed)<=set(order)
        borrowed={int(p):ds for p,ds in source['borrowed'].items()}
        for p,ds in borrowed.items():
            assert p in catalog.by_prime and p not in fixed
            assert len(ds)==len(set(ds)) and all(1<=d<p for d in ds)
            old_id=catalog.by_prime[p]
            for d in ds:
                child=catalog.sources[old_id]['root_inputs'][d-1]
                assert catalog.contains(old_id,child,a,m)
                dependency_checks+=1
        # Independently check the serialized operations against the
        # disjoint exponent-copy/group induction; no coverage counts used.
        nodes=source['nodes']
        assert nodes[0]==dict(kind='unit',grade=1)
        active=[0]
        cursor=1
        for p in order:
            e=fixed.get(p,0)
            ranked=sorted(active,key=lambda i:(-nodes[i]['grade'],i))
            following=active[:]
            for j in range(1,e+1):
                for child in ranked:
                    expected=dict(kind='fixed',prime=p,exponent=j,children=[child],
                                  grade=nodes[child]['grade']*p**j)
                    assert nodes[cursor]==expected
                    following.append(cursor)
                    cursor+=1
            free=sorted(borrowed.get(p,[]))
            digits=[d for d in range(1,p) if d not in free]
            k=len(digits)
            assert k>0
            for j in range(len(ranked)//k):
                children=ranked[j*k:(j+1)*k]
                expected=dict(kind='tail',prime=p,exponent=e+1,children=children,
                              digits=digits,borrowed=free,
                              grade=min(nodes[c]['grade'] for c in children)*p**(e+1))
                assert nodes[cursor]==expected
                following.append(cursor)
                cursor+=1
            active=following
        assert cursor==len(nodes)
        selected=source['root_inputs']
        assert len(selected)==q-1 and len(set(selected))==q-1
        assert set(selected)<=set(active)
        assert all(q*nodes[i]['grade']>=result['minimum'] for i in selected)
        catalog.append(source)
        covered[a::m]=b'\x01'*len(range(a,L,m))
    for derived in result['derived_progressions']:
        p,a,m=derived['source_prime'],derived['residue'],derived['modulus']
        assert p in catalog.by_prime and 0<=a<m and L%m==0 and m%p!=0
        s=catalog.by_prime[p]
        assert all(catalog.contains(s,c,a,m) for c in catalog.sources[s]['root_inputs'])
        covered[a::m]=b'\x01'*len(range(a,L,m))
    gap=covered.count(0)
    assert gap==result['guaranteed_uncovered']
    catalog.contains.cache_clear()
    return dict(status='VERIFIED_RELATIVE_CONSTRUCTION_SCHEMA' if gap==0 else 'VERIFIED_PARTIAL_RELATIVE_SCHEMA',
                guaranteed_uncovered=gap,period=L,sources=len(result['sources']),
                independently_checked_borrowed_inputs=dependency_checks,
                first_projected_gap=covered.index(0) if gap else None,
                gap_scope='Not an uncovered witness for all expanded congruences')


def search(base,largest=113,maximum_modulus=720,per_shape=4):
    start=time.monotonic()
    base,covered=core_base(base)
    L=base['period']
    shapes=[]
    ps=base['core_primes']
    for es in product(*(range(e+1) for e in base['core_exponents'])):
        m=prod(p**e for p,e in zip(ps,es))
        if m<=maximum_modulus:
            shapes.append((m,es))
    catalog=Catalog()
    derived=[]
    trace=[]
    for q in primes_through(largest):
        if q<=7:
            continue
        while True:
            best=None
            for m,es in shapes:
                hist=(~covered).reshape(L//m,m).sum(axis=0,dtype=np.int64)
                for aa in np.argsort(-hist,kind='stable')[:per_shape]:
                    a=int(aa)
                    gain=int(hist[a])
                    if not gain or best is not None and gain<best['gain']:
                        continue
                    borrowed=catalog.borrowed(a,m)
                    already=next((p for p,ds in borrowed.items() if len(ds)==p-1),None)
                    fixed={p:e for p,e in zip(ps,es) if e}
                    good,order=best_order(q,fixed,borrowed,base['minimum']) if already is None else (None,[])
                    if already is None and good<q-1:
                        continue
                    proposal=dict(residue=a,modulus=m,fixed=fixed,borrowed=borrowed,
                                  gain=gain,available=good,order=order,already=already)
                    if best is None or (gain,-m,a)>(best['gain'],-best['modulus'],best['residue']):
                        best=proposal
            if best is None:
                break
            a,m=best['residue'],best['modulus']
            if best['already'] is not None:
                derived.append(dict(source_prime=best['already'],residue=a,modulus=m))
                covered[a::m]=True
                if np.all(covered):
                    break
                continue
            values,nodes=make_supply(best['order'],best['fixed'],best['borrowed'])
            # Prefer broad inputs to retain useful coverage in neighboring
            # branches; the minimum-modulus threshold remains mandatory.
            good=sorted(((g,i) for g,i in values if q*g>=base['minimum']))
            assert len(good)==best['available'] and len(good)>=q-1
            source=dict(root_prime=q,residue=a,modulus=m,
                        fixed_exponents={str(p):e for p,e in best['fixed'].items()},
                        order=best['order'],borrowed={str(p):ds for p,ds in best['borrowed'].items()},
                        nodes=nodes,root_inputs=[i for _,i in good[:q-1]],
                        available_inputs=len(good))
            catalog.append(source)
            covered[a::m]=True
            info=dict(root_prime=q,residue=a,modulus=m,gain=best['gain'],
                      available_inputs=len(good),guaranteed_uncovered=int(np.count_nonzero(~covered)),
                      borrowed_inputs=sum(len(ds) for ds in best['borrowed'].values()),
                      elapsed_seconds=round(time.monotonic()-start,3))
            trace.append(info)
            print(json.dumps(info),flush=True)
            break
        if np.all(covered):
            break
    result=dict(minimum=base['minimum'],method='branch_dependent_relative_arrow_supplies',
                base=base,sources=catalog.sources,derived_progressions=derived,
                trace=trace,guaranteed_uncovered=int(np.count_nonzero(~covered)),
                largest_prime=largest,maximum_target_modulus=maximum_modulus,
                targets_per_shape=per_shape,elapsed_seconds=time.monotonic()-start,
                recipe='Each source is a largest-prime tail with its recorded root inputs in digit order 1,...,q-1. Own fixed/tail nodes use the recorded target prefix. Borrowed digits reference the earlier source with that root prime and that digit; they add no new leaves. Finitize every own tail occurrence with a distinct auxiliary prime larger than largest_prime, allocated by base stages, sources, then preorder. Base and source exponent families are disjoint by largest-prime/exponent and the checked supply induction.')
    catalog.contains.cache_clear()
    result['verification']=validate(result)
    result['status']=result['verification']['status']
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base')
    parser.add_argument('--largest-prime',type=int,default=113)
    parser.add_argument('--maximum-modulus',type=int,default=720)
    parser.add_argument('--per-shape',type=int,default=4)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    result=search(json.loads(Path(args.base).read_text()),args.largest_prime,args.maximum_modulus,args.per_shape)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['verification'],indent=2))
