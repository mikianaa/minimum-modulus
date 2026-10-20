"""Probe branch-specific coverage of previously constructed arrow blocks.

For each target progression, count prior root-prime inputs that already
cover it, using the old explicit input DAGs. These can be omitted from a
new relative arrow: only missing inputs consume new modulus families.
This script records sufficient input counts, not a full covering certificate.
"""
from functools import lru_cache
from itertools import product
import argparse
import json
from math import prod
from pathlib import Path

from graded_arrow_blocks import supply,check,np
from owens_branch_search import primes_through


class OldBlocks:
    def __init__(self,base):
        check(base)
        self.ps=base['core_primes']
        self.sources={}
        for block in base['blocks']:
            if block is None or not block['tail'] or block['exponent']!=1:
                continue
            p=block['root_prime']
            if p in self.ps:
                continue
            fixed={int(q):e for q,e in block['fixed_exponents'].items()}
            current,nodes=supply(block['order'],fixed,True)
            good=[idx for grade,idx in current if grade>=block['input_threshold']]
            groups=[]
            for k,a in enumerate(block['residues']):
                groups.append((a,good[k*(p-1):(k+1)*(p-1)]))
            self.sources[p]=(nodes,groups)

    @lru_cache(None)
    def covered_digits(self,p,a,exponents):
        if p not in self.sources:
            return ()
        nodes,groups=self.sources[p]
        fixed=dict(zip(self.ps,exponents))

        @lru_cache(None)
        def contains(idx,old_a):
            node=nodes[idx]
            if node['kind']=='unit':
                return True
            q,e=node['prime'],node['exponent']
            f=fixed.get(q,0)
            if node['kind']=='fixed':
                return f>=e and (a-old_a)%q**e==0 and contains(node['children'][0],old_a)
            # The target must be inside the arrow's fixed prefix.
            prefix=old_a%q**(e-1)
            if f<e-1 or (a-prefix)%q**(e-1):
                return False
            delta=(a-prefix)%q**f
            if delta==0:
                return all(contains(child,old_a) for child in node['children'])
            while delta%q==0:
                delta//=q
            return contains(node['children'][delta%q-1],old_a)

        return tuple(d for d in range(1,p)
                     if any(contains(children[d-1],old_a) for old_a,children in groups))


def count_relative(order,fixed,missing):
    count=1
    trace=[]
    for p in order:
        k=missing.get(p,p-1)
        if k==0:
            return None,trace
        nxt=(fixed.get(p,0)+1)*count+count//k
        trace.append(dict(prime=p,required_inputs=k,before=count,after=nxt))
        count=nxt
    return count,trace


def probe(base,root=97,maximum_modulus=720,per_shape=4):
    assert root>base['largest_prime'] and root in primes_through(root)
    prior=OldBlocks(base)
    L=base['period']
    covered=np.zeros(L,dtype=bool)
    for block in base['blocks']:
        if block is not None:
            for a in block['residues']:
                covered[a::block['modulus']]=True
    target_shapes=[]
    ps=base['core_primes']
    for es in product(*(range(e+1) for e in base['core_exponents'])):
        m=prod(p**e for p,e in zip(ps,es))
        if m<=maximum_modulus:
            target_shapes.append((m,es))
    earlier=primes_through(root-1)
    attempts=[]
    for m,es in sorted(target_shapes):
        hist=(~covered).reshape(L//m,m).sum(axis=0,dtype=np.int64)
        for aa in np.argsort(-hist,kind='stable')[:per_shape]:
            a=int(aa)
            if not hist[a]:
                continue
            fixed={p:e for p,e in zip(ps,es) if e}
            borrowed={p:prior.covered_digits(p,a,es) for p in prior.sources}
            missing={p:p-1-len(ds) for p,ds in borrowed.items()}
            natural=earlier
            special=[p for p in earlier if p in fixed]
            ordinary=sorted((p for p in earlier if p not in fixed),key=lambda p:(missing.get(p,p-1),p))
            orders=[natural,special+ordinary,list(reversed(special))+ordinary]
            scored=[(count_relative(order,fixed,missing),order) for order in orders]
            best,order=max(scored,key=lambda t:float('inf') if t[0][0] is None else t[0][0])
            count,trace=best
            attempts.append(dict(residue=a,modulus=m,fixed_exponents=list(es),
                                 projected_gain=int(hist[a]),available_sets=count,
                                 required_sets=root-1,order=order,trace=trace,
                                 borrowed_digits={str(p):list(ds) for p,ds in borrowed.items() if ds},
                                 sufficient=count is None or count>=root-1,
                                 already_covered_by_prior_block=count is None))
    successes=[a for a in attempts if a['sufficient']]
    successes.sort(key=lambda a:(-a['projected_gain'],a['modulus'],a['residue']))
    return dict(status='SUFFICIENT_RELATIVE_INPUT_COUNTS_FOUND' if successes else 'NO_SUFFICIENT_COUNT_FOUND',
                root_prime=root,maximum_target_modulus=maximum_modulus,targets_per_shape=per_shape,
                targets_tested=len(attempts),sufficient_targets=len(successes),
                best_targets=successes[:20],attempts=attempts,
                scope='A sufficient-count probe on frozen old blocks, not a full construction certificate. Relative-arrow schema export and an independent dependency audit are still required.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base')
    parser.add_argument('--root',type=int,default=97)
    parser.add_argument('--maximum-modulus',type=int,default=720)
    parser.add_argument('--per-shape',type=int,default=4)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    result=probe(json.loads(Path(args.base).read_text()),args.root,args.maximum_modulus,args.per_shape)
    result['base_file']=args.base
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='attempts'},indent=2))
