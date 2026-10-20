"""Finite one-prime cover frontiers for resource-disjoint input supplies.

Unlike a single tail, a frontier may unite several fixed prime-power
classes and final tail digits. Count vectors allow disjoint packing using
T copies of each exponent resource. Packing here is a constructive greedy
heuristic; coverage and resource capacities are checked exactly.
"""
from functools import lru_cache
from itertools import product
import json
from pathlib import Path


@lru_cache(maxsize=2048)
def frontiers(p,E,states):
    states=tuple(sorted(set(states)))
    assert states and all(0<=s<p**(E-1)*(p-1) for s in states)
    zero=(0,)*(E+1)
    def visit(e,a,needed):
        if not needed:
            return {zero:()}
        result={}
        if e:
            counts=list(zero);counts[e]=1
            result[tuple(counts)]=(('fixed',e,a,0),)
        if e==E-1:
            leafs=tuple(('tail',E,a,s%(p-1)+1) for s in needed)
            counts=list(zero);counts[E]=len(leafs)
            result[tuple(counts)]=leafs
        else:
            combined={zero:()}
            for d in range(p):
                child=a+d*p**e
                child_states=tuple(s for s in needed if (s//(p-1))%p**(e+1)==child)
                options=visit(e+1,child,child_states)
                new={}
                for c1,f1 in combined.items():
                    for c2,f2 in options.items():
                        counts=tuple(x+y for x,y in zip(c1,c2))
                        if counts not in new:
                            new[counts]=f1+f2
                combined=new
            result.update(combined)
        return result
    return tuple(visit(0,0,states).items())


@lru_cache(maxsize=8192)
def pack(p,E,states,T):
    options=frontiers(p,E,tuple(states))
    best=[]
    for policy in range(3):
        remaining=[T]*(E+1)
        chosen=[]
        while True:
            feasible=[(c,f) for c,f in options if all(x<=y for x,y in zip(c,remaining))]
            if not feasible:
                break
            if policy==0:
                key=lambda cf:(sum(c/(r+1) for c,r in zip(cf[0],remaining)),sum(cf[0]),cf[0])
            elif policy==1:
                key=lambda cf:(max(c/(r+1) for c,r in zip(cf[0],remaining)),sum(cf[0]),cf[0])
            else:
                key=lambda cf:(sum(c*c/(r+1) for c,r in zip(cf[0],remaining)),sum(cf[0]),cf[0])
            counts,forest=min(feasible,key=key)
            chosen.append((counts,forest))
            remaining=[r-c for r,c in zip(remaining,counts)]
        if len(chosen)>len(best):
            best=chosen
    # Index the old input sets separately at every exponent; this makes
    # each pair (exponent, old_input) occur at most once.
    used=[0]*(E+1)
    packed=[]
    for counts,forest in best:
        row=[]
        for kind,e,a,d in forest:
            row.append(dict(kind=kind,exponent=e,prefix=a,digit=d,old_input=used[e]))
            used[e]+=1
        packed.append(row)
    assert all(u<=T for u in used)
    return packed


def test():
    cases=0
    covers=0
    for p,E in [(2,3),(3,2),(5,1),(3,3)]:
        n=p**(E-1)*(p-1)
        targets=[tuple(range(n)),tuple(range(0,n,2)),(0,),tuple(range(n//2))]
        for states in targets:
            if not states:
                continue
            for T in [1,2,5,10]:
                packed=pack(p,E,states,T)
                resources=set()
                for forest in packed:
                    for branch in forest:
                        key=branch['exponent'],branch['old_input']
                        assert key not in resources
                        resources.add(key)
                    for s in states:
                        prefix,d0=divmod(s,p-1)
                        def hits(b):
                            if b['kind']=='fixed':
                                return prefix%p**b['exponent']==b['prefix']
                            return prefix==b['prefix'] and d0+1==b['digit']
                        assert sum(hits(b) for b in forest)==1
                    covers+=1
                cases+=1
    # The union 1 mod 3 OR 3 mod 9 has a frontier using one resource
    # of each exponent 1 and 2, unlike the single-tail count rule.
    p,E=3,3
    states=tuple(s for s in range(18) if (s//2)%3==1 or (s//2)%9==3)
    options=frontiers(p,E,states)
    assert (0,1,1,0) in dict(options)
    report=dict(status='PASS',packing_cases=cases,exact_disjoint_frontier_covers=covers,
                two_progression_union_profile=[0,1,1,0],
                scope='Packing primitive used by partial minimum-50 region schemas; this test alone does not certify a covering system.')
    Path('data/prime-forest-supply-tests.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    test()
