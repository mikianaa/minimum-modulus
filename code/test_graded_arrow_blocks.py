"""Exhaustive finite-kernel audits of the graded input-set construction."""
from itertools import product
import json
from math import prod
from pathlib import Path
import random

from construction_kernel import matches,project,sizes
from graded_arrow_blocks import supply,check,search


def run():
    rng=random.Random(5050)
    cases=0
    sets_checked=0
    for ps in [[2,3],[2,3,5],[3,5],[2,5,7]]:
        for trial in range(12):
            fixed={p:rng.randrange(2) for p in ps}
            ds=[fixed[p]+1 for p in ps]
            order=rng.sample(ps,len(ps))
            current,nodes=supply(order,fixed,True)
            grades,_=supply(order,fixed,False)
            assert [g for g,_ in current]==[g for g,_ in grades]
            a=rng.randrange(prod(p**fixed[p] for p in ps))
            def expand(idx):
                node=nodes[idx]
                if node['kind']=='unit':
                    return [dict(exponents=[0]*len(ps),residues=[0]*len(ps))]
                p,e=node['prime'],node['exponent']
                i=ps.index(p)
                out=[]
                for digit,child in enumerate(node['children']):
                    for row in expand(child):
                        assert row['exponents'][i]==0
                        row['exponents'][i]=e
                        row['residues'][i]=(a%p**e if node['kind']=='fixed'
                            else (a%p**(e-1))*(p-1)+digit)
                        out.append(row)
                return out
            all_resources=set()
            target=[x for x in product(*(range(d) for d in sizes(ps,ds)))
                    if all(project(v,p,E,fixed[p])==a%p**fixed[p]
                           for v,p,E in zip(x,ps,ds))]
            for grade,idx in current:
                rows=expand(idx)
                for row in rows:
                    key=tuple(row['exponents'])
                    assert key not in all_resources
                    all_resources.add(key)
                assert min(prod(p**e for p,e in zip(ps,r['exponents'])) for r in rows)==grade
                assert all(any(matches(x,r,ps,ds) for r in rows) for x in target)
                sets_checked+=1
            cases+=1
    control=search(3,13,[3,2,1,0],2)
    assert control['status']=='VERIFIED_CONSTRUCTION_SCHEMA'
    Path('data/construction-graded-control-3.json').write_text(json.dumps(control,indent=2)+'\n')
    bad=json.loads(json.dumps(control))
    occupied=next(b for b in bad['blocks'] if b is not None)
    bad['blocks'].append(occupied)
    try:
        check(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('Repeated root exponent family was accepted')
    result=dict(status='PASS',exhaustive_supply_cases=cases,
                resource_disjoint_target_covers_checked=sets_checked,
                minimum_3_positive_control='VERIFIED_CONSTRUCTION_SCHEMA',
                duplicate_stage_rejected=True)
    Path('data/graded-arrow-block-tests.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    run()
