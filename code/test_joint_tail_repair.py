"""Independent exhaustive checks for projecting and borrowing tail inputs."""
from itertools import product
import json
from pathlib import Path

from construction_kernel import matches
from core_region_arrows import Core,Catalog,joint_source,search as continue_search
from joint_root_repair import target_points,search,validate


def local_state(x,p,E,global_E):
    # Congruence enumeration, independent of the construction's division loop.
    for a in range(p**(E-1)):
        for d in range(1,p):
            if any((x-a-d*p**(k-1))%p**k==0 for k in range(E,global_E+1)):
                return a*(p-1)+d-1
    raise AssertionError('No lower tail state')


def main():
    cert=dict(minimum=3,primes=[2,3],depths=[4,3],rows=[],exact_uncovered_states=144)
    region=dict(exponents=[2,1],residues=[1,0],unions={'1':[1,2]})
    levels=[2,2]
    core,points,needed=target_points(cert,region,levels,True)
    expected=set()
    full=[]
    for s in product(*(range(n) for n in core.dims)):
        if s[0]%4==1 and (s[1]//2)%3 in (1,2):
            x=[int(core.values[i][v]) for i,v in enumerate(s)]
            expected.add(tuple(local_state(v,p,E,G) for v,p,E,G in zip(x,core.ps,levels,core.ds)))
            full.append(s)
    assert expected=={tuple(map(int,p)) for p in points} and len(full)==needed
    simple=dict(exponents=[2,0],residues=[1,0])
    result=search(cert,simple,[2,1],5,10,True)
    assert result['status']=='VERIFIED_RELATIVE_ROOT'
    catalog=Catalog(core,True)
    catalog.append(joint_source(result))
    implications=0
    for es in product(range(5),range(4)):
        from construction_kernel import labels
        for rs in product(*(range(labels(p,E,e)) for p,E,e in zip(core.ps,core.ds,es))):
            pid=core.profile(es,rs)
            selected=[s for s in product(*(range(n) for n in core.dims))
                      if matches(s,dict(exponents=es,residues=rs),core.ps,core.ds)]
            for j,rows in enumerate(result['inputs']):
                if catalog.contains(0,j,pid):
                    for s in selected:
                        x=[int(core.values[i][v]) for i,v in enumerate(s)]
                        point=tuple(local_state(v,p,E,G) for v,p,E,G in zip(x,core.ps,[2,1],core.ds))
                        assert any(matches(point,r,core.ps,[2,1]) for r in rows)
                        implications+=1
    continuation=continue_search(cert,largest=11,max_volume=144,per_shape=4,seconds=10,
                                 branch_sensitive=True,initial_joint=result)
    assert continuation['verification']['checked_joint_roots']==1
    validate(result)
    report=dict(status='PASS',exhaustive_projection_states=144,projected_target_states=needed,
                checked_borrowing_implications=implications,
                integrated_joint_roots=1,continuation_status=continuation['status'])
    Path('data/joint-tail-repair-tests.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
