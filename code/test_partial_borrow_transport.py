"""Check compact union evaluation and reject false transport certificates."""
import json
from pathlib import Path
import random

from construction_kernel import labels
from core_region_arrows import Core,row_union_contains
from partial_borrow_transport import validate


def main():
    core=Core(dict(minimum=3,primes=[2,3],depths=[4,3],rows=[],exact_uncovered_states=144))
    rng=random.Random(1311)
    comparisons=0
    for tail_mode in (False,True):
        levels=[2,2]
        for count in (0,1,2,4,8,12):
            rows=[]
            for _ in range(count):
                es=[rng.randrange(E+1) for E in levels]
                rs=[rng.randrange(labels(p,E,e) if tail_mode else p**e)
                    for p,E,e in zip(core.ps,levels,es)]
                rows.append(dict(exponents=es,residues=rs))
            node=dict(kind='row_union',levels=levels,tail_mode=tail_mode,rows=rows)
            for _ in range(50):
                es=[rng.randrange(E+1) for E in core.ds]
                rs=[rng.randrange(labels(p,E,e)) for p,E,e in zip(core.ps,core.ds,es)]
                pid=core.profile(es,rs)
                literal=row_union_contains(core,pid,node,False,False)
                assert literal==row_union_contains(core,pid,node,False,True)
                assert literal==row_union_contains(core,pid,node,True,True)
                comparisons+=1
    record=json.loads(Path('data/nielsen-prime13-compact.json').read_text())
    checked=validate(record)
    bad=json.loads(json.dumps(record))
    bad['direct_inputs'][0].append(bad['direct_inputs'][0][0].copy())
    try:
        validate(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('Repeated own resource accepted')
    bad=json.loads(json.dumps(record))
    bad['prime11_tail_inputs'][0][0][0]['residues'][0]=1
    try:
        validate(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('Incorrect transported prefix accepted')
    report=dict(status='PASS',cached_union_comparisons=comparisons,
                checked_transport_conditions=checked['checked_state_conditions'],
                repeated_resource_rejected=True,incorrect_prefix_rejected=True)
    Path('data/partial-borrow-transport-tests.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
