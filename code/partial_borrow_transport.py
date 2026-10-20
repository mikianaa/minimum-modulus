"""Checked source DAG for the published partial-borrowing prime-13 step."""
from math import prod
from construction_kernel import verify,labels
from joint_root_repair import target_points
from nielsen_prime13 import coverage,PS,DS


def validate(record):
    core=record['core']
    verify(core)
    assert core['primes']==PS and core['depths']==DS
    levels=record['levels']
    assert record['root_prime']==13 and len(levels)==4 and all(1<=e<=E for e,E in zip(levels,DS))
    assert len(record['direct_inputs'])==len(record['old_prime11_inputs'])==10
    assert len(record['prime11_tail_inputs'])==2
    _,points,needed=target_points(core,record['region'],levels,True)
    seen=set()
    minimum=None
    checks=0
    def check_rows(rows,eleven):
        nonlocal minimum
        for row in rows:
            es,rs=row['exponents'],row['residues']
            assert len(es)==len(rs)==4
            assert all(type(e) is int and 0<=e<=E for e,E in zip(es,levels))
            assert all(type(r) is int and 0<=r<labels(p,E,e) for r,p,E,e in zip(rs,PS,levels,es))
            key=tuple(es)+(eleven,)
            assert key not in seen
            seen.add(key)
            n=13*11**eleven*prod(p**e for p,e in zip(PS,es))
            assert n>=core['minimum']
            minimum=n if minimum is None else min(n,minimum)
    old=[coverage(points,rows,levels) for rows in record['old_prime11_inputs']]
    for rows in record['direct_inputs']:
        check_rows(rows,0)
        assert coverage(points,rows,levels).all()
        checks+=len(points)
    for branches in record['prime11_tail_inputs']:
        assert len(branches)==10
        for d,rows in enumerate(branches):
            check_rows(rows,1)
            assert (coverage(points,rows,levels)|old[d]).all()
            checks+=len(points)
    return dict(status='VERIFIED_PARTIAL_BORROWING_STAGE',root_prime=13,own_resources=len(seen),
                minimum_own_modulus=minimum,checked_state_conditions=checks,target_states=needed,
                scope='Relative prime-13 stage; the older prime-11 rows must be authenticated separately.')


def source(record):
    def grade(rows):
        return min((prod(p**e for p,e in zip(PS,row['exponents'])) for row in rows),default=None)
    levels=record['levels']
    nodes=[dict(kind='row_union',levels=levels,tail_mode=True,rows=rows,grade=grade(rows))
           for rows in record['direct_inputs']]
    selected=list(range(10))
    for branches in record['prime11_tail_inputs']:
        children=[]
        for d,rows in enumerate(branches,1):
            children.append(len(nodes))
            nodes.append(dict(kind='relative_core_union',levels=levels,tail_mode=True,rows=rows,
                              reference_prime=11,reference_input=d-1,grade=grade(rows)))
        selected.append(len(nodes))
        nodes.append(dict(kind='tail',prime=11,exponent=1,prefix=0,digits=list(range(1,11)),
                          borrowed=[],children=children,
                          grade=11*min(nodes[i]['grade'] for i in children if nodes[i]['grade'] is not None)))
    return dict(root_prime=13,region=record['region'],nodes=nodes,root_inputs=selected,
                partial_borrowing_certificate=record)
