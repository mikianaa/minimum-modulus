"""Audit Nielsen's partial-borrowing doubling step at prime 13.

Old input coverage and new own rows are checked jointly for each prime-11
digit. Borrowed old rows are not copied into the new modulus families.
"""
from math import prod
import json
from pathlib import Path

from graded_arrow_blocks import np
from joint_root_repair import target_points
from nielsen_prime11 import inputs,expand,seed_certificate,PS,DS


def coverage(points,rows,depths=DS):
    covered=np.zeros(len(points),dtype=bool)
    for row in rows:
        hit=np.ones(len(points),dtype=bool)
        for i,(p,E,e,r) in enumerate(zip(PS,depths,row['exponents'],row['residues'])):
            hit&=(points[:,i]==r if e==E else points[:,i]//(p-1)%p**e==r)
        covered|=hit
    return covered


def own_copies(depths=DS):
    families=inputs(binary_residue=3,raw=True)
    copies=[]
    for collapse in (False,True):
        branches=[]
        for item in families:
            own=[]
            for original in item:
                if 2 not in original or original[2][1]<2:
                    continue
                family=dict(original)
                if collapse:
                    kind,e,r,d=family.pop(2)
                    if kind=='tail':
                        assert e==3
                        family[2]=('fixed',1,1,0)
                    else:
                        assert e==2
                own.extend(expand(family,depths))
            branches.append(own)
        copies.append(branches)
    return copies


def audit(minimum=42,depths=DS):
    core=seed_certificate(minimum)
    region=dict(exponents=[2,2,0,0],residues=[3,0,0,0],unions={'1':[1,3,4,7]})
    _,points,needed=target_points(core,region,depths,True)
    old=inputs(depths=depths)
    direct=inputs(binary_residue=3,permutation={1:2,2:1},depths=depths)
    copies=own_copies(depths)
    seen=set()
    checks=[]
    minimum_grade=None
    def resources(rows,eleven):
        nonlocal minimum_grade
        for row in rows:
            key=tuple(row['exponents'])+(eleven,)
            assert key not in seen,('Repeated own exponent resource',key)
            seen.add(key)
            grade=13*11**eleven*prod(p**e for p,e in zip(PS,row['exponents']))
            minimum_grade=grade if minimum_grade is None else min(minimum_grade,grade)
            assert grade>=minimum
    for i,rows in enumerate(direct,1):
        resources(rows,0)
        gap=int((~coverage(points,rows,depths)).sum())
        checks.append(dict(input=i,prime11_digit=None,missing=gap))
    old_coverage=[coverage(points,rows,depths) for rows in old]
    for i,branches in enumerate(copies,11):
        for d,rows in enumerate(branches,1):
            resources(rows,1)
            own=coverage(points,rows,depths)
            borrowed=old_coverage[d-1]
            gap=int((~(own|borrowed)).sum())
            checks.append(dict(input=i,prime11_digit=d,missing=gap,
                               states_requiring_borrowed_coverage=int((~own&borrowed).sum())))
    complete=all(c['missing']==0 for c in checks)
    result=dict(core=core,region=region,root_prime=13,levels=depths,
                direct_inputs=direct,prime11_tail_inputs=copies,old_prime11_inputs=old,
                status='VERIFIED_RELATIVE_PRIME13' if complete else 'UNVERIFIED_TRANSCRIPTION',
                audit=dict(target_states=needed,own_resources=len(seen),minimum_own_modulus=minimum_grade,
                           checked_relative_conditions=len(checks),conditions=checks),
                scope='Prime-13 target only, relative to the early core and recorded prime-11 stage; not a full covering system.')
    return result


if __name__=='__main__':
    result=audit()
    Path('data/nielsen-prime13-transcription.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],audit=result['audit']),indent=2))
