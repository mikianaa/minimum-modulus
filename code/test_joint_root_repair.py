"""Positive and deliberately corrupted checks for joint input packing."""
import json
from pathlib import Path
from joint_root_repair import search,validate,validate_capacity


def main():
    certificate=dict(minimum=3,primes=[2,3],depths=[3,2],rows=[],exact_uncovered_states=24)
    region=dict(exponents=[2,0],residues=[1,0])
    result=search(certificate,region,[2,1],5,10)
    assert result['status']=='VERIFIED_RELATIVE_ROOT'
    checked=0
    # Direct congruence evaluation over an entire integer period; this does
    # not use the kernel projection or the validator's membership arrays.
    for rows in result['inputs']:
        for x in range(12):
            if x%4==1:
                assert any(all(x%p**e==r for p,e,r in zip([2,3],row['exponents'],row['residues']))
                           for row in rows)
                checked+=1
    bad=json.loads(json.dumps(result))
    bad['inputs'][1].append(bad['inputs'][0][0].copy())
    try:
        validate(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('Duplicated resource accepted')
    bad=json.loads(json.dumps(result))
    bad['inputs'][0]=[]
    try:
        validate(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('Uncovered input accepted')
    excluded=search(certificate,region,[2,1],7,10)
    assert excluded['status']=='RESTRICTED_MENU_EXCLUDED'
    assert validate_capacity(excluded)['excluded']
    excluded['exact_weighted_capacity']['numerator']-=1
    try:
        validate_capacity(excluded)
    except AssertionError:
        pass
    else:
        raise AssertionError('Incorrect capacity numerator accepted')
    report=dict(status='PASS',direct_integer_input_checks=checked,
                duplicate_resource_rejected=True,uncovered_input_rejected=True,
                exact_capacity_control=True,incorrect_capacity_rejected=True)
    Path('data/joint-root-repair-tests.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
