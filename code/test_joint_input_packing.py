"""Keep partial input packings distinct from complete prime-root claims."""
import json
from pathlib import Path
from joint_root_repair import search,validate,validate_packing,validate_capacity


def main():
    core=dict(minimum=3,primes=[2,3],depths=[3,2],rows=[],exact_uncovered_states=24)
    region=dict(exponents=[2,0],residues=[1,0])
    result=search(core,region,[2,1],5,5,True,3)
    assert result['status']=='VERIFIED_INPUT_PACKING'
    assert validate_packing(result)['input_count']==3
    try:
        validate(result)
    except AssertionError:
        pass
    else:
        raise AssertionError('Partial packing accepted as a whole prime root')
    excluded=search(core,region,[2,1],5,5,True,3,[(0,0),(1,0)],True)
    assert excluded['status']=='RESTRICTED_MENU_EXCLUDED'
    assert validate_capacity(excluded)['excluded']
    bad=json.loads(json.dumps(excluded))
    bad['allowed_resources'].append([99,0])
    try:
        validate_capacity(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('Invalid resource subset accepted')
    report=dict(status='PASS',verified_three_input_packing=True,
                incomplete_root_rejected=True,restricted_capacity_verified=True,
                invalid_resource_subset_rejected=True)
    Path('data/joint-input-packing-tests.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
