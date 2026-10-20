"""Move modulus 90 to repair the odd part of the deleted modulus-45 class.

The replacement 3 mod 90 restores all missing core-45 coverage needed by
the published prime-11 and prime-13 targets, both of which are odd. Its
former even class becomes a new repair obligation, explicitly retained.
"""
from copy import deepcopy
import json
from math import prod
from pathlib import Path

from graded_arrow_blocks import np
from construction_kernel import sizes,verify
from nielsen_prime11 import inputs,PS,DS
from joint_root_repair import validate_packing
from core_region_arrows import Core


def build():
    old=json.loads(Path('data/nielsen-prime11-minimum50-diagnostic.json').read_text())['core']
    core=deepcopy(old)
    rows=[r for r in core['rows'] if r['exponents']==[1,2,1,0]]
    assert len(rows)==1
    before=deepcopy(rows[0])
    rows[0]['residues']=[1,3,3,0]
    states=np.indices(sizes(PS,DS),sparse=True)
    residual=np.ones(sizes(PS,DS),dtype=bool)
    for row in core['rows']:
        hit=True
        for s,p,E,e,r in zip(states,PS,DS,row['exponents'],row['residues']):
            hit=hit & (s==r if e==E else s//(p-1)%p**e==r)
        residual&=~hit
    core['exact_uncovered_states']=int(residual.sum())
    core['relocation']=dict(modulus=90,before=before,after=deepcopy(rows[0]),
                            new_class='3 mod 90',lost_class='78 mod 90',
                            scope='The lost even class remains an explicit repair obligation.')
    core['verification']=verify(core)
    core['status']=core['verification']['status']
    Path('data/construction-50-relocated90-core.json').write_text(json.dumps(core,indent=2)+'\n')
    old42=json.loads(Path('data/nielsen-prime11-transcription.json').read_text())['core']
    before42=Core(old42)
    after50=Core(core)
    odd_index=before42.index([1,2,0,0],[1,0,0,0],{'1':[1,3,4,7]})
    assert np.array_equal(before42.residual[odd_index],after50.residual[odd_index])
    agreement=dict(status='VERIFIED_EXACT_ODD_TARGET_RESTORATION',
                   target='odd integers intersect (1 mod 3 union 3 mod 9)',
                   checked_kernel_states=int(before42.residual[odd_index].size),
                   identical_uncovered_states=int(before42.residual[odd_index].sum()),
                   reason='3 mod 90 restores the odd half of 3 mod 45; the lost 78 mod 90 class is even. The deleted 42 and 48 rows do not meet this target.',
                   scope='Equality on this target only; other new holes remain.')
    Path('data/construction-50-relocated90-odd-agreement.json').write_text(json.dumps(agreement,indent=2)+'\n')
    levels=[3,4,2,1]
    result=dict(core=core,region=dict(exponents=[2,2,0,0],residues=[1,0,0,0],unions={'1':[1,3,4,7]}),
                levels=levels,tail_mode=True,root_prime=11,required_inputs=9,
                inputs=inputs(depths=levels)[1:],method='published_nine_inputs_after_modulus90_relocation')
    result['verification']=validate_packing(result)
    result['status']=result['verification']['status']
    result['scope']='Nine complete admissible inputs; the tenth input still has its modulus-44 gap. Not a full covering system.'
    Path('data/construction-50-relocated90-nine-inputs.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(core_uncovered_states=core['exact_uncovered_states'],
                         packing=result['verification'],relocation=core['relocation']),indent=2))
    return result


if __name__=='__main__':
    build()
