"""Retain the valid higher powers of the forbidden modulus-44 family."""
import json
from pathlib import Path

from graded_arrow_blocks import np
from construction_kernel import verify
from joint_root_repair import target_points,validate_packing


def build():
    packed=json.loads(Path('data/construction-50-relocated90-nine-inputs.json').read_text())
    validate_packing(packed)
    ps=packed['core']['primes']
    levels=packed['levels']
    q=packed['root_prime']
    assert q==11 and len(packed['inputs'])==9
    rows=[]
    for d,classes in enumerate(packed['inputs'],2):
        for row in classes:
            rows.append(dict(exponents=row['exponents']+[1],residues=row['residues']+[d]))
            rows.append(dict(exponents=row['exponents']+[2],residues=row['residues']+[d-1]))
    # The old input 1 used 4 times the whole 11-tail. Delete only its first
    # class (modulus 44), retaining 4 times the 121-tail at digit 1.
    rows.append(dict(exponents=[2,0,0,0,2],residues=[1,0,0,0,0]))
    own=dict(minimum=50,primes=ps+[q],depths=levels+[2],rows=rows)
    verify(own)  # Validates all own resource distinctions and minima.
    _,points,needed=target_points(packed['core'],packed['region'],levels,True)
    root_states=np.arange(q*(q-1))
    root_values=root_states//(q-1)+(root_states%(q-1)+1)*q
    covered=np.zeros((len(points),len(root_states)),dtype=bool)
    for row in rows:
        body=np.ones(len(points),dtype=bool)
        for i,(p,E,e,r) in enumerate(zip(ps,levels,row['exponents'],row['residues'])):
            body&=(points[:,i]==r if e==E else points[:,i]//(p-1)%p**e==r)
        e,r=row['exponents'][-1],row['residues'][-1]
        qhit=(root_states==r if e==2 else root_values%q==r)
        covered|=body[:,None]&qhit[None,:]
    exception=root_values%q==1
    assert covered[:,~exception].all() and not covered[:,exception].any()
    record=dict(core=packed['core'],target_region=packed['region'],own_certificate=own,
                status='VERIFIED_RELATIVE_ROOT_EXCEPT_MODULUS44',
                exception=dict(modulus=44,residue=1,
                               scope='Only the stated target, after the early core, is being audited.'),
                verification=dict(own_resources=len(rows),projected_core_states=len(points),
                    root_states=len(root_states),covered_input_states=int(covered.sum()),
                    remaining_input_states=int((~covered).sum()),core_target_states=needed),
                finite_completion='Use a different auxiliary prime above all structural primes for each own terminal-tail occurrence. The checked exponent families are disjoint.',
                scope='This preserves the higher-power part of the first input; it is not a complete covering system.')
    Path('data/construction-50-relocated90-prime11-truncated.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status=record['status'],verification=record['verification'],exception=record['exception']),indent=2))
    return record


if __name__=='__main__':
    build()
