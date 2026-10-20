"""Exact witness refinement for joint root inputs with implicit row choices.

Every exponent resource chooses a color (root input) and residue labels.
Witness constraints use element relations rather than enumerating all
possible congruence classes. After each solution the entire projected
target is checked. Only a zero-gap independently verified result succeeds.
"""
import argparse
from itertools import product
import json
from math import prod
from pathlib import Path
import time

from joint_root_repair import target_points,validate,validate_packing
from graded_arrow_blocks import np
from construction_kernel import verify,prime
from ortools.sat.python import cp_model


def search(certificate,region,levels,q,seconds=120,input_count=None,seed_inputs=None,fixed_inputs=()):
    started=time.monotonic()
    verify(certificate)
    core,points,needed=target_points(certificate,region,levels,True)
    assert len(points) and prime(q) and q>max(core.ps)
    count=q-1 if input_count is None else input_count
    assert type(count) is int and 1<=count<=q-1
    fixed_inputs=set(fixed_inputs)
    assert not fixed_inputs or seed_inputs is not None
    if seed_inputs is not None:
        assert len(seed_inputs)==count
    assert fixed_inputs<=set(range(count))
    active=[i for i in range(count) if i not in fixed_inputs]
    assert active
    frozen={tuple(r['exponents']) for i in fixed_inputs for r in seed_inputs[i]}
    exponents=[es for es in product(*(range(E+1) for E in levels))
               if q*prod(p**e for p,e in zip(core.ps,es))>=certificate['minimum'] and es not in frozen]
    projections={(i,e):(points[:,i] if e==E else points[:,i]//(p-1)%p**e)
                 for i,(p,E) in enumerate(zip(core.ps,levels)) for e in range(E+1)}
    model=cp_model.CpModel()
    colors=[model.new_int_var(0,len(active)-1,f'color_{j}') for j in range(len(exponents))]
    residue_vars=[]
    primary=colors[:]
    for j,es in enumerate(exponents):
        row=[]
        for i,e in enumerate(es):
            allowed=list(map(int,np.unique(projections[i,e])))
            var=model.new_int_var_from_domain(cp_model.Domain.from_values(allowed),f'r_{j}_{i}')
            row.append(var)
            primary.append(var)
        residue_vars.append(row)
    # Canonical first occurrences of colors remove input permutations.
    if seed_inputs is None:
        model.add(colors[0]==0)
        high=colors[0]
        for j in range(1,len(colors)):
            model.add(colors[j]<=high+1)
            next_high=model.new_int_var(0,len(active)-1,f'high_{j}')
            model.add_max_equality(next_high,[high,colors[j]])
            high=next_high
    else:
        seed_rows={tuple(r['exponents']):(active.index(i),r['residues'])
                   for i in active for r in seed_inputs[i]}
        changes=[]
        for j,es in enumerate(exponents):
            if es in seed_rows:
                color,rs=seed_rows[es]
                changed=model.new_bool_var(f'changed_{j}')
                changes.append(changed)
                model.add(colors[j]==color).only_enforce_if(changed.Not())
                model.add_hint(colors[j],color)
                for var,r in zip(residue_vars[j],rs):
                    model.add_hint(var,r)
                    model.add(var==r).only_enforce_if(changed.Not())
        model.minimize(sum(changes))
    witnesses=set()
    def constrain(index):
        if index in witnesses:
            return
        witnesses.add(index)
        for d in range(len(active)):
            witness=model.new_int_var(0,len(exponents)-1,f'w_{index}_{d}')
            model.add_element(witness,colors,d)
            for i in range(len(core.ps)):
                constants=[int(projections[i,es[i]][index]) for es in exponents]
                value=model.new_int_var(0,max(constants),f'v_{index}_{d}_{i}')
                model.add_element(witness,[row[i] for row in residue_vars],value)
                model.add_element(witness,constants,value)
    for i in sorted({0,len(points)//2,len(points)-1}):
        constrain(i)
    if seed_inputs is not None:
        # Start from every failure of the published pattern, not arbitrary
        # target samples. Subsequent checks still exhaust the whole target.
        for d in active:
            covered=np.zeros(len(points),dtype=bool)
            for row in seed_inputs[d]:
                hit=np.ones(len(points),dtype=bool)
                for i,(e,r) in enumerate(zip(row['exponents'],row['residues'])):
                    hit&=projections[i,e]==r
                covered|=hit
            for index in np.flatnonzero(~covered):
                constrain(int(index))
    trace=[]
    outcome='TIME_LIMIT'
    last=None
    while time.monotonic()-started<seconds:
        solver=cp_model.CpSolver()
        solver.parameters.max_time_in_seconds=min(30,max(0.1,seconds-(time.monotonic()-started)))
        solver.parameters.num_search_workers=4
        status=solver.solve(model)
        info=dict(iteration=len(trace)+1,witnesses=len(witnesses),solver_status=solver.status_name(status),
                  elapsed_seconds=round(time.monotonic()-started,3))
        if seed_inputs is not None and status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
            info['changed_resources']=round(solver.objective_value)
        if status not in (cp_model.OPTIMAL,cp_model.FEASIBLE):
            trace.append(info)
            print(json.dumps(info),flush=True)
            outcome=solver.status_name(status)
            break
        inputs=[list(seed_inputs[i]) if i in fixed_inputs else [] for i in range(count)]
        for j,es in enumerate(exponents):
            inputs[active[solver.value(colors[j])]].append(dict(exponents=list(es),
                residues=[solver.value(var) for var in residue_vars[j]]))
        gaps=[]
        new=[]
        for rows in inputs:
            covered=np.zeros(len(points),dtype=bool)
            for row in rows:
                hit=np.ones(len(points),dtype=bool)
                for i,(e,r) in enumerate(zip(row['exponents'],row['residues'])):
                    hit&=projections[i,e]==r
                covered|=hit
            missed=np.flatnonzero(~covered)
            gaps.append(len(missed))
            if len(missed):
                assert int(missed[0]) not in witnesses
                new.append(int(missed[0]))
        info['exact_input_gaps']=gaps
        trace.append(info)
        print(json.dumps(info),flush=True)
        last=inputs
        if not new:
            outcome='CANDIDATE_COMPLETE'
            break
        model.clear_hints()
        for var in primary:
            model.add_hint(var,solver.value(var))
        for index in sorted(set(new)):
            constrain(index)
    result=dict(core=certificate,region=region,levels=levels,root_prime=q,tail_mode=True,
                required_inputs=count,frozen_input_indices=sorted(fixed_inputs),
                method='implicit_joint_root_witness_refinement',status='NO_VERIFIED_ROOT',
                search_outcome=outcome,projected_points=len(points),core_states_targeted=needed,
                resources=len(exponents),witnesses=len(witnesses),trace=trace,
                elapsed_seconds=time.monotonic()-started,
                scope='Witnesses guide search; every projected point must pass before any relative covering claim.')
    if outcome=='CANDIDATE_COMPLETE':
        result['inputs']=last
        result['verification']=validate_packing(result)
        result['status']=result['verification']['status']
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('core')
    parser.add_argument('region')
    parser.add_argument('--levels',default='5,3,2,2')
    parser.add_argument('--root',type=int,default=11)
    parser.add_argument('--seconds',type=float,default=120)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    result=search(json.loads(Path(args.core).read_text()),json.loads(Path(args.region).read_text()),
                  list(map(int,args.levels.split(','))),args.root,args.seconds)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('core','inputs','trace')},indent=2))
