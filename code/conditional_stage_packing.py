"""Pack a new root against the actual joint coverage of earlier stages.

The prime-11 exceptional head and the prime-13 partial contributions are
evaluated together before projecting. A projection is existential on the
uncovered set: one surviving fine state keeps its entire coarse cell as an
obligation. Thus no incidental coverage is assumed from a sampled point.
"""
import argparse
from itertools import product
import json
from math import prod
from pathlib import Path
import time

from graded_arrow_blocks import np
from construction_kernel import verify, labels, prime
from joint_root_repair import target_points, validate_packing
from nielsen_prime13 import coverage
from partial_borrow_transport import validate as validate_transport
from ortools.sat.python import cp_model


def read(path):
    return json.loads(Path(path).read_text())


def source_record():
    packed=read('data/construction-50-relocated90-nine-inputs.json')
    transport=read('data/construction-50-relocated90-prime13.json')
    validate_packing(packed)
    validate_transport(transport)
    assert packed['core']==transport['core']
    assert transport['old_prime11_inputs']==[[]]+packed['inputs']
    return dict(core=packed['core'],packing=packed,transport=transport)


def prior_rows(record):
    """Actual own rows, with no borrowed row copied into a new family."""
    packed,transport=record['packing'],record['transport']
    assert packed['core']==record['core']==transport['core']
    validate_packing(packed)
    validate_transport(transport)
    assert packed['levels']==transport['levels']
    assert packed['root_prime']==11 and packed['required_inputs']==9
    assert transport['old_prime11_inputs']==[[]]+packed['inputs']
    rows=[]
    def add(body,eleven,thirteen):
        for b in body:
            for e11,r11 in eleven:
                rows.append(dict(exponents=b['exponents']+[e11,thirteen[0]],
                                 residues=b['residues']+[r11,thirteen[1]]))
    for d,body in enumerate(packed['inputs'],2):
        add(body,[(1,d),(2,d-1)],(0,0))
    # Keep all powers >= 2 of the old 4 times 11-tail, dropping only 44.
    add([dict(exponents=[2,0,0,0],residues=[1,0,0,0])],[(2,0)],(0,0))
    for d,body in enumerate(transport['direct_inputs']):
        add(body,[(0,0)],(1,d))
    for d,branches in enumerate(transport['prime11_tail_inputs'],10):
        for digit,body in enumerate(branches,1):
            add(body,[(1,digit),(2,digit-1)],(1,d))
    certificate=dict(minimum=record['core']['minimum'],primes=[2,3,5,7,11,13],
                     depths=packed['levels']+[2,1],rows=rows)
    verify(certificate)  # Includes all resource and minimum checks.
    return certificate


def projected_target(record,region,independent=False):
    """Project joint fine obligations to core tails x (Z/11Z) x 13-tail.

    The core projection is conservative, since early core rows may vary
    inside a fiber; prior rows are constant there. Prime 11 is then reduced
    to ordinary residues only AFTER evaluating every one of its 110 states.
    """
    prior=prior_rows(record)
    levels=record['packing']['levels']
    _,body,fine_needed=target_points(record['core'],region,levels,True)
    survivors=np.ones((len(body),110,12),dtype=bool)
    for row in prior['rows']:
        if independent:
            hit=np.ones(len(body),dtype=bool)
            for i,(p,E,e,r) in enumerate(zip([2,3,5,7],levels,row['exponents'],row['residues'])):
                values=body[:,i]//(p-1)+(body[:,i]%(p-1)+1)*p**(E-1)
                residue=r//(p-1)+(r%(p-1)+1)*p**(E-1) if e==E else r
                hit&=(values-residue)%p**e==0
        else:
            hit=coverage(body,[row],levels)
        if not hit.any():
            continue
        e11,e13=row['exponents'][-2:]
        r11,r13=row['residues'][-2:]
        states=np.arange(110)
        qhit=(states==r11 if e11==2 else states//10%11**e11==r11)
        dhit=np.arange(12)==r13 if e13 else np.ones(12,dtype=bool)
        survivors[np.ix_(np.flatnonzero(hit),np.flatnonzero(qhit),np.flatnonzero(dhit))]=False
    if independent:
        projected=np.stack([survivors[:,list(range(10*a,10*a+10)),:].any(axis=1)
                            for a in range(11)],axis=1)
    else:
        projected=survivors.reshape(len(body),11,10,12).any(axis=2)
    indices=np.argwhere(projected)
    points=np.column_stack([body[indices[:,0]],indices[:,1:]])
    metadata=dict(core_projected_states=len(body),fine_core_target_states=fine_needed,
                  joint_kernel_obligations=int(survivors.sum()),
                  projected_obligations=len(points),prior_own_resources=len(prior['rows']))
    return points,metadata


def projections(points,ps,levels,modes):
    return {(i,e):(points[:,i] if mode=='tail' and e==E else
                  (points[:,i]//(p-1) if mode=='tail' else points[:,i])%p**e)
            for i,(p,E,mode) in enumerate(zip(ps,levels,modes)) for e in range(E+1)}


def row_hit(row,proj,n):
    hit=np.ones(n,dtype=bool)
    for i,(e,r) in enumerate(zip(row['exponents'],row['residues'])):
        if e:
            hit&=proj[i,e]==r
    return hit


def validate(result):
    record=result['prior']
    assert result['primes']==[2,3,5,7,11,13]
    assert result['levels']==record['packing']['levels']+[1,1]
    assert result['modes']==['tail']*4+['fixed','tail']
    points,metadata=projected_target(record,result['region'],True)
    ps,levels,modes=result['primes'],result['levels'],result['modes']
    q=result['root_prime']
    assert prime(q) and q>max(ps)
    assert len(result['inputs'])==q-1
    seen=set()
    gaps=[]
    # The arithmetic congruence test differs from the search projections.
    values=[points[:,i]//(p-1)+(points[:,i]%(p-1)+1)*p**(E-1)
            if mode=='tail' else points[:,i]
            for i,(p,E,mode) in enumerate(zip(ps,levels,modes))]
    for rows in result['inputs']:
        covered=np.zeros(len(points),dtype=bool)
        for row in rows:
            es,rs=row['exponents'],row['residues']
            assert len(es)==len(rs)==len(ps)
            assert all(type(e) is int and 0<=e<=E for e,E in zip(es,levels))
            assert all(type(r) is int and 0<=r<(labels(p,E,e) if mode=='tail' else p**e)
                       for p,E,mode,e,r in zip(ps,levels,modes,es,rs))
            assert tuple(es) not in seen
            seen.add(tuple(es))
            assert q*prod(p**e for p,e in zip(ps,es))>=record['core']['minimum']
            hit=np.ones(len(points),dtype=bool)
            for value,p,E,mode,e,r in zip(values,ps,levels,modes,es,rs):
                residue=r//(p-1)+(r%(p-1)+1)*p**(E-1) if mode=='tail' and e==E else r
                hit&=(value-residue)%p**e==0
            covered|=hit
        gaps.append(int((~covered).sum()))
    assert not any(gaps),gaps
    return dict(status='VERIFIED_CONDITIONAL_RELATIVE_ROOT',root_prime=q,
                own_resources=len(seen),checked_input_points=len(points)*(q-1),**metadata,
                scope='Only this target is covered relative to the authenticated earlier core and stages.')


def search(record,region,q=17,seconds=120,frozen_inputs=()):
    started=time.monotonic()
    points,metadata=projected_target(record,region)
    ps=[2,3,5,7,11,13]
    levels=record['packing']['levels']+[1,1]
    modes=['tail']*4+['fixed','tail']
    assert prime(q) and q>max(ps) and len(points)
    proj=projections(points,ps,levels,modes)
    frozen_resources=set()
    for rows in frozen_inputs:
        covered=np.zeros(len(points),dtype=bool)
        for row in rows:
            key=tuple(row['exponents'])
            assert key not in frozen_resources
            frozen_resources.add(key)
            assert q*prod(p**e for p,e in zip(ps,key))>=record['core']['minimum']
            covered|=row_hit(row,proj,len(points))
        assert covered.all(),int((~covered).sum())
    count=q-1-len(frozen_inputs)
    assert count>0
    exponents=[es for es in product(*(range(E+1) for E in levels))
               if q*prod(p**e for p,e in zip(ps,es))>=record['core']['minimum']
               and es not in frozen_resources]
    result=dict(prior=record,region=region,primes=ps,levels=levels,modes=modes,
                root_prime=q,method='joint_conditional_stage_packing',status='NO_VERIFIED_ROOT',
                resources=len(exponents),target_audit=metadata,frozen_inputs=list(frozen_inputs))
    print(json.dumps(dict(event='target_built',**metadata,resources=len(exponents))),flush=True)
    # An exact uniform capacity screen over every projected obligation.
    maxima=[]
    for es in exponents:
        code=np.zeros(len(points),dtype=np.int64)
        for i,(p,E,mode,e) in enumerate(zip(ps,levels,modes,es)):
            code=code*(labels(p,E,e) if mode=='tail' else p**e)+proj[i,e]
        maxima.append(int(np.bincount(code).max()))
    numerator=sum(maxima)
    result['uniform_capacity']=dict(numerator=numerator,denominator=len(points),
                                    required_new_inputs=count,excluded=numerator<count*len(points))
    print(json.dumps(dict(event='uniform_capacity',**result['uniform_capacity'])),flush=True)
    if result['uniform_capacity']['excluded']:
        result['status']='UNIFORM_CAPACITY_EXCLUDES_THIS_MENU'
        result['elapsed_seconds']=time.monotonic()-started
        return result
    model=cp_model.CpModel()
    colors=[model.new_int_var(0,count-1,f'color_{j}') for j in range(len(exponents))]
    vars=[]
    for j,es in enumerate(exponents):
        vars.append([model.new_int_var_from_domain(cp_model.Domain.from_values(
            list(map(int,np.unique(proj[i,e])))),f'r_{j}_{i}') for i,e in enumerate(es)])
    model.add(colors[0]==0)
    high=colors[0]
    for j in range(1,len(colors)):
        model.add(colors[j]<=high+1)
        nxt=model.new_int_var(0,count-1,f'high_{j}')
        model.add_max_equality(nxt,[high,colors[j]])
        high=nxt
    witnesses=set()
    def constrain(index):
        if index in witnesses:
            return
        witnesses.add(index)
        for d in range(count):
            chosen=model.new_int_var(0,len(exponents)-1,f'w_{index}_{d}')
            model.add_element(chosen,colors,d)
            for i in range(len(ps)):
                constants=[int(proj[i,es[i]][index]) for es in exponents]
                value=model.new_int_var(0,max(constants),f'v_{index}_{d}_{i}')
                model.add_element(chosen,[vs[i] for vs in vars],value)
                model.add_element(chosen,constants,value)
    for index in sorted({0,len(points)//2,len(points)-1}):
        constrain(index)
    trace=[]
    while time.monotonic()-started<seconds:
        solver=cp_model.CpSolver()
        solver.parameters.max_time_in_seconds=min(30,max(.1,seconds-time.monotonic()+started))
        solver.parameters.num_search_workers=4
        status=solver.solve(model)
        info=dict(iteration=len(trace)+1,witnesses=len(witnesses),status=solver.status_name(status))
        trace.append(info)
        if status not in (cp_model.OPTIMAL,cp_model.FEASIBLE):
            print(json.dumps(info),flush=True)
            break
        inputs=[[] for _ in range(count)]
        for j,es in enumerate(exponents):
            inputs[solver.value(colors[j])].append(dict(exponents=list(es),
                residues=[solver.value(v) for v in vars[j]]))
        new=set()
        gaps=[]
        for rows in inputs:
            covered=np.zeros(len(points),dtype=bool)
            for row in rows:
                covered|=row_hit(row,proj,len(points))
            missing=np.flatnonzero(~covered)
            gaps.append(len(missing))
            if len(missing):
                new.add(int(missing[0]))
        info['exhaustive_input_gaps']=gaps
        print(json.dumps(info),flush=True)
        if not new:
            result['inputs']=list(frozen_inputs)+inputs
            result['verification']=validate(result)
            result['status']=result['verification']['status']
            break
        assert not new&witnesses
        for index in sorted(new):
            constrain(index)
        model.clear_hints()
        for var in colors+[v for vs in vars for v in vs]:
            model.add_hint(var,solver.value(var))
    result['trace']=trace
    result['elapsed_seconds']=time.monotonic()-started
    result['scope']='Bounded search of a finite resource menu; no general impossibility follows from failure.'
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target',choices=['odd','even','all'],default='odd')
    parser.add_argument('--root',type=int,default=17)
    parser.add_argument('--seconds',type=float,default=150)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    region=dict(exponents=[0 if args.target=='all' else 1,0,0,0],
                residues=[1 if args.target=='odd' else 0,0,0,0])
    result=search(source_record(),region,args.root,args.seconds)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('prior','inputs','trace')},indent=2))
