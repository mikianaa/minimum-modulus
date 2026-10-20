"""Joint exact packing of all inputs of one new largest-prime arrow.

Unlike sequential supplies, every input may be a union of cylinders in
several primes. All candidate fixed exponent resources are shared across
the inputs and may be used at most once. The universe is an exhaustive
projection of the core residual, not a sample. A partial branch success
does not imply that the complete system covers all integers.
"""
import argparse
from itertools import product
import json
from math import prod
from pathlib import Path
import time
from collections import defaultdict

from core_region_arrows import Core
from construction_kernel import prime,verify,labels
from graded_arrow_blocks import np
from ortools.sat.python import cp_model
from ortools.linear_solver import pywraplp


def target_points(certificate,region,levels,tail_mode=False):
    core=Core(certificate)
    assert len(levels)==len(core.ps)
    assert all(1<=e<=E if tail_mode else 0<=e<E for e,E in zip(levels,core.ds))
    residual=np.zeros(core.dims,dtype=bool)
    index=core.region_index(region)
    residual[index]=core.residual[index]
    if tail_mode:
        marginal=residual
        for i,(p,E) in enumerate(zip(core.ps,levels)):
            mapping=[]
            for x in core.values[i]:
                a=int(x)%p**(E-1)
                y=(int(x)-a)//p**(E-1)
                assert y>0
                while y%p==0:
                    y//=p
                mapping.append(a*(p-1)+(y%p)-1)
            mapping=np.array(mapping)
            marginal=np.stack([np.take(marginal,np.flatnonzero(mapping==s),axis=i).any(axis=i)
                               for s in range(p**(E-1)*(p-1))],axis=i)
    else:
        marginal=core.marginal(residual,levels)
    points=np.argwhere(marginal>0)
    return core,points,int(residual.sum())


def validate(result,input_count=None):
    certificate=result['core']
    verify(certificate)
    tail_mode=result.get('tail_mode',False)
    core,points,needed=target_points(certificate,result['region'],result['levels'],tail_mode)
    q=result['root_prime']
    assert prime(q) and q>max(core.ps)
    count=q-1 if input_count is None else input_count
    assert type(count) is int and 1<=count<=q-1
    assert len(result['inputs'])==count
    seen=set()
    checks=0
    for rows in result['inputs']:
        covered=np.zeros(len(points),dtype=bool)
        for row in rows:
            es,rs=row['exponents'],row['residues']
            assert len(es)==len(rs)==len(core.ps)
            assert all(type(e) is int and 0<=e<=L for e,L in zip(es,result['levels']))
            assert all(type(r) is int and 0<=r<(labels(p,E,e) if tail_mode else p**e)
                       for r,p,E,e in zip(rs,core.ps,result['levels'],es))
            assert tuple(es) not in seen
            seen.add(tuple(es))
            modulus=prod(p**e for p,e in zip(core.ps,es))
            assert q*modulus>=certificate['minimum']
            hit=np.ones(len(points),dtype=bool)
            for i,(p,E,e,r) in enumerate(zip(core.ps,result['levels'],es,rs)):
                if tail_mode:
                    hit&=(points[:,i]==r if e==E else (points[:,i]//(p-1)-r)%p**e==0)
                else:
                    hit&=(points[:,i]-r)%p**e==0
            covered|=hit
        assert covered.all()
        checks+=len(points)
    return dict(status='VERIFIED_RELATIVE_ROOT' if count==q-1 else 'VERIFIED_INPUT_PACKING',
                root_prime=q,input_count=count,projected_points=len(points),
                checked_input_points=checks,resource_count=len(seen),core_states_repaired=needed,
                scope=('Only the specified target region is repaired relative to the saved core.' if count==q-1 else
                       'These resource-disjoint inputs do not fill all root digits and do not certify a complete root.'))


def validate_packing(result):
    return validate(result,result['required_inputs'])


def validate_capacity(result):
    """Rebuild maximal residue weights using Python integer histograms."""
    verify(result['core'])
    tail_mode=result.get('tail_mode',False)
    core,points,needed=target_points(result['core'],result['region'],result['levels'],tail_mode)
    q=result['root_prime']
    assert prime(q) and q>max(core.ps)
    record=result['exact_weighted_capacity']
    weights=record['integer_point_weights']
    assert len(weights)==len(points) and all(type(w) is int and w>=0 for w in weights)
    total=sum(weights)
    assert total>0 and record['denominator']==total
    capacity=0
    resources=0
    allowed=None
    if 'allowed_resources' in result:
        saved=result['allowed_resources']
        allowed={tuple(es) for es in saved}
        assert len(saved)==len(allowed)
        assert all(len(es)==len(core.ps) and all(type(e) is int and 0<=e<=E for e,E in zip(es,result['levels']))
                   and q*prod(p**e for p,e in zip(core.ps,es))>=result['core']['minimum'] for es in allowed)
    for es in product(*(range(e+1) for e in result['levels'])):
        if allowed is not None and es not in allowed:
            continue
        moduli=[p**e for p,e in zip(core.ps,es)]
        if q*prod(moduli)<result['core']['minimum']:
            continue
        histogram=defaultdict(int)
        for point,w in zip(points,weights):
            key=tuple((int(x) if e==E else (int(x)//(p-1))%m) if tail_mode else int(x)%m
                      for x,p,E,e,m in zip(point,core.ps,result['levels'],es,moduli))
            histogram[key]+=w
        capacity+=max(histogram.values())
        resources+=1
    count=result.get('required_inputs',q-1)
    assert type(count) is int and 1<=count<=q-1
    assert capacity==record['numerator'] and record['required_inputs']==count
    assert record['excluded']==(capacity<count*total)
    return dict(status='VERIFIED_RESTRICTED_CAPACITY',numerator=capacity,denominator=total,
                required_inputs=count,excluded=record['excluded'],resources=resources,
                scope='Only the saved exponent/tail resource menu and target are bounded.')


def search(certificate,region,levels,q,seconds=60,tail_mode=False,input_count=None,allowed_resources=None,capacity_only=False,seed_inputs=None):
    started=time.monotonic()
    verify(certificate)
    core,points,needed=target_points(certificate,region,levels,tail_mode)
    assert len(points) and prime(q) and q>max(core.ps)
    count=q-1 if input_count is None else input_count
    assert type(count) is int and 1<=count<=q-1
    assert seed_inputs is None or len(seed_inputs)==count
    allowed=None if allowed_resources is None else {tuple(es) for es in allowed_resources}
    metadata=dict(required_inputs=count)
    if allowed is not None:
        metadata['allowed_resources']=[list(es) for es in sorted(allowed)]
    rows=[]
    groups=[]
    hits=[]
    for es in product(*(range(e+1) for e in levels)):
        if allowed is not None and es not in allowed:
            continue
        moduli=np.array([p**e for p,e in zip(core.ps,es)],dtype=np.int64)
        if q*int(prod(map(int,moduli)))<certificate['minimum']:
            continue
        projected=(np.column_stack([points[:,i] if e==E else points[:,i]//(p-1)%p**e
                                    for i,(p,E,e) in enumerate(zip(core.ps,levels,es))])
                   if tail_mode else points%moduli)
        residues,inverse=np.unique(projected,axis=0,return_inverse=True)
        group=[]
        for j,rs in enumerate(residues):
            group.append(len(rows))
            rows.append(dict(exponents=list(es),residues=list(map(int,rs))))
            hits.append(np.flatnonzero(inverse==j))
        groups.append(group)
    uniform_capacity=sum(max(len(hits[r]) for r in group) for group in groups)
    if uniform_capacity<count*len(points):
        result=dict(core=certificate,region=region,levels=levels,root_prime=q,
                    tail_mode=tail_mode,method='joint_resource_root_packing',status='RESTRICTED_MENU_EXCLUDED',
                    candidate_rows=len(rows),resources=len(groups),projected_points=len(points),
                    core_states_targeted=needed,elapsed_seconds=time.monotonic()-started,
                    exact_weighted_capacity=dict(numerator=uniform_capacity,denominator=len(points),
                        integer_point_weights=[1]*len(points),required_inputs=count,excluded=True),
                    scope='Uniform exact capacity excludes only this exponent/tail menu on the specified target.')
        result.update(metadata)
        result['verification']=validate_capacity(result)
        return result
    # A dual weighting bounds the number of disjoint complete inputs. The
    # numerical optimizer only proposes weights; the saved inequality is
    # recomputed exactly after rounding to nonnegative integers.
    dual=pywraplp.Solver.CreateSolver('GLOP')
    dual.SetTimeLimit(int(seconds*1000))
    weights=[dual.NumVar(0,1,f'w_{x}') for x in range(len(points))]
    bounds=[dual.NumVar(0,1,f'b_{g}') for g in range(len(groups))]
    dual.Add(sum(weights)==1)
    for g,group in enumerate(groups):
        for r in group:
            dual.Add(bounds[g]>=sum(weights[int(x)] for x in hits[r]))
    dual.Minimize(sum(bounds))
    dual_status=dual.Solve()
    integer_weights=[1]*len(points)
    if dual_status==pywraplp.Solver.OPTIMAL:
        proposed=[max(0,round(w.solution_value()*10**9)) for w in weights]
        if sum(proposed)>0:
            integer_weights=proposed
    denominator=sum(integer_weights)
    numerator=sum(max(sum(integer_weights[int(x)] for x in hits[r]) for r in group) for group in groups)
    capacity=dict(numerator=numerator,denominator=denominator,integer_point_weights=integer_weights,
                  required_inputs=count,excluded=numerator<count*denominator)
    if capacity['excluded'] or capacity_only:
        result=dict(core=certificate,region=region,levels=levels,root_prime=q,
                    tail_mode=tail_mode,method='joint_resource_root_packing',
                    status='RESTRICTED_MENU_EXCLUDED' if capacity['excluded'] else 'NOT_EXCLUDED_BY_SAVED_WEIGHTS',
                    candidate_rows=len(rows),resources=len(groups),projected_points=len(points),
                    core_states_targeted=needed,exact_weighted_capacity=capacity,
                    elapsed_seconds=time.monotonic()-started,
                    scope='Exact weighted capacity excludes only this exponent/tail menu on the specified target; richer resource menus and other targets are not excluded.')
        result.update(metadata)
        result['verification']=validate_capacity(result)
        return result
    model=cp_model.CpModel()
    variables=[[model.new_bool_var(f'x_{r}_{d}') for d in range(count)] for r in range(len(rows))]
    for group in groups:
        model.add_at_most_one(variables[r][d] for r in group for d in range(count))
    covering=[[] for _ in points]
    for r,indices in enumerate(hits):
        for x in indices:
            covering[int(x)].append(r)
    for choices in covering:
        for d in range(count):
            model.add_bool_or(variables[r][d] for r in choices)
    # Inputs are interchangeable before assigning root digits. This ordering
    # removes only permutation symmetry, even if two scores happen to tie.
    if seed_inputs is None:
        scores=[sum((r+1)*variables[r][d] for r in range(len(rows))) for d in range(count)]
        for d in range(count-1):
            model.add(scores[d]<=scores[d+1])
    else:
        unchanged={(tuple(row['exponents']),tuple(row['residues']),d)
                   for d,classes in enumerate(seed_inputs) for row in classes}
        changed=[]
        for r,row in enumerate(rows):
            for d in range(count):
                key=tuple(row['exponents']),tuple(row['residues']),d
                if key in unchanged:
                    model.add_hint(variables[r][d],1)
                else:
                    changed.append(variables[r][d])
        model.minimize(sum(changed))
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=seconds
    solver.parameters.num_search_workers=4
    status=solver.solve(model)
    result=dict(core=certificate,region=region,levels=levels,root_prime=q,
                tail_mode=tail_mode,method='joint_resource_root_packing',solver_status=solver.status_name(status),
                candidate_rows=len(rows),resources=len(groups),projected_points=len(points),
                exact_weighted_capacity=capacity,
                core_states_targeted=needed,elapsed_seconds=time.monotonic()-started)
    result.update(metadata)
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        if seed_inputs is not None:
            result['changed_assignments']=round(solver.objective_value)
        result['inputs']=[[row for r,row in enumerate(rows) if solver.boolean_value(variables[r][d])]
                          for d in range(count)]
        result['verification']=validate_packing(result)
        result['status']=result['verification']['status']
    else:
        result['status']='NO_VERIFIED_ROOT'
        result['scope']='Solver result applies only to this finite exponent/tail resource menu; no general impossibility is asserted.'
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('core')
    parser.add_argument('--root',type=int,default=11)
    parser.add_argument('--seconds',type=float,default=60)
    parser.add_argument('--tail-mode',action='store_true')
    parser.add_argument('--levels',default='3,3,1,1')
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    # Nielsen's structural target: 1 mod 4 intersect (1 mod 3 union 3 mod 9).
    region=dict(exponents=[2,2,0,0],residues=[1,0,0,0],unions={'1':[1,3,4,7]})
    result=search(json.loads(Path(args.core).read_text()),region,list(map(int,args.levels.split(','))),args.root,args.seconds,args.tail_mode)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    summary={k:v for k,v in result.items() if k not in ('core','inputs','exact_weighted_capacity')}
    summary['capacity']={k:v for k,v in result['exact_weighted_capacity'].items() if k!='integer_point_weights'}
    print(json.dumps(summary,indent=2))
