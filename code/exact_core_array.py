"""Exhaustive finite-kernel core construction using exact integer marginals."""
import argparse
from fractions import Fraction
import json
from math import prod
from pathlib import Path
import time

from graded_arrow_blocks import np
from construction_kernel import resources,labels,capacity,verify
from structural_cover_search import coordinates,row_from_condition,condition_from_row
from verify_construction_searches import verify_uncovered_progression


def build(minimum,ps,ds,cap=50000,order='mass',focus_prime=None,focus_residue=0,focus_weight=16,seed_rows=None):
    start=time.monotonic()
    domains,blocks=coordinates(ps,ds)
    total=prod(domains)
    residual=np.ones(domains,dtype=bool)
    weights=None
    focus_slice=None
    if focus_prime is not None:
        i=ps.index(focus_prime)
        assert ds[i]>=2 and 0<=focus_residue<focus_prime and focus_weight>=1
        weights=np.ones(domains,dtype=np.int16)
        focus_slice=tuple(focus_residue if v==blocks[i][0] else slice(None) for v in range(len(domains)))
        weights[focus_slice]=focus_weight
    exps=list(resources(ps,ds,minimum,cap))
    all_exps=exps[:]
    seed_rows=seed_rows or []
    if seed_rows:
        verify(dict(minimum=minimum,primes=ps,depths=ds,rows=seed_rows))
        seen={tuple(r['exponents']) for r in seed_rows}
        assert seen<=set(exps)
        exps=[e for e in exps if e not in seen]
    base=lambda es:prod(p**e for p,e in zip(ps,es))
    volume=lambda es:prod(labels(p,E,e) for p,E,e in zip(ps,ds,es))
    if order=='mass':
        exps.sort(key=lambda es:(volume(es),base(es)))
    elif order=='base':
        exps.sort(key=base)
    else:
        exps.sort(key=lambda es:(max(p for p,e in zip(ps,es) if e),base(es)))
    rows=list(seed_rows)
    for row in rows:
        assignment=condition_from_row(row,ps,ds,blocks)
        residual[tuple(assignment.get(v,slice(None)) for v in range(len(domains)))]=False
    before=int(residual.sum())
    trace=[]
    for j,es in enumerate(exps):
        selected=sorted(v for block,e in zip(blocks,es) for v in block[:e])
        summed=tuple(v for v in range(len(domains)) if v not in selected)
        objective=residual if weights is None else residual*weights
        marginal=objective.sum(axis=summed,dtype=np.int64) if summed else objective
        score=int(marginal.max())
        if score:
            values=np.unravel_index(int(marginal.argmax()),marginal.shape)
            assignment=dict(zip(selected,map(int,values)))
            sl=tuple(assignment.get(v,slice(None)) for v in range(len(domains)))
            assert int(objective[sl].sum())==score
            gain=int(residual[sl].sum())
            residual[sl]=False
            rows.append(row_from_condition(es,assignment,ps,ds,blocks))
            before-=gain
        if j%25==0 or before==0:
            info=dict(processed=j+1,rows=len(rows),uncovered=before,
                      fraction=str(Fraction(before,total)),seconds=round(time.monotonic()-start,3))
            trace.append(info)
            print(json.dumps(info),flush=True)
        if before==0:
            break
    assert int(residual.sum())==before
    # Independent reconstruction, separate from placement marginals.
    rebuilt=np.ones(domains,dtype=bool)
    for row in rows:
        assignment=condition_from_row(row,ps,ds,blocks)
        rebuilt[tuple(assignment.get(v,slice(None)) for v in range(len(domains)))]=False
    assert np.array_equal(rebuilt,residual)
    result=dict(minimum=minimum,primes=ps,depths=ds,rows=rows,trace=trace,
                method='exhaustive_core_integer_marginals',order=order,
                maximum_resource_base=cap,resource_count=len(all_exps),capacity=str(capacity(ps,ds,all_exps)),
                seeded_resources=len(seed_rows),
                kernel_states=total,exact_uncovered_states=before,
                exact_uncovered_fraction=str(Fraction(before,total)),
                elapsed_seconds=time.monotonic()-start)
    if focus_prime is not None:
        result['structural_focus']=dict(prime=focus_prime,residue=focus_residue,weight=focus_weight,
                                        remaining_in_progression=int(residual[focus_slice].sum()),
                                        progression_kernel_states=total//focus_prime)
    result['verification']=verify(result)
    result['status']=result['verification']['status']
    if result['status']=='UNCOVERED':
        result['independent_witness_check']=verify_uncovered_progression(result)
    if before==0:
        assert result['status']=='VERIFIED'
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum',type=int,default=50)
    parser.add_argument('--primes',default='2,3,5,7')
    parser.add_argument('--depths',default='6,4,3,2')
    parser.add_argument('--cap',type=int,default=50000)
    parser.add_argument('--order',choices=['mass','base','prime-stage'],default='mass')
    parser.add_argument('--output',required=True)
    parser.add_argument('--focus-prime',type=int)
    parser.add_argument('--focus-residue',type=int,default=0)
    parser.add_argument('--focus-weight',type=int,default=16)
    args=parser.parse_args()
    r=build(args.minimum,list(map(int,args.primes.split(','))),list(map(int,args.depths.split(','))),args.cap,args.order,
            args.focus_prime,args.focus_residue,args.focus_weight)
    Path(args.output).write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['status','kernel_states','exact_uncovered_states','exact_uncovered_fraction','elapsed_seconds']},indent=2))
