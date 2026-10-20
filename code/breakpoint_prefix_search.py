"""Experimental early-coordinate search at exact loss breakpoints.

The late adjoint ranks trials, including conditioned-state weight changes.
The final full evaluation must improve before a candidate is retained.
This is a finite local candidate scan, not a global optimum certificate.
"""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np
from early_exact_states import states,specification
from joint_distortion_search import JointProblem
from global_distortion_search import GlobalProblem


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('witness')
    ap.add_argument('--output',required=True);ap.add_argument('--candidates',type=int,default=16)
    a=ap.parse_args();w=json.loads(Path(a.witness).read_text())
    p=JointProblem(w,K=w['K'],terminal=w['terminal_index'],freeze=150,knots=12)
    x=np.zeros(p.dimension);ne=p.early_count;p.set_prefix(x[:ne])
    base,_,_=GlobalProblem.evaluate(p,x[ne:],True)
    AP,AM,AR,AC,AU,potential=p.prefix_sensitivity
    P0,M0,R0=(v.copy() for v in p.initial);mu0,C0,U0,prob0=p.mu0,p.C0,p.U0,p.sm.prob.copy()
    def score(xx):
        p.set_prefix(xx[:ne]);P,M,R=p.initial
        return base-(p.mu0-mu0)+AP@(P-P0)+AM@(M-M0)+AR@(R-R0)+AC*(p.C0-C0)+AU*(p.U0-U0)+potential@(p.sm.prob-prob0)
    current=base;history=[]
    for j in range(p.small-p.original_start):
        if j in p.fixed_early_indices:continue
        idx=p.original_start+j;q=p.ps[idx];_,den,_=specification(p.K,q)
        critical={Fraction(0),Fraction(1,2)}
        critical.update(Fraction(v,den*(q-1)) for _,_,v in states(p.K,q,p.ps[:idx],q-1))
        critical={v for v in critical if 0<=v<=Fraction(1,2)}
        original=float(p.early_seed[j]+p.early_scales[j]*x[j])
        values=sorted(critical,key=lambda v:abs(math.log((float(v)+1e-12)/(original+1e-12))))[:a.candidates]
        values=list(map(float,values))+[original,0.,min(.5,2*original),original/2]
        best=(current,x[j]);trials=[]
        for d in sorted(set(values)):
            xx=x.copy();xx[j]=(d-p.early_seed[j])/p.early_scales[j]
            value=float(score(xx));trials.append(dict(delta=d,surrogate=value))
            if value<best[0]:best=(value,xx[j])
        current,x[j]=best
        row=dict(prime=q,critical_values=len(critical),selected_delta=float(p.early_seed[j]+p.early_scales[j]*x[j]),surrogate=float(current),trials=trials)
        history.append(row);print(q,row['selected_delta'],current,flush=True)
    actual,record=p.evaluate(x,save=True);accepted=actual<base
    if not accepted:actual,record=p.evaluate(np.zeros(p.dimension),save=True)
    record.update(breakpoint_history=history,breakpoint_source=a.witness,
                  initial_objective=float(base),surrogate_final=float(current),
                  candidate_accepted=bool(accepted),
                  search_note='Finite breakpoint scan with a full final evaluation; not global optimization or a proof.')
    Path(a.output).write_text(json.dumps(record)+'\n')
    print('Full final objective',actual,'accepted',accepted,flush=True)
