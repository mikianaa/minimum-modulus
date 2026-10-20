"""Experimental coordinate searches using an exact fixed-suffix objective.

Only prefix primes outside the conditioned set vary. Thus the conditioned
weights and all suffix distortion parameters are fixed, and the suffix
objective is linear in the prefix distributions and moments. One adjoint
supplies the objective for every subsequent coordinate trial. The final
complete evaluation checks this identity numerically; proof is independent.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from joint_distortion_search import JointProblem
from global_distortion_search import GlobalProblem


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('witness');ap.add_argument('--output',required=True)
    ap.add_argument('--coordinates',type=int,default=12);ap.add_argument('--iterations',type=int,default=12)
    ap.add_argument('--window',type=float,default=.15);a=ap.parse_args()
    w=json.loads(Path(a.witness).read_text());p=JointProblem(w,K=w['K'],terminal=w['terminal_index'],freeze=150)
    assert not p.search_stratified_tail,'The fixed-suffix linearity used here requires a single-moment objective.'
    x=np.zeros(p.dimension);ne=p.early_count;p.set_prefix(x[:ne]);base,_,_=GlobalProblem.evaluate(p,x[ne:],True)
    AP,AM,AR,AC,AU,_=p.prefix_sensitivity
    P0,M0,R0=(v.copy() for v in p.initial);mu0,C0,U0=p.mu0,p.C0,p.U0
    def objective(xx):
        p.set_prefix(xx[:ne]);P,M,R=p.initial
        return base-(p.mu0-mu0)+AP@(P-P0)+AM@(M-M0)+AR@(R-R0)+AC*(p.C0-C0)+AU*(p.U0-U0)
    g=w['history'][-1]['gradient']
    indices=sorted(range(p.small-p.original_start,p.individual),key=lambda j:-abs(g[j]))[:a.coordinates]
    current=base;history=[];ratio=(math.sqrt(5)-1)/2
    for j in indices:
        origin=x[j];lo=max(origin-a.window,-p.early_seed[j]/p.early_scales[j]);hi=min(origin+a.window,(.5-p.early_seed[j])/p.early_scales[j])
        best=(current,origin)
        def trial(t):
            global best
            xx=x.copy();xx[j]=t;f=objective(xx)
            if f<best[0]:best=(f,t)
            return f
        trial(lo);trial(hi)
        u=hi-ratio*(hi-lo);v=lo+ratio*(hi-lo);fu=trial(u);fv=trial(v)
        for _ in range(a.iterations):
            if fu<fv:hi=v;v=u;fv=fu;u=hi-ratio*(hi-lo);fu=trial(u)
            else:lo=u;u=v;fu=fv;v=lo+ratio*(hi-lo);fv=trial(v)
        current,x[j]=best
        history.append(dict(prime=p.ps[p.original_start+j],objective=current,control=float(x[j])))
        print(history[-1],flush=True)
    actual,record=p.evaluate(x,save=True)
    assert abs(actual-current)<2e-6,(actual,current)
    record.update(coordinate_history=history,coordinate_source=a.witness,coordinate_surrogate_discrepancy=actual-current)
    Path(a.output).write_text(json.dumps(record)+'\n');print('Final objective',actual,'mu',record['mu'],flush=True)
