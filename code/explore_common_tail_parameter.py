"""Floating discovery scan of one COMMON future distortion; not a proof.

Rebuild only the product distribution, since finite losses do not depend on
the future distortion. The completed integer report must certify all claims.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from binned_distribution import Binned
from minimum_modulus_experiment import sieve
from reachable_support import reachable_upper_support
from stratified_tail_search import StratifiedObjective


if __name__ == '__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('witness');ap.add_argument('--output',required=True)
    ap.add_argument('--lower',type=int,default=50)
    ap.add_argument('--upper',type=int,default=500)
    ap.add_argument('--step',type=int,default=10)
    a=ap.parse_args();raw=Path(a.witness).read_bytes();w=json.loads(raw)
    n=w['start']+len(w['deltas']);ps=sieve(max(2000000,20*n))[:n]
    assert len(ps)==n and 0<a.lower<=a.upper<=500 and a.step>0
    b=Binned(w['cap'],w['resolution'])
    P=np.zeros(len(b.values));P[1]=1.;M=P.copy();U=1.
    num,den=w['stratified_overflow_orders'][0];s=num/den
    for j,p in enumerate(ps):
        d=0. if j<w['start'] else w['deltas'][j-w['start']]/1e12
        v=1/(1-d);P,M=b.update_pair(P,M,p,v)
        term=1./p;h=1;increment=0.
        while term>1e-30:
            increment+=((h+1)**s-h**s)*term;term/=p;h+=1
        U*=1+v*increment
        if (j+1)%100000==0:print('Experimental product stage',j+1,flush=True)
    upper=reachable_upper_support(tuple(map(int,b.values)));rows=[]
    for d in range(a.lower,a.upper+1,a.step):
        objective=StratifiedObjective(b.values,upper,ps[-1],n,d*10**9,
                                      w['stratified_tail_orders'],(num,den))
        assert not objective.excluded_orders
        tail=objective.evaluate(P,M,U)[0]
        rows.append(dict(tail_delta_numerator=d*10**9,tail=tail))
    result=dict(status='EXPERIMENTAL',source=a.witness,
                source_sha256=hashlib.sha256(raw).hexdigest(),terminal_index=n,
                terminal_prime=ps[-1],best_tested=min(rows,key=lambda row:row['tail']),
                rows=rows,note='No finite losses replayed; no modulus conclusion or global optimality claim.')
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
