"""Experimental fine moment-order scan; no values here certify a bound."""
import argparse,json,math
from pathlib import Path
import numpy as np
from minimum_modulus_experiment import sieve


def evaluate(w,numerator,denominator=32):
    n=w.get('terminal_index',w['start']+len(w['deltas']));s=numerator/denominator
    ps=np.asarray(sieve(max(2000000,20*n))[:n],dtype=float)
    ds=np.asarray([0]*w['start']+w['deltas'][:n-w['start']])/1e12
    nu=1/(1-ds);inc=np.zeros(n);qp=1/ps
    for a in range(1,120):
        inc+=((a+1)**s-a**s)*qp;qp/=ps
    U=math.exp(float(np.log1p(nu*inc).sum()))
    P=float(ps[-1]);ell=math.log(P)-1.1;r=P/(P-1);d=.5
    alpha=(2**s-1)*r**3/(1-d);m=math.ceil(alpha);theta=alpha-(m-1)
    z=alpha*(1/ell-n/P);beta=s-1
    def J(j):
        t=total=1.
        for h in range(1,j+1):t*=(j-h+1)/(beta*ell);total+=t
        return total/(beta*P**beta)
    def Q(j):return s*J(j-1)/ell-j*J(j-2)/ell**2-n/P**s
    kappa=beta**beta/(s**s*d**beta*(1-d))
    tail=kappa*U*r**s/(1-z)*((1-theta)*Q(m-1)+theta*Q(m))
    return dict(numerator=numerator,denominator=denominator,order=s,moment=U,tail=tail)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness');ap.add_argument('--output',required=True);a=ap.parse_args()
    w=json.loads(Path(a.witness).read_text());rows=[evaluate(w,j) for j in range(80,97)]
    result=dict(status='EXPERIMENTAL',witness=a.witness,rows=rows,best=min(rows,key=lambda x:x['tail']))
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['best']))
