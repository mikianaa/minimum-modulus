"""Experimental convex-comparison loss bound. Requires NumPy; NOT a proof."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from minimum_modulus_experiment import sieve,factors_table,initial,moments


def update(dist,q,nu):
    cap=len(dist)-1
    new=(1-nu/q)*dist
    weight=nu*(q-1)/q**2
    h=2
    while h<=cap and weight>1e-30:
        new[h::h]+=weight*dist[1:cap//h+1]
        weight/=q;h+=1
    return new


def run(K,start,cap=10000,limit=2000000,correction=0.,record=False):
    ps=sieve(limit);fs=factors_table((K-1)//ps[start]+1)
    mu,C,T=initial(K,ps,start);nu={}
    dist=np.zeros(cap+1);dist[1]=1.
    for q in ps[:start]:dist=update(dist,q,1.)
    n=np.arange(cap+1,dtype=float)
    rows=[];ds=[];wins={'first':0,'second':0,'stoploss':0}
    for idx,p in enumerate(ps[start:],start):
        if mu<=0:return dict(K=K,start=start,success=False,p=p,mu=mu,wins=wins)
        a=(3*p-1)/(p-1)**2
        m1,m2=moments(K,p,C,T,nu,fs) if p<K else (C/(p-1),T/(p-1)**2)
        power=p
        while power<K:power*=p
        shift=1-p/power
        # Minimize the next T/mu for the first and second moment bounds.
        choices=[]
        if m1<mu:choices.append((T*(1+a)/(mu-m1),0.,m1,'first'))
        if m2<mu:
            d=(1+a)/(1+math.sqrt(1+4*mu*a*(1+a)/m2))
            d*=1-correction/math.sqrt(p)
            d=min(.5,max(1e-10,d));loss=m2/(4*d*(1-d))
            if loss<mu:choices.append((T*(1+a/(1-d))/(mu-loss),d,loss,'second'))
        if p<=2*cap:
            cdf=np.cumsum(dist);weighted=np.cumsum(dist*n)
            # For this bound, T_new/mu_new is linear-fractional between
            # consecutive integer breakpoints, so endpoints suffice.
            d=(n-shift)/(p-1)
            valid=(d>0)&(d<=.5)
            loss=(C-n+n*cdf-weighted)/((p-1)*(1-d))
            valid&=(loss>=0)&(loss<mu)
            values=np.full(cap+1,np.inf)
            values[valid]=T*(1+a/(1-d[valid]))/(mu-loss[valid])
            j=int(np.argmin(values))
            if math.isfinite(values[j]):choices.append((float(values[j]),float(d[j]),float(loss[j]),'stoploss'))
            # Include delta=1/2 even when its breakpoint is nonintegral.
            t=(p-1)*.5+shift;j=int(t)
            if j<=cap:
                loss=(C-t+t*cdf[j]-weighted[j])/((p-1)*.5)
                if 0<=loss<mu:choices.append((T*(1+2*a)/(mu-loss),.5,loss,'stoploss'))
        if not choices:return dict(K=K,start=start,success=False,p=p,mu=mu,wins=wins)
        _,d,loss,kind=min(choices)
        wins[kind]+=1;mu-=loss
        if record:ds.append(d)
        v=1/(1-d);C*=1+v/(p-1);T*=1+a*v;nu[p]=v
        if p<=2*cap:dist=update(dist,p,v)
        i=idx+1;f=T/mu;stop=i*(math.log(i)+math.log(math.log(i))-3)**2
        if i in (100,1000,10000,51000) or f<stop:rows.append(dict(i=i,p=p,mu=mu,f=f,stop=stop,wins=wins.copy()))
        if f<stop:
            r=dict(K=K,start=start,cap=cap,success=True,rows=rows,wins=wins)
            if record:r['deltas']=ds
            return r
    return dict(K=K,start=start,cap=cap,success=False,rows=rows,wins=wins)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--K',type=int,default=40000)
    ap.add_argument('--start',type=int,default=25);ap.add_argument('--cap',type=int,default=10000)
    ap.add_argument('--limit',type=int,default=2000000);ap.add_argument('--correction',type=float,default=0.)
    ap.add_argument('--output');args=vars(ap.parse_args());out=args.pop('output')
    r=run(**args,record=bool(out))
    if out:Path(out).write_text(json.dumps(r))
    print(json.dumps({k:v for k,v in r.items() if k!='deltas'},indent=2))
