"""Experimental restricted stop-loss bound conditioned on small valuations."""
import argparse,json,math
from pathlib import Path
import numpy as np
from minimum_modulus_experiment import sieve,factors_table,initial,moments
from stoploss_experiment import update
from binned_distribution import Binned


class Small:
    def __init__(self,K,p,exponents=(32,22)):
        aa,bb=np.meshgrid(np.arange(exponents[0]+1),np.arange(exponents[1]+1),indexing='ij')
        self.a=aa.ravel();self.b=bb.ravel()
        self.tau=(self.a+1)*(self.b+1)
        self.prob=2.**(-self.a-1)*2*3.**(-self.b-1)
        self.mean=3.
        self.tail=self.mean-float(np.dot(self.prob,self.tau))
        self.count=np.zeros((len(self.a),(K-1)//p+2),dtype=np.int16)
        nmax=self.count.shape[1]-1
        a=0;x=1
        while x<nmax:
            b=0;d=x
            while d<nmax:
                self.count[:,d+1]+=((self.a>=a)&(self.b>=b))
                b+=1;d*=3
            a+=1;x*=2
        self.count=np.cumsum(self.count,axis=1)

    def excluded(self,K,p):
        B=np.zeros(len(self.a));power=p
        while power<K:
            s=(K-1)//power+1
            B+=(p-1)/power*self.count[:,s]
            power*=p
        return B

    def loss(self,d,p,B,Crest,cdf,wcdf,values=None):
        t=(B+d*(p-1))/self.tau
        t=np.minimum(t,len(cdf)-1 if values is None else values[-1])
        j=t.astype(np.int64) if values is None else np.searchsorted(values,t,side='right')-1
        e=Crest-t+t*cdf[j]-wcdf[j]
        return (float(np.dot(self.prob*self.tau,e))+Crest*self.tail)/((p-1)*(1-d))


def run(K=40000,start=20,cap=20000,limit=1000000,record=False,scale=1.,resolution=0,min_stop=0):
    ps=sieve(limit);fs=factors_table((K-1)//ps[start]+1)
    sm=Small(K,ps[start]);mu,C,T=initial(K,ps,start);nu={}
    binned=Binned(cap,resolution) if resolution else None
    n=np.arange(cap+1,dtype=float) if binned is None else binned.values
    upd=update if binned is None else binned.update
    dist=np.zeros(len(n));dist[1]=1.
    for q in ps[2:start]:dist=upd(dist,q,1.)
    rows=[];ds=[];wins={'first':0,'second':0,'conditional':0}
    for idx,p in enumerate(ps[start:],start):
        if mu<=0:return dict(K=K,start=start,success=False,p=p,mu=mu,wins=wins)
        a=(3*p-1)/(p-1)**2
        m1,m2=moments(K,p,C,T,nu,fs) if p<K else (C/(p-1),T/(p-1)**2)
        cdf=np.cumsum(dist);wcdf=np.cumsum(dist*n);B=sm.excluded(K,p)
        def score(d):
            loss=sm.loss(d,p,B,C/3,cdf,wcdf,n if binned else None)
            return T*(1+a/(1-d))/(mu-loss) if 0<=loss<mu else math.inf
        lo=0.;hi=.5;g=(5**.5-1)/2;x=hi-g*(hi-lo);y=lo+g*(hi-lo)
        fx=score(x);fy=score(y)
        for _ in range(22):
            if fx<fy:hi,y,fy=y,x,fx;x=hi-g*(hi-lo);fx=score(x)
            else:lo,x,fx=x,y,fy;y=lo+g*(hi-lo);fy=score(y)
        d=min(.5,(x if fx<fy else y)*scale)
        choices=[]
        for dd in (d,.5):
            loss=sm.loss(dd,p,B,C/3,cdf,wcdf,n if binned else None)
            if 0<=loss<mu:choices.append((T*(1+a/(1-dd))/(mu-loss),dd,loss,'conditional'))
        if m1<mu:choices.append((T*(1+a)/(mu-m1),0.,m1,'first'))
        if m2<mu:
            dd=(1+a)/(1+math.sqrt(1+4*mu*a*(1+a)/m2))
            loss=m2/(4*dd*(1-dd))
            if loss<mu:choices.append((T*(1+a/(1-dd))/(mu-loss),dd,loss,'second'))
        if not choices:return dict(K=K,start=start,success=False,p=p,mu=mu,wins=wins)
        _,d,loss,kind=min(choices);mu-=loss;wins[kind]+=1
        if record:ds.append(round(d*10**12))
        v=1/(1-d);C*=1+v/(p-1);T*=1+a*v;nu[p]=v
        dist=upd(dist,p,v)
        i=idx+1;f=T/mu;stop=i*(math.log(i)+math.log(math.log(i))-3)**2
        if i in (100,1000,10000,51000) or (f<stop and i>=min_stop):rows.append(dict(i=i,p=p,mu=mu,f=f,stop=stop,wins=wins.copy()))
        if f<stop and i>=min_stop:
            r=dict(K=K,start=start,cap=cap,resolution=resolution,success=True,rows=rows,wins=wins,small_exponents=[32,22])
            if record:r['deltas']=ds
            return r
    return dict(K=K,start=start,cap=cap,success=False,rows=rows,wins=wins)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--K',type=int,default=40000)
    ap.add_argument('--start',type=int,default=20);ap.add_argument('--cap',type=int,default=20000)
    ap.add_argument('--limit',type=int,default=1000000);ap.add_argument('--scale',type=float,default=1.)
    ap.add_argument('--resolution',type=int,default=0);ap.add_argument('--min-stop',type=int,default=0)
    ap.add_argument('--output');args=vars(ap.parse_args());out=args.pop('output')
    r=run(**args,record=bool(out))
    if out:Path(out).write_text(json.dumps(r))
    print(json.dumps({k:v for k,v in r.items() if k!='deltas'},indent=2))
