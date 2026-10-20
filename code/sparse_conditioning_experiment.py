"""Experimental sparse conditioning with an unrestricted fallback on the tail.

The finite conditioning set is D_small < state_cap, rather than an exponent
box. Omitted states keep the unrestricted convex bound, not a first-moment
tail. Only a nonnegative, bounded correction is subtracted.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from minimum_modulus_experiment import sieve,initial,moments,factors_table
from binned_distribution import Binned,MomentPrefixes
_cached_states={}


class SparseStates:
    def __init__(self,K,p,small_primes,state_cap,normalize=False):
        # Enumerate smooth integers directly: memory depends on the number of
        # retained states, rather than on the (much larger) numerical cutoff.
        states=[]
        def visit(j,n,tau):
            if j==len(small_primes):
                states.append((n,tau));return
            q=small_primes[j];power=1;a=0
            while n*power<state_cap:
                visit(j+1,n*power,tau*(a+1));power*=q;a+=1
        visit(0,1,1);states.sort()
        self.ns=np.asarray([n for n,t in states],dtype=np.int64)
        self.tau=np.asarray([t for n,t in states],dtype=np.int32)
        self.prob=math.prod(1-1/q for q in small_primes)/self.ns
        if normalize:self.prob/=np.sum(self.prob)
        maximum=(K-1)//p+1
        thresholds=set()
        for q in sieve(K):
            if q<p:continue
            power=q
            while power<K:thresholds.add((K-1)//power+1);power*=q
        self.count={};counts=np.zeros(len(self.ns),dtype=np.int16);d=1
        for s in sorted(thresholds):
            while d<s:counts+=(self.ns%d==0);d+=1
            self.count[s]=counts.copy()
        self.cache={}

    def groups(self,K,p):
        powers=[];power=p
        while power<K:powers.append(power);power*=p
        ss=tuple((K-1)//pw+1 for pw in powers)
        if ss not in self.cache:
            keys,inv=np.unique(np.column_stack([self.tau]+[self.count[s] for s in ss]),axis=0,return_inverse=True)
            self.cache[ss]=(keys,np.bincount(inv,weights=self.prob))
        keys,weights=self.cache[ss];B=np.zeros(len(keys))
        for j,pw in enumerate(powers):B+=(p-1)/pw*keys[:,j+1]
        return keys[:,0],B,weights


class DynamicSparseStates:
    def __init__(self,K,ps,start,small,state_cap):
        static=SparseStates(K,ps[start],ps[:small],state_cap)
        self.ns=static.ns;self.tau=static.tau;self.largest=np.ones(len(self.ns),dtype=np.int32)
        for q in ps[:small]:self.largest[self.ns%q==0]=q
        self.coeff=1./self.ns;self.normalizer=1.;self.last=ps[small-1]
        thresholds=set()
        for p in ps[start:]:
            if p>=K:break
            power=p
            while power<K:thresholds.add((K-1)//power+1);power*=p
        self.count={s:static.count[s] for s in thresholds}
        self.cache={};self.weight_cache={}
        for q in ps[:start]:self.update(q,1.)

    def update(self,q,nu):
        if q>self.last:return
        self.weight_cache.clear()
        self.coeff[self.ns%q==0]*=nu*(q-1)/(q-nu)
        self.normalizer*=1-nu/q

    def groups(self,K,p):
        powers=[];power=p
        while power<K:powers.append(power);power*=p
        ss=tuple((K-1)//pw+1 for pw in powers)
        cachekey=(min(p,self.last+1),ss)
        if cachekey in self.weight_cache:
            keys,weights=self.weight_cache[cachekey]
            B=np.zeros(len(keys))
            for j,pw in enumerate(powers):B+=(p-1)/pw*keys[:,j+1]
            return keys[:,0],B,weights
        active=self.largest<p
        if cachekey not in self.cache:
            keys,inv=np.unique(np.column_stack([self.tau[active]]+[self.count[s][active] for s in ss]),axis=0,return_inverse=True)
            self.cache[cachekey]=(keys,inv)
        keys,inv=self.cache[cachekey]
        weights=np.bincount(inv,weights=self.normalizer*self.coeff[active])
        if p>self.last:self.weight_cache[cachekey]=(keys,weights)
        B=np.zeros(len(keys))
        for j,pw in enumerate(powers):B+=(p-1)/pw*keys[:,j+1]
        return keys[:,0],B,weights


def run(K=10000,start=12,small=12,state_cap=1000000,cap=1000000,resolution=256,
        policy=1.,limit=2000000,max_steps=100000,min_stop=0,record=False,dynamic=False,
        support_mode='uniform',joint_presence=False):
    ps=sieve(limit);assert 2<=small and start>=2 and (dynamic or small<=start)
    key=(K,start,small,state_cap,dynamic)
    if not dynamic and key not in _cached_states:
        _cached_states.clear();_cached_states[key]=SparseStates(K,ps[start],ps[:small],state_cap)
    sm=DynamicSparseStates(K,ps,start,small,state_cap) if dynamic else _cached_states[key]
    mu,C,T=initial(K,ps,start);initial_mu=mu;nu={};fs=factors_table((K-1)//ps[start]+1)
    bins=Binned(cap,resolution);values=bins.values
    dist=np.zeros(len(values));dist[1]=1.;moment=dist.copy();rest=dist.copy()
    for j,q in enumerate(ps[:start]):
        dist,moment=bins.update_pair(dist,moment,q,1.)
        if j>=small:rest=bins.update(rest,q,1.)
    ds=[];rows=[];best_ratio=math.inf
    joint=None
    for idx,p in enumerate(ps[start:],start):
        if joint_presence and idx==small:
            from bernoulli_conditioning_experiment import JointStates
            if dynamic:sm.prob=sm.normalizer*sm.coeff
            joint=JointStates(sm,K,p,p)
        a=(3*p-1)/(p-1)**2;mp=MomentPrefixes(values,dist,moment,resolution,support_mode)
        P=np.cumsum(rest);W=np.cumsum(rest*values)
        tau,B,weights=sm.groups(K,p)
        def loss(d):
            A=(p-1)*d;t=min(A,cap)
            baseline=C-t+float(mp.complement(np.asarray([t]))[0])
            t0=A/tau;t1=(A+B)/tau
            j=np.minimum(np.searchsorted(values,t0,side='right')-1,len(values)-1)
            k=np.minimum(np.searchsorted(values,t1,side='right')-1,len(values)-1)
            gain=tau*(W[k]-W[j])-A*(P[k]-P[j])+B*(P[-1]-P[k])
            correction=float(np.sum(weights*np.maximum(0.,gain)))
            if joint is not None and p*ps[small]<K:
                correction=max(correction,joint.correction(p,d,ps[small:idx],[nu[q] for q in ps[small:idx]]))
            return (baseline-correction)/((p-1)*(1-d))
        def score(d):
            ll=loss(d)
            return math.log1p(a/(1-d))-policy*math.log1p(-ll/mu) if 0<=ll<mu else math.inf
        lo=0.;hi=.5;g=(5**.5-1)/2;x=hi-g*hi;y=g*hi;fx=score(x);fy=score(y)
        for _ in range(20):
            if fx<fy:hi,y,fy=y,x,fx;x=hi-g*(hi-lo);fx=score(x)
            else:lo,x,fx=x,y,fy;y=lo+g*(hi-lo);fy=score(y)
        candidates=[(score(d),d,loss(d)) for d in (0.,x if fx<fy else y,.5)]
        if p<K:
            m1,_=moments(K,p,C,T,nu,fs)
            if 0<=m1<mu:candidates.append((math.log1p(a)-policy*math.log1p(-m1/mu),0.,m1))
        val,d,ll=min(candidates)
        if not math.isfinite(val):
            return dict(K=K,start=start,small=small,state_cap=state_cap,cap=cap,
                resolution=resolution,support_mode=support_mode,dynamic=dynamic,
                joint_presence=joint_presence,policy=policy,success=False,p=p,mu=mu,
                rows=rows,best_ratio=best_ratio,deltas=ds)
        mu-=ll
        if record:ds.append(round(d*10**12))
        v=1/(1-d);C*=1+v/(p-1);T*=1+a*v;nu[p]=v
        dist,moment=bins.update_pair(dist,moment,p,v)
        if idx>=small:rest=bins.update(rest,p,v)
        if dynamic:sm.update(p,v)
        i=idx+1;f=T/mu;stop=i*(math.log(i)+math.log(math.log(i))-3)**2
        best_ratio=min(best_ratio,f/stop)
        if i in (100,1000,10000,50000,100000) or (f<stop and i>=min_stop):
            rows.append(dict(i=i,p=p,mu=mu,f=f,stop=stop,ratio=f/stop))
        if (f<stop and i>=min_stop) or i>=max_steps:break
    r=dict(K=K,start=start,small=small,state_cap=state_cap,cap=cap,resolution=resolution,
           policy=policy,dynamic=dynamic,support_mode=support_mode,joint_presence=joint_presence,success=f<stop,initial_mu=initial_mu,rows=rows,
           terminal=dict(i=i,p=p,mu=mu,f=f,stop=stop,ratio=f/stop),
           best_ratio=best_ratio,states=len(sm.ns))
    if record:r['deltas']=ds
    return r


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name,default in [('K',10000),('start',12),('small',12),('state_cap',1000000),
                         ('cap',1000000),('resolution',256),('limit',2000000),
                         ('max_steps',100000),('min_stop',0)]:
        ap.add_argument('--'+name.replace('_','-'),type=int,default=default)
    ap.add_argument('--policy',type=float,default=1.);ap.add_argument('--output')
    ap.add_argument('--dynamic',action='store_true')
    ap.add_argument('--joint-presence',action='store_true')
    ap.add_argument('--support-mode',choices=['uniform','reachable'],default='uniform')
    args=vars(ap.parse_args());out=args.pop('output');t=time.time();r=run(**args,record=bool(out))
    r['seconds']=time.time()-t
    if out:Path(out).write_text(json.dumps(r))
    print(json.dumps({k:v for k,v in r.items() if k!='deltas'},indent=2),flush=True)
