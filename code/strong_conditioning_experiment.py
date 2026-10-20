"""Stronger conditioning and parameter-policy experiments. NOT a certificate."""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from minimum_modulus_experiment import sieve,factors_table,initial,moments
from binned_distribution import Binned,MomentPrefixes


class Conditioning:
    def __init__(self,K,p,primes,exponents):
        self.primes=primes;self.exponents=exponents
        states=np.indices(tuple(e+1 for e in exponents)).reshape(len(primes),-1).T
        self.states=states
        self.tau=np.prod(states+1,axis=1)
        self.weights=self.tau.astype(float)
        self.mean=math.prod(q/(q-1) for q in primes)
        for j,q in enumerate(primes):self.weights*=(q-1)*float(q)**(-states[:,j]-1)
        self.tail=max(0.,self.mean-math.fsum(self.weights))
        n=(K-1)//p+1
        self.count=np.zeros((len(states),n+1),dtype=np.int16)
        def visit(j,m,mask):
            if j==len(primes):
                self.count[:,m+1]+=mask
                return
            q=primes[j];a=0
            while m<n and a<=exponents[j]:
                visit(j+1,m,mask&(states[:,j]>=a));a+=1;m*=q
        visit(0,1,np.ones(len(states),dtype=bool))
        self.count=np.cumsum(self.count,axis=1,dtype=np.int16)
        self.cache={}

    def association_data(self,K,p,nu):
        """Experimental positive-association correction for omitted divisors."""
        base=np.zeros(len(self.tau));maximum=base.copy();mean=base.copy()
        power=p
        while power<K:
            s=(K-1)//power+1;mult=(p-1)/power
            base+=mult*self.count[:,s]
            # Sum over all forbidden divisors, separating the small part.
            fs=factors_table(s-1)
            grouped={}
            for m in range(1,s):
                if fs[m] and fs[m][-1][0]>=p:continue
                es=[0]*len(self.primes);prob=1.;has_rest=False
                for q,a in fs[m]:
                    if q in self.primes:es[self.primes.index(q)]=a
                    else:prob*=nu.get(q,1.)/q**a;has_rest=True
                if not has_rest:continue
                key=tuple(es)
                count,weight=grouped.get(key,(0,0.))
                grouped[key]=(count+1,weight+prob)
            for es,(count,weight) in grouped.items():
                mask=np.all(self.states>=np.asarray(es),axis=1)
                maximum+=mult*count*mask;mean+=mult*weight*mask
            power*=p
        return self.tau,base,self.weights,maximum,mean

    def grouped(self,K,p):
        powers=[];power=p
        while power<K:powers.append(power);power*=p
        ss=tuple((K-1)//pw+1 for pw in powers)
        if ss not in self.cache:
            mat=np.column_stack([self.tau]+[self.count[:,s] for s in ss])
            keys,inverse=np.unique(mat,axis=0,return_inverse=True)
            weights=np.bincount(inverse,weights=self.weights)
            self.cache[ss]=(keys,weights)
        keys,weights=self.cache[ss]
        B=np.zeros(len(keys))
        for j,pw in enumerate(powers):B+=(p-1)/pw*keys[:,j+1]
        return keys[:,0],B,weights


def run(K=10000,start=12,cap=100000,resolution=128,limit=1000000,
        small_exponents=(26,17,11,9),scale=1.,policy=1.,min_stop=0,record=False,
        max_steps=0,switch=False,association=False,moment_bins=False):
    ps=sieve(limit);small_primes=ps[:len(small_exponents)]
    assert start>=len(small_exponents)
    sm=Conditioning(K,ps[start],small_primes,small_exponents)
    fs=factors_table((K-1)//ps[start]+1)
    mu,C,T=initial(K,ps,start);initial_mu=mu;nu={}
    binned=Binned(cap,resolution);values=binned.values
    dist=np.zeros(len(values));dist[1]=1.;moment=dist.copy()
    for q in ps[len(small_primes):start]:
        if moment_bins:dist,moment=binned.update_pair(dist,moment,q,1.)
        else:dist=binned.update(dist,q,1.)
    ds=[];rows=[];wins={'first':0,'second':0,'conditional':0}
    best_ratio=math.inf;switched=False
    for idx,p in enumerate(ps[start:],start):
        use_association=association and p*ps[len(small_primes)]<K
        if switch and p*2>=K and not switched:
            dist=np.zeros(len(values));dist[1]=1.;moment=dist.copy()
            for q in ps[:idx]:
                if moment_bins:dist,moment=binned.update_pair(dist,moment,q,nu.get(q,1.))
                else:dist=binned.update(dist,q,nu.get(q,1.))
            switched=True
        a=(3*p-1)/(p-1)**2
        m1,m2=moments(K,p,C,T,nu,fs) if p<K else (C/(p-1),T/(p-1)**2)
        cdf=np.cumsum(dist);wcdf=np.cumsum(dist*values)
        mp=MomentPrefixes(values,dist,moment,resolution) if moment_bins else None
        if switched:
            tau=np.asarray([1]);B=np.asarray([(p-1)/p if p<K else 0.]);weights=np.asarray([1.])
            crest=C;tail=0.
        else:
            tau,B,weights=sm.grouped(K,p);crest=C/sm.mean;tail=sm.tail
            if use_association:
                tau,B,weights,extra_max,extra_mean=sm.association_data(K,p,nu)
        def loss(d):
            t=np.minimum((B+d*(p-1))/tau,cap)
            j=np.searchsorted(values,t,side='right')-1
            e=crest-t+(mp.complement(t) if mp else t*cdf[j]-wcdf[j])
            if use_association and not switched:
                # Pointwise min(y,x_+) >= y/M min(M,x_+), followed by
                # positive association of increasing functions of the E_q.
                t2=(B+d*(p-1)+extra_max)/tau
                j2=np.searchsorted(values,np.minimum(t2,cap),side='right')-1
                comp=t*cdf[j]-wcdf[j]
                comp2=t2*cdf[j2]-wcdf[j2]
                gain=np.maximum(0.,extra_max+tau*(comp-comp2))
                gain=np.where(t2<=cap,gain,0.)
                factor=np.divide(extra_mean,extra_max,out=np.zeros_like(extra_mean),where=extra_max>0)
                e-=factor*gain/tau
            return (float(np.sum(weights*e))+crest*tail)/((p-1)*(1-d))
        # policy>1 gives greater weight to survivor preservation in the
        # local objective; this is an experimental policy, never a proof input.
        def objective(d,l):
            return math.log1p(a/(1-d))-policy*math.log1p(-l/mu) if 0<=l<mu else math.inf
        def score(d):return objective(d,loss(d))
        lo=0.;hi=.5;g=(5**.5-1)/2;x=hi-g*hi;y=g*hi;fx=score(x);fy=score(y)
        for _ in range(18):
            if fx<fy:hi,y,fy=y,x,fx;x=hi-g*(hi-lo);fx=score(x)
            else:lo,x,fx=x,y,fy;y=lo+g*(hi-lo);fy=score(y)
        d=min(.5,(x if fx<fy else y)*scale)
        choices=[]
        for dd in (0.,d,.5):
            ll=loss(dd)
            if 0<=ll<mu:choices.append((objective(dd,ll),dd,ll,'conditional'))
        if 0<=m1<mu:choices.append((objective(0.,m1),0.,m1,'first'))
        if 0<m2<mu:
            dd=(1+a)/(1+math.sqrt(1+4*mu*a*(1+a)/m2))
            ll=m2/(4*dd*(1-dd))
            if ll<mu:choices.append((objective(dd,ll),dd,ll,'second'))
        if not choices:
            return dict(K=K,start=start,success=False,p=p,mu=mu,rows=rows,best_ratio=best_ratio)
        _,d,ll,kind=min(choices);mu-=ll;wins[kind]+=1
        if record:ds.append(round(d*10**12))
        v=1/(1-d);C*=1+v/(p-1);T*=1+a*v;nu[p]=v
        if moment_bins:dist,moment=binned.update_pair(dist,moment,p,v)
        else:dist=binned.update(dist,p,v)
        i=idx+1;f=T/mu;stop=i*(math.log(i)+math.log(math.log(i))-3)**2
        best_ratio=min(best_ratio,f/stop)
        if i in (100,1000,10000,50000,100000) or (f<stop and i>=min_stop):
            rows.append(dict(i=i,p=p,mu=mu,f=f,stop=stop,ratio=f/stop))
        if (f<stop and i>=min_stop) or (max_steps and i>=max_steps):break
    result=dict(K=K,start=start,cap=cap,resolution=resolution,success=f<stop,
        rows=rows,terminal=dict(i=i,p=p,mu=mu,f=f,stop=stop,ratio=f/stop),
        small_exponents=list(small_exponents),small_primes=small_primes,
        policy=policy,scale=scale,switch=switch,association=association,moment_bins=moment_bins,initial_mu=initial_mu,best_ratio=best_ratio,
        small_tail=sm.tail,wins=wins)
    if record:result['deltas']=ds
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name,default in [('K',10000),('start',12),('cap',100000),('resolution',128),
                         ('limit',1000000),('min_stop',0),('max_steps',0)]:
        ap.add_argument('--'+name.replace('_','-'),type=int,default=default)
    ap.add_argument('--small-exponents',default='26,17,11,9')
    ap.add_argument('--scale',type=float,default=1.);ap.add_argument('--policy',type=float,default=1.)
    ap.add_argument('--association',action='store_true')
    ap.add_argument('--switch',action='store_true');ap.add_argument('--moment-bins',action='store_true')
    ap.add_argument('--output');args=vars(ap.parse_args());out=args.pop('output')
    args['small_exponents']=tuple(map(int,args['small_exponents'].split(',')))
    t=time.time();r=run(**args,record=bool(out));r['seconds']=time.time()-t
    if out:Path(out).write_text(json.dumps(r))
    print(json.dumps({k:v for k,v in r.items() if k!='deltas'},indent=2),flush=True)
