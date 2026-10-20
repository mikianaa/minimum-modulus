"""Exploratory joint conditioning on large-prime presence indicators.

For early stages, forbidden divisors contain at most one prime outside the
initial small-prime set. Track its additive divisor contribution jointly
with the number of present large primes. This is a lower bound for the
bounded correction, not a replacement for the exact full-product mean.
"""
import argparse,json,time,math
from pathlib import Path
import numpy as np
from sparse_conditioning_experiment import SparseStates
from minimum_modulus_experiment import sieve
from binned_distribution import Binned


def presence_distribution(primes,nus,coefficients,count_cap=12):
    width=1+sum(coefficients);a=np.zeros((count_cap+1,width));a[0,0]=1.
    used=0
    for q,v,c in zip(primes,nus,coefficients):
        theta=v/q;new=(1-theta)*a
        new[1:,c:c+used+1]+=theta*a[:-1,:used+1]
        new[-1,c:c+used+1]+=theta*a[-1,:used+1]
        a=new;used+=c
    return a


class JointStates:
    def __init__(self,small,K,p0,q0):
        self.sm=small;self.K=K;self.q0=q0;self.cache={}
        # Only stages satisfying the one-remaining-prime identity use this
        # mask. Starting at the first possible such stage avoids irrelevant
        # high bits for smaller explicitly conditioned prime sets.
        while (K-1)//p0>=q0*q0 or p0*p0*q0<K:p0+=1
        self.maxd=(K-1)//(p0*q0)
        assert self.maxd<=32
        self.mask=np.zeros(len(small.ns),dtype=np.int32)
        for d in range(2,self.maxd+1):self.mask|=(small.ns%d==0).astype(np.int32)<<(d-2)

    def groups(self,p):
        powers=[];power=p
        while power<self.K:powers.append(power);power*=p
        thresholds=tuple((self.K-1)//power+1 for power in powers)
        maxd=(self.K-1)//(p*self.q0)
        mask=self.mask&((1<<max(0,maxd-1))-1)
        key=(thresholds,maxd)
        if key not in self.cache:
            keys,inv=np.unique(np.column_stack([mask,self.sm.tau]+[self.sm.count[s] for s in thresholds]),axis=0,return_inverse=True)
            self.cache[key]=keys,inv
        keys,inv=self.cache[key]
        weights=np.bincount(inv,weights=self.sm.prob);B=np.zeros(len(keys))
        for j,pw in enumerate(powers):B+=(p-1)/pw*keys[:,j+2]
        return keys[:,0],keys[:,1],B,weights

    def correction(self,p,d,primes,nus,rest=None,bins=None):
        assert (self.K-1)//p<self.q0**2
        assert p*p*self.q0>=self.K
        if getattr(self,'prepared_p',None)!=p:
            masks,tau,B,weights=self.groups(p);self.terms=[]
            for mask in np.unique(masks):
                ix=masks==mask;coeff=[]
                for q in primes:
                    mx=(self.K-1)//(p*q)
                    coeff.append(int(mx>=1)+sum((int(mask)>>(j-2))&1 for j in range(2,mx+1)))
                dist=presence_distribution(primes,nus,coeff)
                aa,jj=np.nonzero(dist);ww=dist[aa,jj]
                self.terms.append((tau[ix],B[ix],weights[ix],aa,jj,ww))
            self.prepared_p=p
        A=(p-1)*d;total=0.
        for ts,bs,ws,aa,jj,ww in self.terms:
            # Work in chunks to avoid a large group-by-state tensor.
            for lo in range(0,len(ts),256):
                t=ts[lo:lo+256,None];b=bs[lo:lo+256,None]
                gains=np.minimum(b+(p-1)/p*jj,np.maximum(0.,t*2.**aa-A))
                expectation=gains@ww
                if rest is not None:
                    rp=np.cumsum(rest);rw=np.cumsum(rest*bins)
                    tt=ts[lo:lo+256];bb=bs[lo:lo+256]
                    j=np.searchsorted(bins,A/tt,side='right')-1
                    k=np.searchsorted(bins,(A+bb)/tt,side='right')-1
                    old=tt*(rw[k]-rw[j])-A*(rp[k]-rp[j])+bb*(rp[-1]-rp[k])
                    expectation=np.maximum(expectation,old)
                total+=float(np.sum(ws[lo:lo+256]*expectation))
        return total


def gain(w,K=None,state_cap=None):
    K=K or w['K'];ps=sieve(K);start=w['start'];small=w['small'];assert start==small
    sm=SparseStates(K,ps[start],ps[:small],state_cap or w['state_cap'],w.get('normalize_conditioning',False))
    joint=JointStates(sm,K,ps[start],ps[small]);bins=Binned(w['cap'],w['resolution'])
    rest=np.zeros(len(bins.values));rest[1]=1.;pr=[];nus=[];total=0.;rows=[]
    for idx,p in enumerate(ps[start:],start):
        if p*ps[small]>=K:break
        d=w['deltas'][idx-start]/1e12;v=1/(1-d);A=(p-1)*d
        tau,B,weights=sm.groups(K,p);P=np.cumsum(rest);W=np.cumsum(rest*bins.values)
        j=np.searchsorted(bins.values,A/tau,side='right')-1
        k=np.searchsorted(bins.values,(A+B)/tau,side='right')-1
        old=float(np.sum(weights*np.maximum(0.,tau*(W[k]-W[j])-A*(P[k]-P[j])+B*(P[-1]-P[k]))))
        new=joint.correction(p,d,pr,nus);inc=max(0.,new-old)/((p-1)*(1-d));total+=inc
        rows.append(dict(p=p,old=old,new=new,gain=inc,cumulative=total))
        rest=bins.update(rest,p,v);pr.append(p);nus.append(v)
    return dict(K=K,states=len(sm.ns),improvement=total,rows=rows)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness');ap.add_argument('--K',type=int)
    ap.add_argument('--state-cap',type=int);ap.add_argument('--output',required=True)
    args=ap.parse_args();w=json.loads(Path(args.witness).read_text());t=time.time()
    r=gain(w,args.K,args.state_cap);r['seconds']=time.time()-t
    Path(args.output).write_text(json.dumps(r,indent=2));print(json.dumps({k:v for k,v in r.items() if k!='rows'}),flush=True)
