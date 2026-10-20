"""Binary64 discovery backend; never used to certify a moment or cutoff."""
import math
from functools import lru_cache
import numpy as np
from minimum_modulus_experiment import factors_table


@lru_cache(maxsize=8)
def geometry(K,cap):
    limit=max(cap,(K-1)//2);fs=factors_table(limit)
    largest=np.array([0]+[f[-1][0] if f else 1 for f in fs[1:]])
    divisors=[[] for _ in range(cap+1)];mobius=[[] for _ in range(cap+1)]
    for n in range(1,cap+1):
        ds=[1];ms=[(1,1)]
        for p,a in fs[n]:
            old=ds[:]
            for b in range(1,a+1):ds.extend(d*p**b for d in old)
            ms += [(d*p,-s) for d,s in ms[:]]
        divisors[n]=np.array(ds,dtype=int)
        mobius[n]=[(n//d,s) for d,s in ms]
    return limit,fs,largest,divisors,mobius


class FloatRecursiveAtoms:
    def __init__(self,K,cap,density_degree=1):
        assert density_degree in (1,2)
        self.K,self.cap,self.density_degree=K,cap,density_degree
        self.limit,self.fs,self.largest,self.divisors,self.mobius=geometry(K,cap)
        self.B=np.zeros(self.limit+1);self.B[1]=1.
        self.F=np.zeros(cap+1);self.F[1]=1.;self.nu={};self.C=1.;self.T=1.;self.last=1

    def coefficients(self,p,limit):
        ns=np.arange(1,limit+1);out=np.zeros(limit+1);out[1:]=1/(p-1)
        power=p
        while power<self.K:
            out[1:][ns*power<self.K]/=p
            power*=p
        return out

    def update(self,p,delta):
        v=1/(1-delta);h=self.coefficients(p,max(self.cap,(self.K-1)//p))
        olds=np.flatnonzero((self.largest[:self.cap+1]<p)&(np.arange(self.cap+1)>0))
        small=np.flatnonzero((self.largest[:(self.K-1)//p+1]<p)&(np.arange((self.K-1)//p+1)>0))
        new={}
        for d in olds[olds*p<=self.cap]:
            value=v*self.F[d]
            if delta:
                H=self.C*self.B[d]/(p-1)
                for q,a in self.fs[d]:H*=(a+1+1/(q-1))/(1+self.nu[q]/(q-1))
                H-=self.B[d]*np.sum((1/(p-1)-h[small])*self.B[small]/self.B[np.gcd(d,small)])
                ns=d*olds//np.gcd(d,olds);mask=ns<=self.cap
                H-=np.sum(h[olds[mask]]*np.maximum(0,self.B[ns[mask]]-self.F[ns[mask]]))
                value=min(value,self.F[d]+v*max(0,H))
                if self.density_degree==2:
                    H2=self.T*self.B[d]/(p-1)**2
                    for q,a in self.fs[d]:
                        H2*=((a+1)**2+(2*a+3)/(q-1)+2/(q-1)**2)/(1+self.nu[q]*(3*q-1)/(q-1)**2)
                    value=min(value,self.F[d]+max(0,H)+v*H2)
            power=p
            while d*power<=self.cap:new[d*power]=value/power;power*=p
        for n in range(p,self.limit+1,p):
            if self.largest[n]!=p:continue
            d=n
            while d%p==0:d//=p
            self.B[n]=self.B[d]*v*d/n
        for n,value in new.items():self.F[n]=min(value,self.B[n])
        self.nu[p]=v;self.C*=1+v/(p-1);self.last=p
        self.T*=1+v*(3*p-1)/(p-1)**2

    def gains(self,p):
        h=self.coefficients(p,self.cap);z=np.zeros(self.cap+1);ans=[0.,0.]
        for n in range(1,self.cap+1):
            if self.largest[n]>self.last:continue
            z[n]=np.sum(h[self.divisors[n]])
            diff=max(0,self.B[n]-self.F[n])
            ans[0]+=diff*h[n]
            c=sum(s*z[d]**2 for d,s in self.mobius[n])
            ans[1]+=diff*max(0,c)
        return ans

    def loss(self,p,delta,m1,m2):
        a,b=self.gains(p)
        return min(max(0,m1-a),max(0,m2-b)/(4*delta*(1-delta)) if delta else math.inf)
