"""Floating evaluator of the finite exact early complement (not a proof)."""
import math
import numpy as np
from early_exact_states import states,specification


class EarlyExact:
    def __init__(self,K):self.K=K;self.cache={}

    def loss(self,p,d,primes,nus,m1):
        T=math.floor(2*(p-1)*d);key=(p,T)
        if key not in self.cache:
            groups={}
            for n,mask,w in states(self.K,p,primes,T):
                kk=(mask,w);groups[kk]=groups.get(kk,0.)+1./n
            keys=list(groups);_,den,_=specification(self.K,p)
            self.cache[key]=(np.asarray([a for a,b in keys],dtype=np.int64),
                np.asarray([b/den for a,b in keys]),np.asarray([groups[k] for k in keys]))
        masks,values,base=self.cache[key];prob=base.copy()
        for j,q in enumerate(primes):
            v=nus.get(q,1.)
            prob*=np.where((masks&(1<<j))!=0,v*(q-1)/q,1-v/q)
        complement=prob@np.maximum(0.,(p-1)*d-values)
        return (m1-d+complement/(p-1))/(1-d)
