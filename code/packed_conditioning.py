"""Equivalent floating conditioning grouped by prime-presence mask.

After all small primes are exposed, parameter dependence of a state's
probability is through its presence mask; 1/n is a fixed coefficient.
This representation avoids rescanning every state for each derivative.
"""
import math
import numpy as np


class PackedConditioning:
    def __init__(self,source,primes):
        self.ns=source.ns;self.tau=source.tau;self.count=source.count
        self.last=max(primes);self.primes=list(primes);self.bits={q:1<<j for j,q in enumerate(primes)}
        self.mask=np.zeros(len(self.ns),dtype=np.int64)
        for q,bit in self.bits.items():self.mask|=(self.ns%q==0).astype(np.int64)*bit
        self.base=1./self.ns;self.size=1<<len(primes);self.indices=np.arange(self.size)
        self.cache={};self.compact={};self.weight_cache={};self.reset()

    def reset(self):
        self.factors=np.ones(self.size);self.normalizer=1.;self.weight_cache.clear()

    @property
    def coeff(self):return self.base*self.factors[self.mask]

    def update(self,q,nu):
        if q not in self.bits:return
        self.factors[(self.indices&self.bits[q])!=0]*=nu*(q-1)/(q-nu)
        self.normalizer*=1-nu/q;self.weight_cache.clear()

    def groups(self,K,p):
        assert p>self.last
        powers=[];pw=p
        while pw<K:powers.append(pw);pw*=p
        ss=tuple((K-1)//pw+1 for pw in powers);key=(self.last+1,ss)
        if key not in self.cache:
            keys,inv=np.unique(np.column_stack([self.tau]+[self.count[s] for s in ss]),axis=0,return_inverse=True)
            self.cache[key]=keys,inv
            pairs,iv=np.unique(inv*self.size+self.mask,return_inverse=True)
            self.compact[key]=(pairs//self.size,pairs%self.size,np.bincount(iv,weights=self.base))
        keys,inv=self.cache[key]
        if key not in self.weight_cache:
            groups,masks,base=self.compact[key]
            self.weight_cache[key]=np.bincount(groups,weights=self.normalizer*base*self.factors[masks],minlength=len(keys))
        B=np.zeros(len(keys))
        for j,pw in enumerate(powers):B+=(p-1)/pw*keys[:,j+1]
        return keys[:,0],B,self.weight_cache[key]
