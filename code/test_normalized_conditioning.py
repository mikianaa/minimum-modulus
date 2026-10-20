"""Independent finite-product checks for conditioning on a down-set."""
from fractions import Fraction as F
from itertools import product
from math import prod
from verify_sparse_conditioning import SparseSmallBound,LowerProductBound
from verify_minimum_modulus import S
from verify_stoploss import L


def check():
    primes=[2,3,5];states=[]
    for es in product(range(7),repeat=3):
        n=prod(q**e for q,e in zip(primes,es));tau=prod(e+1 for e in es)
        weight=prod(F(1,q**6) if e==6 else F(q-1,q**(e+1)) for q,e in zip(primes,es))
        states.append((n,tau,weight))
    assert sum(w for n,t,w in states)==1
    cases=0;rest=LowerProductBound(100,8);p=13
    for K,N,dd in product((100,500,2000),(10,20,40),(L//10,L//3,L//2)):
        A=F(dd*(p-1),L);full=F(0);retained=F(0);mass=F(0)
        for n,tau,w in states:
            B=F(0);power=p
            while power<K:
                B+=F(p-1,power)*sum(n%d==0 for d in range(1,(K-1)//power+1));power*=p
            gain=min(B,max(F(0),tau-A));full+=w*gain
            if n<N:retained+=w*gain;mass+=w
        assert retained/mass<=full
        sm=SparseSmallBound(primes,N,(K-1)//p+1)
        raw=sm.correction_lower(K,p,dd,rest)
        normalized=F(raw*S//sm.probability_upper,S)
        assert normalized<=retained/mass<=full
        cases+=1
    print('PASS:',cases,'exact normalized-conditioning comparisons.')


if __name__=='__main__':check()
