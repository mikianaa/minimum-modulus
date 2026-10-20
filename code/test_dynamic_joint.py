"""Independent exact probability and joint-correction checks."""
from fractions import Fraction as F
from itertools import product
from math import prod
from verify_minimum_modulus import S
from verify_stoploss import L
from verify_sparse_conditioning import DynamicSmallBound,SparseSmallBound
from joint_presence_bound import joint_correction_lower,presence_distribution_lower


def weight(n,primes,ds):
    result=F(1)
    for q,D in zip(primes,ds):
        a=0
        while n%q==0:n//=q;a+=1
        nu=F(L,L-D)
        result*=1-nu/q if not a else nu*F(q-1,q**(a+1))
    return result if n==1 else F(0)


def check():
    probabilities=0;comparisons=0
    for N in (20,100,500):
        small=DynamicSmallBound([2,3,5],N,100)
        ds=[0,L//4,L//3]
        for j,(q,D) in enumerate(zip([2,3,5],ds)):
            small.update(q,D);exact=sum(weight(n,[2,3,5][:j+1],ds[:j+1]) for n in range(1,N))
            assert F(small.probability_lower,S)<=exact<=F(small.probability_upper,S)
            for n,mask,c,d in zip(small.ns,small.masks,small.coeff_lo,small.coeff_hi):
                if mask>=1<<(j+1):continue
                w=weight(n,[2,3,5][:j+1],ds[:j+1])
                assert F(small.norm_lo*c//S,S)<=w<=F(-(-small.norm_hi*d//S),S)
                probabilities+=1
    outside=[5,7,11];params={5:L//4,7:L//3,11:L//5};p=13
    for K,N,dd,dynamic in product((80,150),(20,100),(0,L//4,L//2),(False,True)):
        small=(DynamicSmallBound if dynamic else SparseSmallBound)([2,3],N,(K-1)//p+1)
        ds=[0,L//4] if dynamic else [0,0]
        if dynamic:
            for q,d in zip([2,3],ds):small.update(q,d)
        exact=F(0);A=F((p-1)*dd,L)
        for n in range(1,N):
            wn=weight(n,[2,3],ds)
            if not wn:continue
            divisors=[d for d in range(1,n+1) if n%d==0];tau=len(divisors)
            B=F(0);pw=p
            while pw<K:
                B+=F(p-1,pw)*sum(d*pw<K for d in divisors);pw*=p
            for bits in product((0,1),repeat=3):
                prob=wn;extra=0
                for q,b in zip(outside,bits):
                    theta=F(L,q*(L-params[q]));prob*=theta if b else 1-theta
                    extra+=b*sum(d*p*q<K for d in divisors)
                exact+=prob*min(B+F(p-1,p)*extra,max(F(0),tau*2**sum(bits)-A))
        got=F(joint_correction_lower(small,K,p,dd,outside,params,5),S)
        assert got<=exact and exact-got<F(1,10**30),(got,exact)
        comparisons+=1
    dist=presence_distribution_lower(outside,params,[1,2,0],count_cap=2)
    for key,mass in dist.items():
        exact=F(0)
        for bits in product((0,1),repeat=3):
            if (min(2,sum(bits)),bits[0]+2*bits[1])!=key:continue
            exact+=prod(F(L,q*(L-params[q])) if b else 1-F(L,q*(L-params[q])) for q,b in zip(outside,bits))
        assert F(mass,S)<=exact
    structural=0
    for K,n,exponents in product((80,150),(1,2,3,6,12,18),product(range(3),repeat=3)):
        small_divisors=[d for d in range(1,n+1) if n%d==0]
        full_divisors=small_divisors
        for q,e in zip(outside,exponents):
            full_divisors=[d*q**a for d in full_divisors for a in range(e+1)]
        full_B=base_B=F(0);power=p
        while power<K:
            full_B+=F(p-1,power)*sum(d*power<K for d in full_divisors)
            base_B+=F(p-1,power)*sum(d*power<K for d in small_divisors)
            power*=p
        extra=sum(sum(d*p*q<K for d in small_divisors) for q,e in zip(outside,exponents) if e)
        assert full_B==base_B+F(p-1,p)*extra
        y_lower=len(small_divisors)*2**sum(e>0 for e in exponents)
        assert y_lower<=len(full_divisors)
        for delta in (F(0),F(1,4),F(1,2)):
            a=(p-1)*delta
            assert min(full_B,max(F(0),y_lower-a))<=min(full_B,max(F(0),len(full_divisors)-a))
        structural+=1
    print('PASS:',probabilities,'dynamic atom enclosures;',comparisons,'exact joint corrections; capped presence law;',structural,'full-divisor identities.')


if __name__=='__main__':check()
