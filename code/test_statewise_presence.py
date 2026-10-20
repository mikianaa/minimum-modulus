"""Exact independent statewise-max comparisons on truncated product laws."""
from fractions import Fraction as F
from itertools import product
from math import prod
from verify_minimum_modulus import S
from verify_stoploss import L
from verify_sparse_conditioning import DynamicSmallBound,LowerProductBound
from joint_presence_bound import joint_correction_lower
from test_dynamic_joint import weight


if __name__=='__main__':
    qs=[5,7,11];params={5:L//4,7:L//3,11:L//5};p=13;cases=0
    law={}
    for exponents in product(range(4),repeat=3):
        prob=F(1)
        for q,e in zip(qs,exponents):
            nu=F(L,L-params[q]);prob*=1-nu/q if e==0 else nu*F(q-1,q**(e+1))
        z=prod(e+1 for e in exponents);law[z]=law.get(z,F(0))+prob
    rest=LowerProductBound(1000,256);rest.mass=[0]*len(rest.bins)
    for z,prob in law.items():rest.mass[rest.bins.index(z)]=prob.numerator*S//prob.denominator
    for K,N,D in product((80,150),(20,100),(0,L//4,L//2)):
        small=DynamicSmallBound([2,3],N,(K-1)//p+1)
        small.update(2,0);small.update(3,L//4);exact=F(0);A=F((p-1)*D,L)
        for n in range(1,N):
            wn=weight(n,[2,3],[0,L//4])
            if not wn:continue
            div=[d for d in range(1,n+1) if n%d==0];tau=len(div);B=F(0);pw=p
            while pw<K:B+=F(p-1,pw)*sum(d*pw<K for d in div);pw*=p
            old=sum(prob*min(B,max(F(0),tau*z-A)) for z,prob in law.items())
            new=F(0)
            for bits in product((0,1),repeat=3):
                prob=F(1);extra=0
                for q,b in zip(qs,bits):
                    theta=F(L,q*(L-params[q]));prob*=theta if b else 1-theta
                    extra+=b*sum(d*p*q<K for d in div)
                new+=prob*min(B+F(p-1,p)*extra,max(F(0),tau*2**sum(bits)-A))
            exact+=wn*max(old,new)
        got=F(joint_correction_lower(small,K,p,D,qs,params,5,rest),S)
        assert got<=exact and exact-got<F(1,10**30),(got,exact)
        cases+=1
    print('PASS:',cases,'exact statewise maximum comparisons.')
