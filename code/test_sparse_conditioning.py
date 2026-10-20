"""Independent Fraction checks of the bounded-correction certificate."""
from fractions import Fraction as F
from math import prod
from verify_sparse_conditioning import SparseSmallBound,LowerProductBound
from moment_product_bound import MomentProductBound
from verify_stoploss import L
from verify_minimum_modulus import S,ceildiv
from test_stoploss import exact_product,exact_excess


def check():
    cases=0
    for small,rest in (([2,3],[5,7]),([2,3,5],[7,11])):
        ds=[L//4]*len(rest)
        rd,rm=exact_product(rest,ds,100)
        fd,fm=exact_product(small+rest,[0]*len(small)+ds,100)
        lb=LowerProductBound(100,8);fb=MomentProductBound(100,8)
        for q in small:fb.update(q,0)
        for q,D in zip(rest,ds):lb.update(q,D);fb.update(q,D)
        for N in (10,50,120):
            independent=[]
            for n in range(1,N):
                remainder=n
                for q in small:
                    while remainder%q==0:remainder//=q
                if remainder==1:
                    divisors=[d for d in range(1,n+1) if n%d==0]
                    independent.append((n,len(divisors),divisors))
            for K in (14,100,2000):
                p=13;sm=SparseSmallBound(small,N,(K-1)//p+1)
                assert len(sm.states)==len(independent)
                for D in (0,L//10,L//3,L//2):
                    A=F(D*(p-1),L);correction=F(0)
                    for n,tau,divisors in independent:
                        prob=prod(F(q-1,q) for q in small)/n
                        B=F(0);power=p
                        while power<K:
                            B+=F(p-1,power)*sum(d*power<K for d in divisors);power*=p
                        gain=tau*(exact_excess(A/tau,rd,rm)-exact_excess((A+B)/tau,rd,rm))
                        correction+=prob*gain
                    lower=F(sm.correction_lower(K,p,D,lb),S)
                    assert 0<=lower<=correction
                    baseline=F(fb.excess_upper(A.numerator*S//A.denominator,
                        ceildiv(fm.numerator*S,fm.denominator),fb.prefixes()),S)
                    assert baseline-lower>=exact_excess(A,fd,fm)-correction
                    cases+=1
    print(f'PASS: {cases} independent exact sparse-correction and combined-loss checks.')


if __name__=='__main__':check()
