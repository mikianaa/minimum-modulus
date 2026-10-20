"""Independent exact rational checks of the joint product complement."""
from fractions import Fraction as F
from itertools import product
from math import prod
from bisect import bisect_left
from joint_product_bound import product_distribution_upper,joint_product_correction_lower
from verify_minimum_modulus import S
from verify_stoploss import L


def law(q,a,parameters):
    v=F(L,L-parameters.get(q,0))
    return 1-v/q if not a else v*F(q-1,q**(a+1))


if __name__=='__main__':
    checks=0
    for H,coeff,tilt in product((3,6,9),((0,1),(1,2)),(False,True)):
        ps=[5,7];params={5:L//3,7:L//5} if tilt else {}
        exact={}
        for aa in product(range(H),repeat=2):
            z=prod(a+1 for a in aa)
            if z>H:continue
            j=sum(c for c,a in zip(coeff,aa) if a)
            exact[z,j]=exact.get((z,j),F(0))+prod(law(q,a,params) for q,a in zip(ps,aa))
        got=product_distribution_upper(ps,params,coeff,H)
        for key,value in exact.items():
            assert F(got[key],S)>=value and F(got[key],S)-value<F(1,10**35)
        checks+=1
    class Small:
        def __init__(self,n):self.div=[d for d in range(1,n+1) if n%d==0]
        def groups(self,thresholds):return [((len(self.div),)+tuple(bisect_left(self.div,s) for s in thresholds),S)]
    for K,n,dd,tilt in product((70,100,200,250),(1,2,6,12),(0,L//7,2*L//5),(False,True)):
        p=11;ps=[5,7];params={5:L//3,7:L//5} if tilt else {}
        sm=Small(n);tau=len(sm.div);A=F((p-1)*dd,L)
        B0=F(0);pw=p
        while pw<K:B0+=F(p-1,pw)*sum(d*pw<K for d in sm.div);pw*=p
        coeff=[sum(d*p*q<K for d in sm.div) for q in ps]
        mean=B0+F(p-1,p)*sum(c*F(L,q*(L-params.get(q,0))) for q,c in zip(ps,coeff))
        H=(A.numerator+A.denominator-1)//A.denominator+1+2
        complement=F(0)
        for aa in product(range(H),repeat=2):
            z=prod(a+1 for a in aa)
            B=B0+F(p-1,p)*sum(c for c,a in zip(coeff,aa) if a)
            h=min(B,max(F(0),A+B-tau*z))
            if z>H:assert not h
            complement+=h*prod(law(q,a,params) for q,a in zip(ps,aa))
        exact=mean-complement
        got=joint_product_correction_lower(sm,K,p,dd,ps,params,5)
        assert 0<=F(got,S)<=exact and exact-F(got,S)<F(1,10**30),(n,dd,exact,got)
        checks+=1
    print('PASS:',checks,'independent rational joint-product distribution/correction checks.')
