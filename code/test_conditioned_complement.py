"""Independent Fraction checks of retained and omitted complement bounds."""
from fractions import Fraction as F
from itertools import product
from math import prod
from conditioned_complement_bound import loss_upper
from verify_sparse_conditioning import DynamicSmallBound
from verify_minimum_modulus import I,S
from verify_stoploss import L


def law(q,a,params):
    nu=F(L,L-params.get(q,0))
    return 1-nu/q if not a else nu*F(q-1,q**(a+1))


if __name__=='__main__':
    count=0
    for K,D,tilt in product((70,100,200),(0,L//7,2*L//5),(False,True)):
        p=11;ps=[2,3,5,7];N=400
        params={2:L//5,3:L//4,5:L//3,7:L//5} if tilt else {}
        sm=DynamicSmallBound(ps[:2],N,(K-1)//p+1)
        for q in ps[:2]:sm.update(q,params.get(q,0))
        A=F((p-1)*D,L);C=prod(1+F(L,(q-1)*(L-params.get(q,0))) for q in ps)
        # Direct finite sum for the exact first moment, by factoring every
        # forbidden divisor. The unrestricted tau mean is the Euler product.
        EB=F(0);pw=p
        while pw<K:
            for d in range(1,(K-1)//pw+1):
                n=d;prob=F(1)
                for q in ps:
                    a=0
                    while n%q==0:n//=q;a+=1
                    if a:prob*=F(L,(L-params.get(q,0))*q**a)
                if n==1:EB+=F(p-1,pw)*prob
            pw*=p
        m1=(C-EB)/(p-1)
        # The pairing lemma confines positive complements to n<L0^2 or
        # tau<=2A. These exponent ranges contain both finite sets.
        L0=(K-1)//p+1;ranges=[]
        for q in ps:
            a=0;power=1
            while power<L0*L0:power*=q;a+=1
            ranges.append(range(max(a,int(2*A))))
        exact=F(0);outside_exact=F(0)
        for aa in product(*ranges):
            n=prod(q**a for q,a in zip(ps,aa));tau=prod(a+1 for a in aa)
            divisors=[1]
            for q,a in zip(ps,aa):divisors=[d*q**b for d in divisors for b in range(a+1)]
            B=F(0);pw=p
            while pw<K:B+=F(p-1,pw)*sum(d*pw<K for d in divisors);pw*=p
            complement=max(F(0),A+B-tau)
            prob=prod(law(q,a,params) for q,a in zip(ps,aa))
            exact+=prob*complement
            if prod(q**a for q,a in zip(ps[:2],aa[:2]))>=N:outside_exact+=prob*complement
        expected=(m1-F(D,L)+exact/(p-1))/(1-F(D,L))
        got,outside=loss_upper(sm,K,p,D,ps[2:],params,5,N,I.rat(m1.numerator,m1.denominator))
        assert F(got,S)>=expected,(K,D,tilt,got,expected)
        assert F(outside,S)>=outside_exact
        count+=1
    print('PASS:',count,'independent rational conditioned-complement enclosures, including omitted configurations.')
