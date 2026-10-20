"""Exhaustive independent checks of the finite early-complement reduction."""
from fractions import Fraction as F
from itertools import product
from math import prod
from early_exact_states import states,loss_upper
from verify_minimum_modulus import I,S
from verify_stoploss import L


if __name__=='__main__':
    cases=0
    for K,p,D,tilt in product((20,50,100),(5,7),(0,L//8,3*L//10,L//2),(False,True)):
        ps=[q for q in (2,3,5) if q<p];params={3:L//5,5:L//7} if tilt else {}
        nus={q:F(L,L-params.get(q,0)) for q in ps};d=F(D,L);A=(p-1)*d
        L0=(K-1)//p+1;N=L0**2;T=int(2*A);bounds=[]
        for q in ps:
            a=0;qp=1
            while qp*q<N:a+=1;qp*=q
            bounds.append(max(a,T-1))
        complement=F(0);positive=set()
        for exps in product(*(range(a+1) for a in bounds)):
            n=prod(q**a for q,a in zip(ps,exps));tau=prod(a+1 for a in exps)
            divisors=[1];prob=F(1)
            for q,a in zip(ps,exps):
                divisors=[v*q**j for v in divisors for j in range(a+1)]
                prob*=1-nus[q]/q if a==0 else nus[q]*F(q-1,q**(a+1))
            B=F(0);pw=p
            while pw<K:B+=F(p-1,pw)*sum(v*pw<K for v in divisors);pw*=p
            if A>tau-B:positive.add(n)
            complement+=prob*max(F(0),A-tau+B)
        generated=list(states(K,p,ps,T));ns=[n for n,mask,w in generated]
        assert len(ns)==len(set(ns)) and positive<=set(ns)
        C=prod(1+nus[q]/(q-1) for q in ps);mean=C/(p-1);pw=p
        while pw<K:
            for m in range(1,(K-1)//pw+1):
                rem=m;prob=F(1)
                for q in ps:
                    a=0
                    while rem%q==0:rem//=q;a+=1
                    if a:prob*=nus[q]/q**a
                if rem==1:mean-=prob/pw
            pw*=p
        exact=(mean-d+complement/(p-1))/(1-d)
        bound,_=loss_upper(K,p,D,ps,params,I.rat(mean.numerator,mean.denominator))
        assert F(bound,S)>=exact>=0 and F(bound,S)-exact<F(1,10**30)
        cases+=1
    print('PASS:',cases,'exhaustive early-complement support and loss checks.')
