"""Independent rational checks for arbitrary-prime conditioning and grouping."""
from fractions import Fraction as F
from itertools import product
from math import prod
from verify_strong_conditioning import StrongSmallBound
from verify_stoploss import ProductBound,L
from verify_minimum_modulus import S,ceildiv
from test_stoploss import exact_product,exact_excess


def check():
    cases=0
    for primes,es,rest in (([2,3],[3,2],[5,7]),([2,3,5],[3,2,2],[7,11]),
                           ([2,3,5,7],[2,2,1,1],[11])):
        dist,mean=exact_product(rest,[L//4]*len(rest),80)
        b=ProductBound(80,1000)
        for q in rest:b.update(q,L//4)
        mean_upper=ceildiv(mean.numerator*S,mean.denominator)
        total_small=prod(F(q,q-1) for q in primes)
        for K in (14,40,100,2000):
            p=13;sm=StrongSmallBound(primes,es,(K-1)//p+1)
            for D in (0,L//10,L//3,L//2):
                total=F(0);keptmean=F(0)
                for exponents in product(*(range(e+1) for e in es)):
                    tau=prod(e+1 for e in exponents)
                    probability=prod(F(q-1,q**(e+1)) for q,e in zip(primes,exponents))
                    keptmean+=tau*probability
                    divisors=[prod(q**a for q,a in zip(primes,ee))
                              for ee in product(*(range(e+1) for e in exponents))]
                    excluded=F(0);power=p
                    while power<K:
                        excluded+=F(p-1,power)*sum(m*power<K for m in divisors)
                        power*=p
                    threshold=(F((p-1)*D,L)+excluded)/tau
                    total+=probability*tau*exact_excess(threshold,dist,mean)
                total+=(total_small-keptmean)*mean
                exact=total/((p-1)*(1-F(D,L)))
                got=F(sm.loss_upper(K,p,D,mean_upper,b),S)
                assert got>=exact,(primes,K,D,got-exact)
                assert got-exact<F(1,10**25)
                cases+=1
    print(f'PASS: {cases} independent exact-rational generalized conditioning checks.')


if __name__=='__main__':check()
