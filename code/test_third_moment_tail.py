"""Independent exact-factor and high-precision tail-enclosure checks."""
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
from math import factorial
from third_moment_tail import third_factor,tail_bound
from verify_minimum_modulus import S,ONE,primes_to
from verify_stoploss import L


def check():
    count=0
    for p in (2,3,5,13,59,1009):
        for dd in (0,L//7,L//3,L//2):
            r=F(1,p);nu=F(L,L-dd)
            # Independent ordinary generating function for sum h^3 r^h.
            total=r*(1+4*r+r*r)/(1-r)**4
            exact=1-nu/p+nu*(p-1)*(total-r)
            got=third_factor(p,dd)
            assert F(got.lo,S)<=exact<=F(got.hi,S)
            count+=1
    tails=0
    with localcontext() as ctx:
        ctx.prec=75
        for P in (100003,611953,1299709):
            n=len(primes_to(P));ell=D(P).ln()
            for dd in (L//3,450000000000,L//2):
                bound,info=tail_bound(P,n,ONE,dd)
                c=D(info['prime_count_constant']);m=info['growth_exponent']
                d=D(dd)/L;r=D(P)/(P-1);a=7*r**3/(1-d)
                z=a*(c/ell-D(n)/P)
                def integral(j):
                    return sum(D(factorial(j))/factorial(j-h)/(2*ell)**h for h in range(j+1))/(2*D(P)**2)
                bracket=c*(3*integral(m-1)/ell-m*integral(m-2)/ell**2)-D(n)/P**3
                exact=4/(27*d*d*(1-d))*r**3*z.exp()*bracket
                assert D(bound)/S>=exact>0
                assert (D(bound)/S)/exact<D('1.001')
                tails+=1
    print('PASS:',count,'exact third-moment factors;',tails,'75-digit tail comparisons.')
    for P in (100003,611953,1299709):
        n=len(primes_to(P))
        for dd in (L//3,450000000000,L//2):
            old,_=tail_bound(P,n,ONE,dd)
            new,info=tail_bound(P,n,ONE,dd,True)
            assert 0<new<=old
            # Convexity in the exponent is the new proof step. Verify its
            # numerical enclosure independently at 75 digits.
            with localcontext() as ctx:
                ctx.prec=75;ell=D(P).ln();c=D(info['prime_count_constant']);m=info['growth_exponent']
                d=D(dd)/L;r=D(P)/(P-1);a=7*r**3/(1-d);theta=D(info['growth_exponent_upper'])-(m-1)
                z=a*(c/ell-D(n)/P)
                def jj(j):
                    return sum(D(factorial(j))/factorial(j-h)/(2*ell)**h for h in range(j+1))/(2*D(P)**2)
                def prime_sum(j):return c*(3*jj(j-1)/ell-j*jj(j-2)/ell**2)-D(n)/P**3
                exact=4/(27*d*d*(1-d))*r**3*z.exp()*((1-theta)*prime_sum(m-1)+theta*prime_sum(m))
                assert D(new)/S>=exact>0
    print('PASS: 9 interpolated tail enclosures and monotone improvements.')
    shifted=0
    for P in (100003,611953,1299709,2015177):
        n=len(primes_to(P))
        for dd in (L//3,450000000000,L//2):
            bound,info=tail_bound(P,n,ONE,dd,True,True)
            with localcontext() as ctx:
                ctx.prec=75;ell=D(P).ln()-D('1.1');m=info['growth_exponent']
                d=D(dd)/L;r=D(P)/(P-1);a=7*r**3/(1-d)
                theta=D(info['growth_exponent_upper'])-(m-1);z=a*(1/ell-D(n)/P)
                def jj(j):
                    return sum(D(factorial(j))/factorial(j-h)/(2*ell)**h for h in range(j+1))/(2*D(P)**2)
                def ps(j):return 3*jj(j-1)/ell-j*jj(j-2)/ell**2-D(n)/P**3
                exact=4/(27*d*d*(1-d))*r**3*z.exp()*((1-theta)*ps(m-1)+theta*ps(m))
                assert D(bound)/S>=exact>0
                assert (D(bound)/S)/exact<D('1.001')
            shifted+=1
    print('PASS:',shifted,'independent shifted-log tail enclosures.')


if __name__=='__main__':check()
