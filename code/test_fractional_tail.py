"""Independent 80-digit moment and prime-tail comparisons."""
from decimal import Decimal as D,localcontext
from math import factorial
from fractional_moment_tail import increment,fractional_tail_bound
from verify_minimum_modulus import S,ONE,primes_to


if __name__=='__main__':
    factors=tails=0
    with localcontext() as ctx:
        ctx.prec=80
        # Include every order in the dense certificate and one finer dyadic.
        for numerator,denominator in [(j,16) for j in range(17,81)]+[(87,32)]:
            s=D(numerator)/denominator;beta=s-1
            for q in (2,3,7,59,1009):
                exact=sum(((D(a+1)**s-D(a)**s)/D(q)**a for a in range(1,501)),D(0))
                got=increment(q,numerator,denominator)
                # Decimal summation can end just below an exact integer
                # (the q=2 fourth increment is exactly 149).
                slack=D('1e-70') if s in (4,5) else D(0)
                assert D(got.lo)/S-slack<=exact<=D(got.hi)/S+slack,(q,numerator,exact,got.lo,got.hi)
                factors+=1
            for P in ((100003,611953,2015177) if s<=4 else (2015177,2750159,10000019)):
                n=len(primes_to(P));ell=D(P).ln()-D('1.1')
                for dd in (300000000000,450000000000,500000000000):
                    bound,r=fractional_tail_bound(P,n,ONE,dd,numerator,denominator)
                    d=D(dd)/10**12;ratio=D(P)/(P-1);alpha=(2**s-1)*ratio**3/(1-d)
                    if s>3:alpha=(2**s-1)/15*(15+D(5)/P+D(5)/D(P)**2)*ratio**4/(1-d)
                    if s>4:alpha=(2**s-1)/31*(31+D(56)/P+D(36)/D(P)**2+D(1)/D(P)**4)*ratio**5/(1-d)
                    z=alpha*(1/ell-D(n)/P);m=r['growth_exponent'];theta=D(r['growth_exponent_upper'])-(m-1)
                    assert D(r['growth_exponent_upper'])>=alpha and 0<=theta<=1
                    def J(j):
                        return sum(D(factorial(j))/factorial(j-h)/(beta*ell)**h for h in range(j+1))/(beta*D(P)**beta)
                    def B(j):return s*J(j-1)/ell-D(j)*J(j-2)/ell**2-D(n)/D(P)**s
                    exact=beta**beta/s**s/d**beta/(1-d)*ratio**s*z.exp()*((1-theta)*B(m-1)+theta*B(m))
                    assert 0<exact<=D(bound)/S and (D(bound)/S)/exact<D('1.001')
                    tails+=1
    print('PASS:',factors,'80-digit fractional moment increments;',tails,'tail enclosures.')
