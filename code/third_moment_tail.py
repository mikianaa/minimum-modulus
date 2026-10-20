"""Integer enclosure of an infinite third-moment distortion tail.

External analytic input: Dusart (2010), Theorem 6.9, inequality (6.6):
pi(x) <= x/(log(x)-1.1), x >= 60184.
All evaluated quantities use the existing outward integer interval class.
"""
from verify_minimum_modulus import I,ONE,S,ceildiv,log_rational,fixed_text
from verify_stoploss import L


def third_factor(p,D):
    v=I.rat(L,L-D)
    return ONE+v*I.rat(7*p*p-2*p+1,(p-1)**3)


def third_moment(primes,start,ds):
    u=ONE
    for j,p in enumerate(primes[:start+len(ds)]):
        u=u*third_factor(p,0 if j<start else ds[j-start])
    return u


def tail_bound(P,n,U,D=450000000000,interpolate=False,shifted=False):
    assert P>=60184 and n>=1 and 0<D<=L//2
    lp=log_rational(P);d=I.rat(D,L);v=ONE/(ONE-d)
    # Choose a rational c above log(P)/(log(P)-1.1).
    cu=(lp/(lp-I.rat(11,10))).hi
    c=I.rat(ceildiv(cu*10**6,S),10**6)
    if shifted:
        # Use pi(x) <= x/(log(x)-1.1) directly. In both partial
        # summations the integration variable is log(x)-1.1.
        lp=lp-I.rat(11,10);c=ONE
    r=I.rat(P,P-1);r3=r*r*r;alpha=I.rat(7)*v*r3
    exponent_upper=(alpha*c).hi;m=ceildiv(exponent_upper,S)
    assert m>=2 and 3*lp.lo>m*S
    z=alpha*(c/lp-I.rat(n,P))
    assert 0<=z.lo<=z.hi<S
    exponential_upper=ONE/(ONE-z) # exp(z) <= 1/(1-z), 0<=z<1.
    def integral(r):
        term=total=ONE
        for j in range(1,r+1):
            term=term*I.rat(r-j+1)/(I.rat(2)*lp)
            total=total+term
        return total*I.rat(1,2*P*P)
    def prime_sum(j):
        return c*(I.rat(3)*integral(j-1)/lp-I.rat(j)*integral(j-2)/(lp*lp))-I.rat(n,P**3)
    prime_tail=prime_sum(m)
    if interpolate:
        assert m>=3
        theta=I(exponent_upper-(m-1)*S)
        assert 0<=theta.lo==theta.hi<=S
        prime_tail=(ONE-theta)*prime_sum(m-1)+theta*prime_tail
    assert prime_tail.lo>0
    coefficient=I.rat(4,27)/(d*d*(ONE-d))
    tail=coefficient*U*r3*exponential_upper*prime_tail
    return tail.hi,dict(third_moment_upper=fixed_text(U.hi),tail_delta_numerator=D,
        prime_count_constant=fixed_text(c.hi),growth_exponent=m,
        interpolate_tail=interpolate,shifted_tail=shifted,growth_exponent_upper=fixed_text(exponent_upper),
        exponential_argument_upper=fixed_text(z.hi),tail_loss_upper=fixed_text(tail.hi))
