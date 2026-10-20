"""Outward integer certificate for dyadic-order prime moment tails.

Uses the same convex comparison and Dusart input as third_moment_tail.py.
Square roots are enclosed by integer square roots at scale 10**40.
"""
from math import isqrt
from functools import lru_cache
from verify_minimum_modulus import I,ONE,S,ceildiv,log_rational,fixed_text
from verify_stoploss import L
from fractions import Fraction


def root(x):
    assert x.lo>=0
    lo=isqrt(x.lo*S);hi=isqrt(x.hi*S)
    if hi*hi<x.hi*S:hi+=1
    return I(lo,hi)


def quarter_power(x,numerator):
    return dyadic_power(x,numerator,4)


def dyadic_power(x,numerator,denominator=4):
    assert denominator in (4,8,16,32) and numerator>=0
    whole,part=divmod(numerator,denominator);result=ONE
    for _ in range(whole):result=result*x
    if part:
        inside=ONE
        for _ in range(part):inside=inside*x
        for _ in range(denominator.bit_length()-1):inside=root(inside)
        result=result*inside
    return result


@lru_cache(maxsize=None)
def power(a,numerator=10,denominator=4):
    return dyadic_power(I.rat(a),numerator,denominator)


def increment(q,numerator=10,denominator=4):
    assert denominator in (4,8,16,32) and denominator<numerator<=5*denominator
    if numerator==5*denominator:
        return I.rat(31*q**4+56*q**3+36*q*q-4*q+1,(q-1)**5)
    if numerator==4*denominator:
        return I.rat(15*q**3+5*q*q+5*q-1,(q-1)**4)
    # For 1 < s <= 3, (a+1)^s-a^s <= (2^s-1) a(a+1)/2.
    H=1;qp=q
    while qp<10**35:H+=1;qp*=q
    total=I.rat(0);den=q
    for a in range(1,H+1):
        total=total+(power(a+1,numerator,denominator)-power(a,numerator,denominator))*I.rat(1,den);den*=q
    first=I.rat((H+1)*(H+2),2*den)
    ratio=I.rat(H+3,q*(H+1))
    b=power(2,numerator,denominator)-ONE
    # Nonnegative omitted series; only its upper bound is added.
    tail=(b*first/(ONE-ratio)).hi
    if numerator>3*denominator:
        a=H+1
        first=I.rat(4*a**3+6*a*a+4*a+1,den)
        ratio=I.rat((a+1)**3,q*a**3)
        tail=(b*I.rat(1,15)*first/(ONE-ratio)).hi
    if numerator>4*denominator:
        a=H+1
        first=I.rat(5*a**4+10*a**3+10*a*a+5*a+1,den)
        ratio=I.rat((a+1)**4,q*a**4)
        tail=(b*I.rat(1,31)*first/(ONE-ratio)).hi
    return I(total.lo,total.hi+tail)


def fractional_moment(primes,start,ds,numerator=10,denominator=4):
    u=ONE
    for j,q in enumerate(primes[:start+len(ds)]):
        D=0 if j<start else ds[j-start]
        u=u*(ONE+I.rat(L,L-D)*increment(q,numerator,denominator))
    return u


def fractional_tail_bound(P,n,U,D=500000000000,numerator=10,denominator=4):
    assert P>=60184 and n>=1 and 0<D<=L//2
    assert denominator in (4,8,16,32) and denominator<numerator<=5*denominator
    ell=log_rational(P)-I.rat(11,10);d=I.rat(D,L);nu=ONE/(ONE-d)
    r=I.rat(P,P-1);r3=r*r*r;rs=dyadic_power(r,numerator,denominator)
    alpha=nu*(power(2,numerator,denominator)-ONE)*r3
    if numerator>3*denominator:
        alpha=nu*(power(2,numerator,denominator)-ONE)*I.rat(1,15)*(I.rat(15)+I.rat(5,P)+I.rat(5,P*P))*r3*r
    if numerator>4*denominator:
        alpha=nu*(power(2,numerator,denominator)-ONE)*I.rat(1,31)*(I.rat(31)+I.rat(56,P)+I.rat(36,P*P)+I.rat(1,P**4))*r3*r*r
    exponent=max(alpha.hi,2*S);m=max(3,ceildiv(exponent,S))
    order=I.rat(numerator,denominator);beta=I.rat(numerator-denominator,denominator)
    assert m>=3 and (order*ell).lo>m*S
    z=alpha*(ONE/ell-I.rat(n,P));assert 0<=z.lo<=z.hi<S
    Pbeta=dyadic_power(I.rat(P),numerator-denominator,denominator);Porder=dyadic_power(I.rat(P),numerator,denominator)
    def integral(j):
        term=total=ONE
        for h in range(1,j+1):
            term=term*I.rat(j-h+1)/(beta*ell);total=total+term
        return total/(beta*Pbeta)
    def prime_sum(j):
        return order*integral(j-1)/ell-I.rat(j)*integral(j-2)/(ell*ell)-I.rat(n)/Porder
    theta=I(exponent-(m-1)*S);assert 0<=theta.lo<=S
    prime_tail=(ONE-theta)*prime_sum(m-1)+theta*prime_sum(m)
    assert prime_tail.lo>0
    coefficient=dyadic_power(beta,numerator-denominator,denominator)/(dyadic_power(order,numerator,denominator)*dyadic_power(d,numerator-denominator,denominator)*(ONE-d))
    tail=coefficient*U*rs/(ONE-z)*prime_tail
    return tail.hi,dict(fractional_moment_upper=fixed_text(U.hi),moment_order=str(Fraction(numerator,denominator)),moment_numerator=numerator,moment_denominator=denominator,
        tail_delta_numerator=D,shifted_tail=True,interpolate_tail=True,
        growth_exponent=m,growth_exponent_upper=fixed_text(exponent),
        exponential_argument_upper=fixed_text(z.hi),tail_loss_upper=fixed_text(tail.hi))
