"""Independent exact-rational checks of the finite moment reduction."""
from fractions import Fraction as F
from functools import lru_cache
from math import prod
from verify_minimum_modulus import I,S,moments,factorizations,primes_to


def enclosed(q,x):
    assert F(x.lo,S)<=q<=F(x.hi,S),(q,x.lo,x.hi)


def categories(s,p):
    out=[(0,0,s)];j=1;power=p
    while power<s:
        out.append((j,j,(s-1)//power+1));power*=p;j+=1
    out.append((j,None,1))
    return out


def weight(a,b,p):
    j,jend,_=a;k,kend,_=b;r=F(1,p)
    if jend is not None and kend is not None:return r**max(j,k)
    if jend is None and kend is not None:return weight(b,a,p)
    if jend is not None:
        if j<k:return r**k/(1-r)
        return (j-k+1)*r**j+r**(j+1)/(1-r)
    low,high=sorted((j,k))
    return r**high*((high-low+1)/(1-r)+2*r/(1-r)**2)


def direct_theta(primes,nu):
    @lru_cache(None)
    def theta(i,s,t):
        if i==0:return F(int(s<=1 and t<=1))
        p=primes[i-1];v=nu[p]
        return sum((1 if a[0]==b[0]==0 else v)*weight(a,b,p)*theta(i-1,a[2],b[2])
                   for a in categories(s,p) for b in categories(t,p))
    return theta


def check():
    arithmetic=0
    for a in range(-7,8):
        for b in range(1,8):
            for c in range(-3,4):
                q=F(a,b);r=F(c,5);x=I.rat(a,b);y=I.rat(c,5)
                enclosed(q+r,x+y);enclosed(q-r,x-y);enclosed(q*r,x*y)
                arithmetic+=3
                if c>0:enclosed(q/r,x/y);arithmetic+=1
    cases=0;fs=factorizations(100)
    ps=primes_to(20)
    for count in (1,2,3,4):
        prev=ps[:count];p=ps[count]
        nu={q:F(10+j,9+j) for j,q in enumerate(prev)}
        C=prod(1+nu[q]/(q-1) for q in prev)
        T=prod(1+nu[q]*F(3*q-1,(q-1)**2) for q in prev)
        theta=direct_theta(prev,nu)
        for K in (p,p+1,p*p,p*p+1,37,83):
            # The outer p-adic sum has coefficient p^(-j-k), unlike Theta.
            cats=categories(K,p)[1:]
            zs=[(a[2], F(1,p**a[0]) if a[1] is not None else F(1,p**(a[0]-1)*(p-1))) for a in cats]
            exact=sum(z*w*theta(count,s,t) for s,z in zs for t,w in zs)
            _,bound=moments(K,p,I.rat(C.numerator,C.denominator),I.rat(T.numerator,T.denominator),
                            {q:I.rat(v.numerator,v.denominator) for q,v in nu.items()},fs)
            enclosed(exact,bound);cases+=1
    print(f'PASS: {arithmetic} interval arithmetic checks; {cases} independent rational moment checks.')


if __name__=='__main__':check()
