"""Rigorous BBMST certificate checker. Python standard library only.

All numerical bounds use outward-rounded integer fixed-point intervals.
The supplied deltas are rational witnesses, not trusted numerical results.
No floating-point computation participates in the verification.
"""
import argparse
import gzip
import hashlib
import json
import math
import time
from pathlib import Path

S = 10**40


def ceildiv(a,b):
    return -((-a)//b)


class I:
    __slots__ = ('lo','hi')

    def __init__(self,lo,hi=None):
        self.lo=lo
        self.hi=lo if hi is None else hi
        assert self.lo<=self.hi

    @staticmethod
    def rat(a,b=1):
        assert b>0
        return I(a*S//b,ceildiv(a*S,b))

    def __add__(self,o):
        return I(self.lo+o.lo,self.hi+o.hi)

    def __sub__(self,o):
        return I(self.lo-o.hi,self.hi-o.lo)

    def __mul__(self,o):
        v=(self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi)
        return I(min(v)//S,ceildiv(max(v),S))

    def __truediv__(self,o):
        assert o.lo>0
        v=(self.lo*S//o.lo,self.lo*S//o.hi)
        w=(ceildiv(self.hi*S,o.lo),ceildiv(self.hi*S,o.hi))
        return I(min(v),max(w))


ZERO=I(0)
ONE=I(S)
TWO=I(2*S)
FOUR=I(4*S)


def isum(values):
    lo=hi=0
    for v in values:
        lo+=v.lo;hi+=v.hi
    return I(lo,hi)


def primes_to(n):
    """Independent elementary sieve for the certificate checker."""
    marked=bytearray(n+1)
    for p in range(2,math.isqrt(n)+1):
        if not marked[p]:
            marked[p*p:n+1:p]=b'\x01'*((n-p*p)//p+1)
    return [p for p in range(2,n+1) if not marked[p]]


def factorizations(n):
    smallest=list(range(n+1))
    for p in range(2,math.isqrt(n)+1):
        if smallest[p]==p:
            for m in range(p*p,n+1,p):
                if smallest[m]==m:smallest[m]=p
    result=[[] for _ in range(n+1)]
    for m in range(2,n+1):
        q=m
        while q>1:
            p=smallest[q];a=0
            while q%p==0:q//=p;a+=1
            result[m].append((p,a))
    return result


def initial(K,primes,start):
    # The forward recurrence generates exactly the p_start-smooth numbers.
    smooth=bytearray(K);smooth[1]=1
    for p in primes[:start]:
        for m in range(p,K,p):
            if smooth[m//p]:smooth[m]=1
    finite=isum(I.rat(1,m) for m in range(1,K) if smooth[m])
    C=T=ONE
    for p in primes[:start]:
        C=C*I.rat(p,p-1)
        T=T*(ONE+I.rat(3*p-1,(p-1)**2))
    mu=(ONE-C+finite).lo
    assert mu>0
    return mu,C,T


def moments(K,p,C,T,nu,factors):
    thresholds=[];power=p
    while power<K:
        thresholds.append(((K-1)//power+1,I.rat(1,power)))
        power*=p
    n=thresholds[0][0]-1 if thresholds else 0
    w=[ZERO]*(n+1);r=[ZERO]*(n+1);J=[ZERO]*(n+1)
    if n:w[1]=r[1]=J[1]=ONE
    for m in range(2,n+1):
        fs=factors[m]
        if fs[-1][0]>=p:continue
        wm=rm=jm=ONE
        for q,a in fs:
            v=nu.get(q,ONE);qa=q**a
            wm=wm*v*I.rat(1,qa)
            rm=rm*v*I.rat(1,qa)*(I.rat(a+1)+I.rat(1,q-1))/(ONE+v*I.rat(1,q-1))
            jm=jm*((I.rat(q)/v-ONE) if a==1 else I.rat(qa-qa//q)/v)
        assert jm.lo>=0
        w[m],r[m],J[m]=wm,rm,jm
    U=isum(z*isum(r[1:s]) for s,z in thresholds)*C
    v=[ZERO]*(n+1)
    for s,z in thresholds:
        for m in range(1,s):v[m]=v[m]+z*w[m]
    parts=[]
    for d in range(1,n+1):
        v_sum=isum(v[d::d])
        parts.append(J[d]*v_sum*v_sum)
    V=isum(parts)
    M2=T*I.rat(1,(p-1)**2)-TWO*U*I.rat(1,p-1)+V
    M1=C*I.rat(1,p-1)-isum(z*isum(w[1:s]) for s,z in thresholds)
    assert M1.lo>=0 and M2.lo>=0
    return M1,M2


def log_unit(x):
    """log(x), 1<=x<=2, from the atanh series with rigorous tail."""
    assert S<=x.lo<=x.hi<=2*S
    z=(x-ONE)/(x+ONE);z2=z*z
    power=z;terms=[]
    for j in range(80):
        terms.append(power*I.rat(2,2*j+1))
        power=power*z2
    # power is now z**161. Tail <= 2*z**161/(161*(1-z*z)).
    partial=isum(terms)
    tail=TWO*power/(I.rat(161)*(ONE-z2))
    return I(partial.lo,partial.hi+tail.hi)


def log_rational(a,b=1):
    assert a>=b>0
    e=0
    while a>=2*b:b*=2;e+=1
    return I.rat(e)*log_unit(TWO)+log_unit(I.rat(a,b))


def stopping_bound(k):
    L=log_rational(k)
    # monotonicity permits using the certified lower endpoint here.
    LL=log_rational(L.lo,S)
    lam=I(L.lo)+I(LL.lo)-I.rat(3)
    assert lam.lo>0
    return (I.rat(k)*lam*lam).lo


def fixed_text(v):
    sign='-' if v<0 else ''
    v=abs(v)
    return f'{sign}{v//S}.{v%S:040d}'


def verify(witness_path,output):
    raw=Path(witness_path).read_bytes()
    if str(witness_path).endswith('.gz'):raw=gzip.decompress(raw)
    witness=json.loads(raw)
    K=witness['K'];start=witness['start'];ds=witness['deltas']
    assert isinstance(K,int) and 2<K<616000
    assert isinstance(start,int) and start>0
    assert ds and all(isinstance(d,int) and 0<d<=5*10**11 for d in ds)
    k=start+len(ds)
    # A guessed sieve range is enlarged until enough primes are available;
    # there is no assumed unverified nth-prime bound.
    limit=max(K,20*k)
    primes=primes_to(limit)
    while len(primes)<k:
        limit*=2;primes=primes_to(limit)
    factors=factorizations((K-1)//primes[start]+1)
    mu,C,T=initial(K,primes,start)
    initial_mu=mu
    nu={};checkpoints=[];first_moment_wins=0
    for idx,D in enumerate(ds,start):
        p=primes[idx];i=idx+1;delta=I.rat(D,10**12)
        assert delta.lo==delta.hi
        if p<K:M1,M2=moments(K,p,C,T,nu,factors)
        else:M1,M2=C*I.rat(1,p-1),T*I.rat(1,(p-1)**2)
        second=M2/(FOUR*delta*(ONE-delta))
        loss=min(M1.hi,second.hi)
        first_moment_wins+=M1.hi<second.hi
        mu-=loss
        assert mu>0, ('nonpositive survivor lower bound',i)
        v=ONE/(ONE-delta)
        C=C*(ONE+v*I.rat(1,p-1))
        T=T*(ONE+v*I.rat(3*p-1,(p-1)**2))
        nu[p]=v
        if i in (100,1000,10000,51000,100000) or i==k:
            fu=ceildiv(T.hi*S,mu)
            checkpoints.append(dict(i=i,p=p,mu_lower=fixed_text(mu),T_upper=fixed_text(T.hi),f_upper=fixed_text(fu)))
    f_upper=ceildiv(T.hi*S,mu)
    stop_lower=stopping_bound(k)
    assert k>=10 and f_upper<stop_lower, ('stopping criterion failed',k)
    result=dict(status='VERIFIED',cutoff=K,minimum_modulus_at_most=K-1,
                arithmetic='integer fixed-point outward intervals',decimal_places=40,
                witness_sha256=hashlib.sha256(raw).hexdigest(),
                initial_mu_lower=fixed_text(initial_mu),
                terminal_index=k,terminal_prime=primes[k-1],
                f_upper=fixed_text(f_upper),stopping_lower=fixed_text(stop_lower),
                certified_gap=fixed_text(stop_lower-f_upper),
                first_moment_wins=first_moment_wins,checkpoints=checkpoints)
    Path(output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('witness',nargs='?',default='data/minimum-modulus-certificate.json.gz')
    ap.add_argument('--output',default='data/minimum-modulus-certificate-verified.json')
    args=ap.parse_args();t=time.time()
    print(json.dumps(verify(args.witness,args.output),indent=2))
    print('seconds',time.time()-t)
