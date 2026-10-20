"""Finite support of the exact early restricted-comparison complement.

If n >= L0**2, divisors below L0 occupy at most half its divisor pairs.
Thus W=tau(n)-B(n) >= tau(n)/2. Only n<L0**2 or tau(n)<=2A can
contribute to (A-W)+. Both sets have finite exact enumerations.
"""
import bisect


def specification(K,p):
    powers=[];pw=p
    while pw<K:powers.append(pw);pw*=p
    assert powers
    thresholds=[(K-1)//pw+1 for pw in powers]
    return thresholds,powers[-1],[(p-1)*(powers[-1]//pw) for pw in powers]


def states(K,p,primes,T):
    thresholds,bden,multipliers=specification(K,p);L0=thresholds[0];N=L0*L0
    def finish(n,tau,divisors,mask):
        div=sorted(divisors)
        w=tau*bden-sum(m*bisect.bisect_left(div,s) for m,s in zip(multipliers,thresholds))
        assert w>=0
        if 2*w<(T+1)*bden:yield n,mask,w
    def small(j,n,tau,divisors,mask):
        if j==len(primes):yield from finish(n,tau,divisors,mask);return
        q=primes[j];qp=1;a=0;div=[]
        while n*qp<N:
            div.extend(d*qp for d in divisors if d*qp<L0)
            yield from small(j+1,n*qp,tau*(a+1),div.copy(),mask|((1<<j) if a else 0))
            qp*=q;a+=1
    def few(j,n,tau,divisors,mask):
        if j==len(primes):
            if n>=N:yield from finish(n,tau,divisors,mask)
            return
        q=primes[j];qp=1;div=[]
        for a in range(T//tau):
            div.extend(d*qp for d in divisors if d*qp<L0)
            yield from few(j+1,n*qp,tau*(a+1),div.copy(),mask|((1<<j) if a else 0))
            qp*=q
    yield from small(0,1,1,[1],0)
    if T:yield from few(0,1,1,[1],0)


def loss_upper(K,p,D,primes,parameters,M1):
    from verify_minimum_modulus import I,ONE,S,ceildiv
    from verify_stoploss import L
    T=2*(p-1)*D//L;_,bden,_=specification(K,p)
    norm=ONE;factors=[]
    for q in primes:
        nu=I.rat(L,L-parameters.get(q,0))
        assert nu.hi<q*S
        norm=norm*(ONE-nu*I.rat(1,q))
        factors.append(nu*I.rat(q-1)/(I.rat(q)-nu))
    cache={0:norm};total=0;count=0
    def factor(mask):
        if mask not in cache:
            bit=mask&-mask;cache[mask]=factor(mask-bit)*factors[bit.bit_length()-1]
        return cache[mask]
    for n,mask,w in states(K,p,primes,T):
        count+=1;gain=max(0,D*(p-1)*bden-w*L)
        prob=ceildiv(factor(mask).hi,n)
        total+=ceildiv(prob*gain,L*bden)
    d=I.rat(D,L)
    result=(M1-d+I.rat(total,S*(p-1)))/(ONE-d)
    assert result.hi>=0
    return result.hi,count
