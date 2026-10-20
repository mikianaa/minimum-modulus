"""Exploratory BBMST Sections 6--8 calculation (binary64, NOT a proof)."""
import argparse
import gzip
import json
import math
import time
from pathlib import Path


def sieve(n):
    a = bytearray(b'\x01') * (n + 1)
    a[:2] = b'\x00\x00'
    for p in range(2, math.isqrt(n) + 1):
        if a[p]:
            a[p*p::p] = b'\x00' * ((n-p*p)//p+1)
    return [p for p in range(2, n+1) if a[p]]


def factors_table(n):
    table = [[] for _ in range(n+1)]
    for p in sieve(n):
        for m in range(p,n+1,p):
            q, a = m, 0
            while q % p == 0:
                q //= p
                a += 1
            table[m].append((p,a))
    return table


def initial(K, primes, start):
    assert start>=0
    if start==0:return 1.,1.,1.
    smooth = bytearray(b'\x01')*K
    for p in sieve(K-1):
        if p > primes[start-1]:
            smooth[p::p] = b'\x00'*((K-1)//p)
    total = math.prod(p/(p-1) for p in primes[:start])
    mu = 1-total+math.fsum(1/m for m in range(1,K) if smooth[m])
    return mu, total, math.prod(1+(3*p-1)/(p-1)**2 for p in primes[:start])


def moments(K,p,C,T,nu,factors,first=False):
    """Inclusion/exclusion; finite pair sum via a gcd convolution."""
    thresholds=[]
    power=p
    while power<K:
        thresholds.append(((K-1)//power+1,1/power))
        power*=p
    n=thresholds[0][0]-1 if thresholds else 0
    w=[0.]*(n+1); r=[0.]*(n+1); J=[0.]*(n+1)
    if n:
        w[1]=r[1]=J[1]=1.
    for m in range(2,n+1):
        fs=factors[m]
        if fs[-1][0]>=p:
            continue
        wm,rm,jm=1.,1.,1.
        for q,a in fs:
            v=nu.get(q,1.)
            qa=q**a
            wm*=v/qa
            rm*=v/qa*(a+1+1/(q-1))/(1+v/(q-1))
            jm*=(q/v-1) if a==1 else (qa-qa//q)/v
        w[m],r[m],J[m]=wm,rm,jm
    if first:
        return C/(p-1)-math.fsum(z*math.fsum(w[1:s]) for s,z in thresholds),None
    # Weighted finite one-coordinate complement.
    U=math.fsum(z*math.fsum(r[1:s]) for s,z in thresholds)*C
    # sum_j z_j 1_{m<s_j} enters squared in the gcd convolution.
    v=[0.]*(n+1)
    for s,z in thresholds:
        for m in range(1,s):
            v[m]+=z*w[m]
    V=math.fsum(J[d]*math.fsum(v[d::d])**2 for d in range(1,n+1))
    M2=T/(p-1)**2-2*U/(p-1)+V
    M1=C/(p-1)-math.fsum(z*math.fsum(w[1:s]) for s,z in thresholds)
    return M1,M2


def run(K=616000,start=51,scale=1.,correction=1.,limit=2000000,record=False,beta=.5,ramp=0.,blocks=None,min_stop=0):
    primes=sieve(limit)
    factors=factors_table((K-1)//primes[start]+1)
    mu,C,T=initial(K,primes,start)
    initial_mu=mu
    nu={}
    deltas=[]
    checkpoints=[]
    for idx in range(start,len(primes)):
        p=primes[idx]; i=idx+1
        a=(3*p-1)/(p-1)**2
        M1,M2=moments(K,p,C,T,nu,factors) if p<K else (C/(p-1),T/(p-1)**2)
        if mu<=0 or M2>=mu:
            return dict(K=K,start=start,scale=scale,correction=correction,success=False,i=i,p=p,mu=mu)
        factor=1-correction/(p**beta)
        if ramp:factor*=1-math.exp(-(p-primes[start-1])/ramp)
        if blocks:
            for upper,mult in blocks:
                if p<=upper:
                    factor*=mult
                    break
        d=factor*scale*(1+a)/(1+math.sqrt(1+4*mu*a*(1+a)/M2))
        d=max(1e-12,min(.5,d))
        # Quantize the witness; verification will not need square roots.
        D=round(d*10**12); d=D/10**12
        mu-=min(M1,M2/(4*d*(1-d)))
        C*=1+1/((1-d)*(p-1))
        T*=1+a/(1-d)
        nu[p]=1/(1-d)
        if record: deltas.append(D)
        f=T/mu if mu>0 else math.inf
        stop=i*(math.log(i)+math.log(math.log(i))-3)**2
        if i in (51,100,1000,10000,51000,100000) or (f<stop and i>=min_stop):
            checkpoints.append(dict(i=i,p=p,mu=mu,f=f,stop=stop))
        if f<stop and i>=min_stop:
            ans=dict(K=K,start=start,scale=scale,correction=correction,beta=beta,ramp=ramp,blocks=blocks,success=True,initial_mu=initial_mu,checkpoints=checkpoints)
            if record: ans['deltas']=deltas
            return ans
    ans=dict(K=K,start=start,scale=scale,correction=correction,beta=beta,ramp=ramp,blocks=blocks,success=False,initial_mu=initial_mu,checkpoints=checkpoints)
    if record: ans['deltas']=deltas
    return ans


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--K',type=int,default=616000)
    ap.add_argument('--start',type=int,default=51)
    ap.add_argument('--scale',type=float,default=1.)
    ap.add_argument('--correction',type=float,default=1.)
    ap.add_argument('--beta',type=float,default=.5)
    ap.add_argument('--ramp',type=float,default=0.)
    ap.add_argument('--blocks',help='JSON file containing a blocks field')
    ap.add_argument('--limit',type=int,default=2000000)
    ap.add_argument('--min-stop',type=int,default=0)
    ap.add_argument('--output')
    args=vars(ap.parse_args()); out=args.pop('output')
    if args['blocks']:args['blocks']=json.loads(Path(args['blocks']).read_text())['blocks']
    t=time.time(); result=run(**args,record=bool(out))
    if out:
        raw=json.dumps(result).encode()
        Path(out).write_bytes(gzip.compress(raw,mtime=0) if out.endswith('.gz') else raw)
    print(json.dumps({k:v for k,v in result.items() if k!='deltas'},indent=2))
    print('seconds',time.time()-t)
