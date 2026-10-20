"""Floating counterpart of joint_product_bound, for search only."""
import math
import numpy as np
from functools import lru_cache


def product_distribution(primes,nus,coefficients,H):
    return _product_distribution(tuple(primes),tuple(nus),tuple(coefficients),H)


@lru_cache(maxsize=1024)
def _product_distribution(primes,nus,coefficients,H):
    # Dense vector operations implement the same finite recurrence. No
    # probability cutoff is used: all products z <= H are retained.
    width=1+sum(coefficients)
    dist=np.zeros((H+1,width));dist[1,0]=1.;used=0
    for q,nu,c in zip(primes,nus,coefficients):
        new=(1-nu/q)*dist
        source,multiplier,target=_transitions(H,width,c,used)
        prob=nu*(q-1)*np.power(float(q),-multiplier)
        weights=(dist[source,:used+1]*prob[:,None]).ravel()
        new+=np.bincount(target,weights=weights,minlength=new.size).reshape(new.shape)
        used+=c;dist=new
    z,j=np.nonzero(dist)
    return {(int(a),int(b)):float(dist[a,b]) for a,b in zip(z,j)}


@lru_cache(maxsize=512)
def _transitions(H,width,c,used):
    source=np.asarray([z for m in range(2,H+1) for z in range(1,H//m+1)])
    multiplier=np.asarray([m for m in range(2,H+1) for z in range(1,H//m+1)])
    source=source.astype(int);multiplier=multiplier.astype(int)
    target=((source*multiplier)[:,None]*width+np.arange(c,c+used+1)).ravel()
    return source,multiplier,target


def reference_product_distribution(primes,nus,coefficients,H):
    """Original sparse recurrence, retained as an independent search check."""
    dist={(1,0):1.}
    for q,nu,c in zip(primes,nus,coefficients):
        new={}
        for (z,j),mass in dist.items():
            key=(z,j);new[key]=new.get(key,0.)+mass*(1-nu/q)
            prob=nu*(q-1)/q**2
            for a in range(1,H//z):
                key=(z*(a+1),j+c);new[key]=new.get(key,0.)+mass*prob
                prob/=q
        dist=new
    return dist


def complementary_values(dist,A,tau,B,step,subtract_unshifted=False):
    """Evaluate positive parts by cumulative mass and first moment in z."""
    columns={}
    for (z,j),mass in dist.items():columns.setdefault(j,[]).append((z,mass))
    result=np.zeros(len(tau))
    for j,entries in columns.items():
        entries.sort();z=np.asarray([v[0] for v in entries]);mass=np.asarray([v[1] for v in entries])
        P=np.r_[0.,np.cumsum(mass)];M=np.r_[0.,np.cumsum(mass*z)]
        threshold=A+B+step*j;k=np.searchsorted(z,threshold/tau,side='right')
        value=threshold*P[k]-tau*M[k]
        if subtract_unshifted:
            k=np.searchsorted(z,A/tau,side='right')
            value-=A*P[k]-tau*M[k]
        result+=np.maximum(0.,value)
    return result


def correction(joint,p,d,primes,nus):
    assert (joint.K-1)//p<joint.q0**2 and p*p*joint.q0>=joint.K
    masks,tau,B,weights=joint.groups(p);A=(p-1)*d;total=0.
    H=math.ceil(A)+1+sum(p*q<joint.K for q in primes)
    for mask in np.unique(masks):
        ix=masks==mask;coeff=[]
        for q in primes:
            mx=(joint.K-1)//(p*q)
            coeff.append(int(mx>=1)+sum((int(mask)>>(j-2))&1 for j in range(2,mx+1)))
        dist=product_distribution(primes,nus,coeff,H)
        extra=(p-1)/p*sum(c*v/q for q,v,c in zip(primes,nus,coeff))
        ts=tau[ix];bs=B[ix];ws=weights[ix]
        complement=complementary_values(dist,A,ts,bs,(p-1)/p,True)
        total+=float(ws@np.maximum(0.,bs+extra-complement))
    return total


def conditioned_loss(joint,p,d,primes,nus,small_primes,small_nus,state_cap,m1):
    """Exact retained complement and a divisor-pair bound outside it.

    Search only. The independent integer counterpart certifies parameters.
    """
    assert state_cap>=((joint.K-1)//p+1)**2
    assert (joint.K-1)//p<joint.q0**2 and p*p*joint.q0>=joint.K
    masks,tau,B,weights=joint.groups(p);A=(p-1)*d;total=0.
    H=math.ceil(A)+1+sum(p*q<joint.K for q in primes)
    for mask in np.unique(masks):
        ix=masks==mask;coeff=[]
        for q in primes:
            mx=(joint.K-1)//(p*q)
            coeff.append(int(mx>=1)+sum((int(mask)>>(j-2))&1 for j in range(2,mx+1)))
        dist=product_distribution(primes,nus,coeff,H)
        ts=tau[ix];bs=B[ix];ws=weights[ix]
        total+=float(ws@complementary_values(dist,A,ts,bs,(p-1)/p))
    T=math.floor(2*A)
    if T:
        small_dist=product_distribution(small_primes,small_nus,[0]*len(small_primes),T)
        retained=np.bincount(tau.astype(int),weights=weights,minlength=T+1)
        rest=product_distribution(primes,nus,[0]*len(primes),T)
        zz=np.asarray([z for z,j in rest]);ww=np.asarray(list(rest.values()))
        for (t,_),mass in small_dist.items():
            missing=max(0.,mass-retained[t])
            total+=missing*float(np.maximum(0.,A-t*zz/2)@ww)
    return (m1-d+total/(p-1))/(1-d)
