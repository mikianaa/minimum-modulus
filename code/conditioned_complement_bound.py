"""Outward integer enclosure of the restricted-loss complement.

Retained small-prime configurations are integrated exactly by the joint
product DP. Outside them, divisor pairing bounds W >= tau/2. Only a finite
divisor-count distribution is needed to bound that remaining complement.
"""
import bisect
from joint_product_bound import product_distribution_upper
from verify_minimum_modulus import I,ONE,S,ceildiv
from verify_stoploss import L


def upper_groups(small,thresholds):
    assert small.processed==len(small.primes)
    if not hasattr(small,'complement_upper_cache'):small.complement_upper_cache={}
    if thresholds not in small.complement_upper_cache:
        groups={}
        for (tau,_,divisors),coefficient in zip(small.original_states,small.coeff_hi):
            key=(tau,)+tuple(bisect.bisect_left(divisors,s) for s in thresholds)
            weight=ceildiv(small.norm_hi*coefficient,S)
            groups[key]=groups.get(key,0)+weight
        small.complement_upper_cache[thresholds]=list(groups.items())
    return small.complement_upper_cache[thresholds]


def loss_upper(small,K,p,D,primes,parameters,q0,state_cap,M1):
    assert state_cap>=((K-1)//p+1)**2
    assert p>=q0 and (K-1)//p<q0*q0 and p*p*q0>=K
    powers=[];pw=p
    while pw<K:powers.append(pw);pw*=p
    base=tuple((K-1)//pw+1 for pw in powers)
    qs=tuple((K-1)//(p*q)+1 for q in primes)
    thresholds=tuple(dict.fromkeys(base+qs));positions={s:j+1 for j,s in enumerate(thresholds)}
    bden=powers[-1];cache={};total=0
    H=ceildiv((p-1)*D,L)+1+sum(p*q<K for q in primes)
    for key,weight in upper_groups(small,thresholds):
        coefficients=tuple(key[positions[s]] for s in qs)
        if coefficients not in cache:
            cache[coefficients]=product_distribution_upper(primes,parameters,coefficients,H)
        bnum=sum((p-1)*(bden//pw)*key[positions[s]] for pw,s in zip(powers,base))
        complement=0
        for (z,j),mass in cache[coefficients].items():
            btotal=bnum+(p-1)*(bden//p)*j
            value=max(0,(D*(p-1)-key[0]*z*L)*bden+btotal*L)
            complement+=ceildiv(mass*value,bden*L)
        total+=ceildiv(weight*complement,S)
    T=2*(p-1)*D//L
    outside=0
    if T:
        full=product_distribution_upper(small.primes,parameters,[0]*len(small.primes),T)
        retained={}
        for tau,weight,_ in small.states:retained[tau]=retained.get(tau,0)+weight
        rest=product_distribution_upper(primes,parameters,[0]*len(primes),T)
        for (tau,_),mass in full.items():
            missing=max(0,mass-retained.get(tau,0));complement=0
            for (z,_),weight in rest.items():
                value=max(0,2*(p-1)*D-tau*z*L)
                complement+=ceildiv(weight*value,2*L)
            outside+=ceildiv(missing*complement,S)
    d=I.rat(D,L)
    value=(M1-d+I.rat(total+outside,S*(p-1)))/(ONE-d)
    assert value.hi>=0
    return value.hi,outside
