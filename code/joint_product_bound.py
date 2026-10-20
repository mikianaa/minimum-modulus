"""Certified joint product/forbidden-divisor correction.

The bounded complement vanishes beyond an explicit product cutoff, so
all exponent tails are accounted for by the exact mean of B.
"""
from verify_minimum_modulus import I,S,ceildiv
from verify_stoploss import L


def product_distribution_upper(primes,parameters,coefficients,H):
    dist={(1,0):S}
    for q,c in zip(primes,coefficients):
        D=parameters.get(q,0);den=q*(L-D);assert den>=L
        new={}
        for (z,j),mass in dist.items():
            key=(z,j);new[key]=new.get(key,0)+ceildiv(mass*(den-L),den)
            qp=q*q
            for a in range(1,H//z):
                key=(z*(a+1),j+c)
                new[key]=new.get(key,0)+ceildiv(mass*L*(q-1),(L-D)*qp)
                qp*=q
        dist=new
    return dist


def joint_product_correction_lower(small,K,p,D,primes,parameters,q0):
    assert p>=q0 and (K-1)//p<q0*q0 and p*p*q0>=K
    powers=[];pw=p
    while pw<K:powers.append(pw);pw*=p
    base=tuple((K-1)//pw+1 for pw in powers)
    qs=tuple((K-1)//(p*q)+1 for q in primes)
    thresholds=tuple(dict.fromkeys(base+qs));positions={s:j+1 for j,s in enumerate(thresholds)}
    bden=powers[-1];cache={};total=0
    H=ceildiv((p-1)*D,L)+1+sum(p*q<K for q in primes)
    for key,weight in small.groups(thresholds):
        coefficients=tuple(key[positions[s]] for s in qs)
        if coefficients not in cache:
            dist=product_distribution_upper(primes,parameters,coefficients,H)
            extra=I.rat(0)
            for q,c in zip(primes,coefficients):
                extra=extra+I.rat((p-1)*c*L,p*q*(L-parameters.get(q,0)))
            cache[coefficients]=dist,extra
        dist,extra=cache[coefficients]
        bnum=sum((p-1)*(bden//pw)*key[positions[s]] for pw,s in zip(powers,base))
        complement=0
        for (z,j),mass in dist.items():
            btotal=bnum+(p-1)*(bden//p)*j
            value=min(btotal*L,max(0,(D*(p-1)-key[0]*z*L)*bden+btotal*L))
            complement+=ceildiv(mass*value,bden*L)
        gain=max(0,(I.rat(bnum,bden)+extra).lo-complement)
        total+=weight*gain//S
    return total
