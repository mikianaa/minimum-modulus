"""Exact lower bounds for joint large-prime presence corrections."""
from verify_minimum_modulus import S
from verify_stoploss import L


def presence_distribution_lower(primes,parameters,coefficients,count_cap=12):
    dist={(0,0):S}
    for q,c in zip(primes,coefficients):
        D=parameters.get(q,0);den=q*(L-D);assert 0<=L<=den
        new={}
        for (a,j),mass in dist.items():
            key=(a,j);new[key]=new.get(key,0)+mass*(den-L)//den
            key=(min(count_cap,a+1),j+c);new[key]=new.get(key,0)+mass*L//den
        dist=new
    return dist


def joint_correction_lower(small,K,p,D,primes,parameters,q0,rest=None):
    assert p>=q0 and (K-1)//p<q0*q0 and p*p*q0>=K
    powers=[];power=p
    while power<K:powers.append(power);power*=p
    base=tuple((K-1)//pp+1 for pp in powers)
    qs=tuple((K-1)//(p*q)+1 for q in primes)
    thresholds=tuple(dict.fromkeys(base+qs));positions={s:j+1 for j,s in enumerate(thresholds)}
    bden=powers[-1];cache={};total=0
    prefix=rest.prefixes() if rest is not None else None
    for key,weight in small.groups(thresholds):
        coefficients=tuple(key[positions[s]] for s in qs)
        if coefficients not in cache:
            cache[coefficients]=presence_distribution_lower(primes,parameters,coefficients)
        bnum=sum((p-1)*(bden//pp)*key[positions[s]] for pp,s in zip(powers,base))
        gain=0
        for (a,j),mass in cache[coefficients].items():
            btotal=bnum+(p-1)*(bden//p)*j
            excess=max(0,(key[0]*(2**a)*L-D*(p-1))*bden)
            gain+=mass*min(btotal*L,excess)//(bden*L)
        if rest is not None:
            gain=max(gain,rest.gain_lower(key[0],D*(p-1),L,bnum,bden,prefix))
        total+=weight*gain//S
    return total
