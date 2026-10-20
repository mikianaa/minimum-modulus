"""Integer-only verification with any finite initial set of conditioned primes.

The witness supplies only integers. Grouping states changes no inequality.
See the general conditioning proposition in the accompanying LaTeX report.
"""
import argparse,bisect,hashlib,itertools,json,math,time
from pathlib import Path
from verify_stoploss import ProductBound,L
from moment_product_bound import MomentProductBound
from verify_minimum_modulus import (S,I,ONE,FOUR,ceildiv,primes_to,
    factorizations,initial,moments,stopping_bound,fixed_text)


class StrongSmallBound:
    def __init__(self,primes,exponents,max_threshold):
        assert len(primes)==len(exponents) and primes
        assert all(isinstance(e,int) and 0<=e<=100 for e in exponents)
        self.mean_num=math.prod(primes);self.mean_den=math.prod(q-1 for q in primes)
        numerator=self.mean_den
        self.states=[];mean_lower=0;divisor_cache={}
        clips=[]
        for q in primes:
            e=0;power=q
            while power<max_threshold:e+=1;power*=q
            clips.append(e)
        for es in itertools.product(*(range(e+1) for e in exponents)):
            tau=math.prod(e+1 for e in es)
            den=math.prod(q**(e+1) for q,e in zip(primes,es))
            low=numerator*tau*S//den;up=ceildiv(numerator*tau*S,den)
            mean_lower+=low
            clipped=tuple(min(e,c) for e,c in zip(es,clips))
            if clipped not in divisor_cache:
                divisors=[1]
                for q,e in zip(primes,clipped):
                    divisors=[d*q**a for d in divisors for a in range(e+1)
                              if d*q**a<max_threshold]
                divisor_cache[clipped]=sorted(divisors)
            self.states.append((tau,up,divisor_cache[clipped]))
        self.tail_mean_upper=ceildiv(self.mean_num*S,self.mean_den)-mean_lower
        assert self.tail_mean_upper>=0
        self.cache={}

    def groups(self,thresholds):
        if thresholds not in self.cache:
            groups={}
            for tau,weight,divisors in self.states:
                key=(tau,)+tuple(bisect.bisect_left(divisors,s) for s in thresholds)
                groups[key]=groups.get(key,0)+weight
            self.cache[thresholds]=list(groups.items())
        return self.cache[thresholds]

    def loss_upper(self,K,p,D,mean_upper,product):
        assert 0<=D<=L//2
        prefix=product.prefixes();powers=[];power=p
        while power<K:powers.append(power);power*=p
        thresholds=tuple((K-1)//pp+1 for pp in powers)
        bden=powers[-1] if powers else 1
        multipliers=[(p-1)*(bden//pp) for pp in powers]
        total=ceildiv(mean_upper*self.tail_mean_upper,S)
        for key,weight in self.groups(thresholds):
            tau=key[0];bnum=sum(m*c for m,c in zip(multipliers,key[1:]))
            tlo=((D*(p-1)*bden+L*bnum)*S)//(L*bden*tau)
            e=product.excess_upper(tlo,mean_upper,prefix)
            total+=ceildiv(weight*e,S)
        return ceildiv(total*L,(p-1)*(L-D))


def verify(path,output):
    if not __debug__:raise RuntimeError('Assertions must be enabled.')
    raw=Path(path).read_bytes();w=json.loads(raw)
    K=w['K'];start=w['start'];ds=w['deltas'];cap=w['cap'];resolution=w['resolution']
    es=w['small_exponents'];r=len(es)
    assert isinstance(K,int) and 2<K<40000
    assert isinstance(start,int) and start>=r>=2
    assert ds and all(isinstance(d,int) and 0<=d<=L//2 for d in ds)
    n=start+len(ds);limit=max(K,20*n);ps=primes_to(limit)
    while len(ps)<n:limit*=2;ps=primes_to(limit)
    assert ps[:r]==w['small_primes']
    fs=factorizations((K-1)//ps[start]+1)
    mu,C,T=initial(K,ps,start);initial_mu=mu
    product_class=MomentProductBound if w.get('moment_bins',False) else ProductBound
    prod=product_class(cap,resolution)
    small=StrongSmallBound(ps[:r],es,(K-1)//ps[start]+1)
    for q in ps[r:start]:prod.update(q,0)
    nu={};rows=[];wins={'first':0,'second':0,'conditional':0,'unconditional':0}
    switched=False
    for idx,D in enumerate(ds,start):
        p=ps[idx];i=idx+1;d=I.rat(D,L)
        if w.get('switch',False) and p*2>=K and not switched:
            prod=product_class(cap,resolution)
            for j,q in enumerate(ps[:idx]):prod.update(q,0 if j<start else ds[j-start])
            switched=True
        if p<K:m1,m2=moments(K,p,C,T,nu,fs)
        else:m1,m2=C*I.rat(1,p-1),T*I.rat(1,(p-1)**2)
        candidates=[(m1.hi,'first')]
        if D:candidates.append(((m2/(FOUR*d*(ONE-d))).hi,'second'))
        if switched:
            # At 2p>=K the only forbidden old divisor is 1 (if p<K).
            shift_num=p-1 if p<K else 0
            tlo=((D*(p-1)*p+L*shift_num)*S)//(L*p)
            e=prod.excess_upper(tlo,C.hi,prod.prefixes())
            candidates.append((ceildiv(e*L,(p-1)*(L-D)),'unconditional'))
        else:
            crest=ceildiv(C.hi*small.mean_den,small.mean_num)
            candidates.append((small.loss_upper(K,p,D,crest,prod),'conditional'))
        loss,kind=min(candidates);wins[kind]+=1;mu-=loss
        assert mu>0,('nonpositive survivor bound',i,p)
        v=ONE/(ONE-d)
        C=C*(ONE+v*I.rat(1,p-1));T=T*(ONE+v*I.rat(3*p-1,(p-1)**2));nu[p]=v
        prod.update(p,D)
        if i in (100,1000,10000,25000,50000) or i==n:
            row=dict(i=i,p=p,mu_lower=fixed_text(mu),T_upper=fixed_text(T.hi),
                f_upper=fixed_text(ceildiv(T.hi*S,mu)),wins=wins.copy())
            rows.append(row);print('Verified stage',i,'prime',p,flush=True)
    fu=ceildiv(T.hi*S,mu);stop=stopping_bound(n)
    assert n>=10 and fu<stop,('analytic stopping criterion failed',n,fixed_text(fu),fixed_text(stop))
    result=dict(status='VERIFIED',cutoff=K,minimum_modulus_at_most=K-1,
        method='general small-prime conditioning with downward product bins',decimal_places=40,
        initial_mu_lower=fixed_text(initial_mu),terminal_index=n,terminal_prime=ps[n-1],
        f_upper=fixed_text(fu),stopping_lower=fixed_text(stop),certified_gap=fixed_text(stop-fu),
        mu_lower=fixed_text(mu),bins=len(prod.bins),cap=cap,resolution=resolution,
        omitted_probability_upper=fixed_text(prod.error),small_tail_mean_upper=fixed_text(small.tail_mean_upper),
        small_primes=ps[:r],small_exponents=es,
        witness_sha256=hashlib.sha256(raw).hexdigest(),wins=wins,checkpoints=rows)
    Path(output).write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness');ap.add_argument('--output',required=True)
    args=ap.parse_args();t=time.time();print(json.dumps(verify(args.witness,args.output),indent=2))
    print('seconds',time.time()-t)
