"""Integer-arithmetic certificate for the conditional convex-loss bound.

Standard library only. The new comparison lemma and rounding argument are
proved in src/conditional-stoploss-improvement.tex. Experimental metadata and
floating-point estimates are never trusted. Zero distortion is allowed.
"""
import argparse,bisect,gzip,hashlib,json,time
from pathlib import Path
from verify_minimum_modulus import (S,I,ONE,FOUR,ceildiv,primes_to,
    factorizations,initial,moments,stopping_bound,fixed_text)

L=10**12


def make_bins(cap,resolution):
    assert cap>=2 and resolution>=1
    result=[0,1]
    while result[-1]<cap:
        result.append(min(cap,result[-1]+max(1,result[-1]//resolution)))
    return result


class ProductBound:
    def __init__(self,cap,resolution):
        self.bins=make_bins(cap,resolution)
        self.mass=[0]*len(self.bins);self.mass[1]=S
        self.error=0
        self.maps={}
        self.tail_cutoff=2*10**32

    def update(self,q,D):
        # The exact exponent law has P(E>=a)=nu/q**a, nu=L/(L-D).
        assert q>=2 and 0<=D<=L//2 and L<=q*(L-D)
        den0=q*(L-D);num0=den0-L
        new=[ceildiv(u*num0,den0) if u else 0 for u in self.mass]
        H=1;qH=q
        while qH<self.tail_cutoff:
            H+=1;qH*=q
        self.error+=ceildiv(L*S,(L-D)*qH)
        power=q*q;num=L*(q-1)
        for h in range(2,H+1):
            if h not in self.maps:
                mapping=[]
                for j,b in enumerate(self.bins):
                    if j==0:continue
                    if b*h>self.bins[-1]:break
                    mapping.append((j,bisect.bisect_right(self.bins,b*h)-1))
                self.maps[h]=mapping
            den=(L-D)*power
            for j,target in self.maps[h]:
                if self.mass[j]:new[target]+=ceildiv(self.mass[j]*num,den)
            power*=q
        self.mass=new

    def prefixes(self):
        probs=[];weighted=[];v=w=0
        for b,u in zip(self.bins,self.mass):
            v+=u;w+=b*u;probs.append(v);weighted.append(w)
        return probs,weighted

    def excess_upper(self,tlo,mean_upper,prefixes):
        # True E(Z-t)_+ decreases in t. Binning rounds Z DOWN only in
        # E(t-Z)_+; E Z is independently bounded by its exact Euler product.
        tlo=min(tlo,self.bins[-1]*S)
        assert tlo>=0
        j=bisect.bisect_right(self.bins,tlo//S)-1
        probs,weighted=prefixes
        comp=ceildiv(tlo*probs[j]-S*weighted[j],S)
        omitted=ceildiv(tlo*self.error,S)
        upper=mean_upper-tlo+comp+omitted
        assert upper>=0
        return upper


class SmallBound:
    def __init__(self,amax,bmax):
        assert 0<=amax<=100 and 0<=bmax<=100
        self.states=[];mean_lower=0
        for a in range(amax+1):
            for b in range(bmax+1):
                tau=(a+1)*(b+1);den=2**(a+1)*3**(b+1)
                self.states.append((a,b,tau,den))
                mean_lower+=2*tau*S//den
        self.tail_mean_upper=3*S-mean_lower
        assert self.tail_mean_upper>=0
        self.count_cache={}

    def counts(self,s):
        if s not in self.count_cache:
            result=[]
            for a,b,_,_ in self.states:
                count=0;x=1
                for aa in range(a+1):
                    if x>=s:break
                    y=x
                    for bb in range(b+1):
                        if y>=s:break
                        count+=1;y*=3
                    x*=2
                result.append(count)
            self.count_cache[s]=result
        return self.count_cache[s]

    def loss_upper(self,K,p,D,mean_upper,product):
        assert 0<D<=L//2
        prefix=product.prefixes()
        thresholds=[];power=p
        while power<K:
            thresholds.append((power,self.counts((K-1)//power+1)))
            power*=p
        bden=thresholds[-1][0] if thresholds else 1
        bnums=[0]*len(self.states)
        for pp,counts in thresholds:
            mult=(p-1)*(bden//pp)
            for j,c in enumerate(counts):bnums[j]+=mult*c
        total=ceildiv(mean_upper*self.tail_mean_upper,S)
        for (_,_,tau,den),bnum in zip(self.states,bnums):
            # Exact threshold is ((p-1)*delta + B_small)/tau.
            tlo=((D*(p-1)*bden+L*bnum)*S)//(L*bden*tau)
            e=product.excess_upper(tlo,mean_upper,prefix)
            total+=ceildiv(2*tau*e,den)
        return ceildiv(total*L,(p-1)*(L-D))


def verify(path,output):
    if not __debug__:raise RuntimeError('Do not run certificate checks with Python -O.')
    raw=Path(path).read_bytes()
    if str(path).endswith('.gz'):raw=gzip.decompress(raw)
    w=json.loads(raw)
    K=w['K'];start=w['start'];ds=w['deltas'];cap=w['cap'];resolution=w['resolution']
    assert isinstance(K,int) and 2<K<615300
    assert isinstance(start,int) and start>=2
    assert ds and all(isinstance(d,int) and 0<=d<=L//2 for d in ds)
    amax,bmax=w['small_exponents']
    n=start+len(ds);limit=max(K,20*n);ps=primes_to(limit)
    while len(ps)<n:limit*=2;ps=primes_to(limit)
    fs=factorizations((K-1)//ps[start]+1)
    mu,C,T=initial(K,ps,start);initial_mu=mu
    prod=ProductBound(cap,resolution);small=SmallBound(amax,bmax)
    for q in ps[2:start]:prod.update(q,0)
    nu={};rows=[];wins={'first':0,'second':0,'conditional':0}
    for idx,D in enumerate(ds,start):
        p=ps[idx];i=idx+1;d=I.rat(D,L)
        if p<K:m1,m2=moments(K,p,C,T,nu,fs)
        else:m1,m2=C*I.rat(1,p-1),T*I.rat(1,(p-1)**2)
        candidates=[(m1.hi,'first')]
        if D:
            candidates.append(((m2/(FOUR*d*(ONE-d))).hi,'second'))
            crest_upper=ceildiv(C.hi,3)
            candidates.append((small.loss_upper(K,p,D,crest_upper,prod),'conditional'))
        loss,kind=min(candidates);wins[kind]+=1;mu-=loss
        assert mu>0,('nonpositive survivor bound',i,p)
        v=ONE/(ONE-d)
        C=C*(ONE+v*I.rat(1,p-1));T=T*(ONE+v*I.rat(3*p-1,(p-1)**2));nu[p]=v
        prod.update(p,D)
        if i in (100,1000,10000,25000) or i==n:
            row=dict(i=i,p=p,mu_lower=fixed_text(mu),T_upper=fixed_text(T.hi),
                     f_upper=fixed_text(ceildiv(T.hi*S,mu)),wins=wins.copy())
            rows.append(row);print('Verified stage',i,'prime',p,flush=True)
    fu=ceildiv(T.hi*S,mu);stop=stopping_bound(n)
    assert n>=10 and fu<stop,('analytic stopping criterion failed',n)
    result=dict(status='VERIFIED',cutoff=K,minimum_modulus_at_most=K-1,
        method='conditional convex loss with downward product bins',decimal_places=40,
        initial_mu_lower=fixed_text(initial_mu),terminal_index=n,terminal_prime=ps[n-1],
        f_upper=fixed_text(fu),stopping_lower=fixed_text(stop),certified_gap=fixed_text(stop-fu),
        mu_lower=fixed_text(mu),bins=len(prod.bins),cap=cap,resolution=resolution,
        omitted_probability_upper=fixed_text(prod.error),small_tail_mean_upper=fixed_text(small.tail_mean_upper),
        witness_sha256=hashlib.sha256(raw).hexdigest(),wins=wins,checkpoints=rows)
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness',nargs='?',default='data/stoploss-certificate-40000.json')
    ap.add_argument('--output',default='data/stoploss-verified-40000.json');args=ap.parse_args()
    t=time.time();r=verify(args.witness,args.output)
    print(json.dumps(r,indent=2));print('seconds',time.time()-t)
