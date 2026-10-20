"""Integer certificate for finite conditioning with an unrestricted fallback."""
import argparse,bisect,hashlib,json,math,time
from pathlib import Path
from verify_stoploss import ProductBound,L
from moment_product_bound import MomentProductBound
from verify_minimum_modulus import (S,I,ONE,FOUR,ceildiv,primes_to,
    factorizations,initial,moments,stopping_bound,fixed_text)


class LowerProductBound(ProductBound):
    def update(self,q,D):
        assert q>=2 and 0<=D<=L//2
        den0=q*(L-D);num0=den0-L
        new=[u*num0//den0 for u in self.mass]
        H=1;qH=q
        while qH<self.tail_cutoff:H+=1;qH*=q
        power=q*q;num=L*(q-1)
        for h in range(2,H+1):
            if h not in self.maps:
                self.maps[h]=[(j,bisect.bisect_right(self.bins,b*h)-1)
                              for j,b in enumerate(self.bins) if j and b*h<=self.bins[-1]]
            den=(L-D)*power
            for j,target in self.maps[h]:
                if self.mass[j]:new[target]+=self.mass[j]*num//den
            power*=q
        self.mass=new

    def gain_lower(self,tau,anum,aden,bnum,bden,prefix):
        # g(z)=min(B,(tau*z-A)_+) is increasing. Retained lower bin masses
        # and downward representatives give a lower expectation directly.
        P,W=prefix
        j=bisect.bisect_right(self.bins,anum//(aden*tau))-1
        k=bisect.bisect_right(self.bins,(anum*bden+bnum*aden)//(aden*bden*tau))-1
        gain=(tau*aden*bden*(W[k]-W[j])-anum*bden*(P[k]-P[j])
              +bnum*aden*(P[-1]-P[k]))//(aden*bden)
        assert gain>=0
        return gain


class SparseSmallBound:
    def __init__(self,primes,state_cap,max_threshold):
        assert 2<=state_cap and max_threshold>=1
        self.states=[];self.ns=[];self.masks=[];numerator=math.prod(q-1 for q in primes);denominator=math.prod(primes)
        self.probability_lower=0;self.probability_upper=0
        def visit(j,n,tau,divisors,mask):
            if j==len(primes):
                weight=numerator*S//(denominator*n)
                self.states.append((tau,weight,sorted(divisors)))
                self.ns.append(n);self.masks.append(mask)
                self.probability_lower+=weight
                self.probability_upper+=ceildiv(numerator*S,denominator*n)
                return
            q=primes[j];power=1;a=0;newdivisors=[]
            while n*power<state_cap:
                newdivisors.extend(d*power for d in divisors if d*power<max_threshold)
                visit(j+1,n*power,tau*(a+1),newdivisors.copy(),mask|((1<<j) if a else 0))
                power*=q;a+=1
        visit(0,1,1,[1],0);self.cache={}
        assert self.probability_lower<=S

    def groups(self,thresholds):
        if thresholds not in self.cache:
            groups={}
            for tau,weight,divisors in self.states:
                key=(tau,)+tuple(bisect.bisect_left(divisors,s) for s in thresholds)
                groups[key]=groups.get(key,0)+weight
            self.cache[thresholds]=list(groups.items())
        return self.cache[thresholds]

    def correction_lower(self,K,p,D,rest):
        if p>=K:return 0
        prefix=rest.prefixes();powers=[];power=p
        while power<K:powers.append(power);power*=p
        thresholds=tuple((K-1)//pp+1 for pp in powers)
        bden=powers[-1];multipliers=[(p-1)*(bden//pp) for pp in powers]
        total=0
        for key,weight in self.groups(thresholds):
            tau=key[0];bnum=sum(m*c for m,c in zip(multipliers,key[1:]))
            gain=rest.gain_lower(tau,D*(p-1),L,bnum,bden,prefix)
            total+=weight*gain//S
        return total


class DynamicSmallBound(SparseSmallBound):
    def __init__(self,primes,state_cap,max_threshold):
        super().__init__(primes,state_cap,max_threshold)
        self.primes=primes;self.original_states=self.states
        self.coeff_lo=[S//n for n in self.ns];self.coeff_hi=[ceildiv(S,n) for n in self.ns]
        self.norm_lo=self.norm_hi=S;self.processed=0;self.refresh()

    def refresh(self):
        active=[];lo=hi=0;mask_limit=1<<self.processed
        for state,mask,c,d in zip(self.original_states,self.masks,self.coeff_lo,self.coeff_hi):
            if mask>=mask_limit:continue
            w=self.norm_lo*c//S;u=ceildiv(self.norm_hi*d,S)
            active.append((state[0],w,state[2]));lo+=w;hi+=u
        self.states=active;self.probability_lower=lo;self.probability_upper=min(S,hi);self.cache={}

    def update(self,q,D):
        if q>self.primes[-1]:return
        assert self.primes[self.processed]==q
        den=q*(L-D)-L;assert den>0
        if D:
            num=L*(q-1);bit=1<<self.processed
            for j,mask in enumerate(self.masks):
                if mask&bit:
                    self.coeff_lo[j]=self.coeff_lo[j]*num//den
                    self.coeff_hi[j]=ceildiv(self.coeff_hi[j]*num,den)
        self.norm_lo=self.norm_lo*den//(q*(L-D))
        self.norm_hi=ceildiv(self.norm_hi*den,q*(L-D));self.processed+=1;self.refresh()


def verify(path,output,finite_only=False):
    if not __debug__:raise RuntimeError('Assertions must be enabled.')
    raw=Path(path).read_bytes();w=json.loads(raw)
    K=w['K'];start=w['start'];r=w['small'];ds=w['deltas'];cap=w['cap'];resolution=w['resolution']
    state_cap=w['state_cap']
    dynamic=w.get('dynamic',False);joint_presence=w.get('joint_presence',False)
    statewise_presence=w.get('statewise_presence',False)
    joint_product=w.get('joint_product',False);assert isinstance(joint_product,bool)
    assert isinstance(statewise_presence,bool)
    assert isinstance(dynamic,bool) and isinstance(joint_presence,bool)
    assert isinstance(K,int) and 2<K<30000
    assert isinstance(start,int) and isinstance(r,int) and start>=0 and r>=2
    assert dynamic or start>=r
    assert ds and all(isinstance(d,int) and 0<=d<=L//2 for d in ds)
    n=start+len(ds);limit=max(K,20*n);ps=primes_to(limit)
    while len(ps)<n:limit*=2;ps=primes_to(limit)
    fs=factorizations((K-1)//ps[start]+1)
    mu,C,T=initial(K,ps,start);initial_mu=mu
    support_mode=w.get('support_mode','uniform')
    normalize_conditioning=w.get('normalize_conditioning',False)
    assert isinstance(normalize_conditioning,bool)
    prefix_resolution=w.get('prefix_resolution',resolution)
    prefix_index=w.get('prefix_index',0)
    assert isinstance(prefix_resolution,int) and prefix_resolution>=resolution
    assert isinstance(prefix_index,int) and (prefix_index==0 or start<prefix_index<n)
    coarse=None;coarse_rest=None
    if prefix_index:
        coarse=MomentProductBound(cap,resolution,support_mode)
        coarse_rest=LowerProductBound(cap,resolution)
    full=MomentProductBound(cap,prefix_resolution if prefix_index else resolution,support_mode)
    rest=LowerProductBound(cap,prefix_resolution if prefix_index else resolution)
    product_tail_power=w.get('product_tail_power',32)
    assert isinstance(product_tail_power,int) and 12<=product_tail_power<=32
    for distribution in (full,rest,coarse,coarse_rest):
        if distribution is not None:distribution.tail_cutoff=2*10**product_tail_power
    small=(DynamicSmallBound if dynamic else SparseSmallBound)(ps[:r],state_cap,(K-1)//ps[start]+1)
    for j,q in enumerate(ps[:start]):
        full.update(q,0)
        if coarse is not None:coarse.update(q,0)
        if j>=r:
            rest.update(q,0)
            if coarse_rest is not None:coarse_rest.update(q,0)
        if dynamic:small.update(q,0)
    print('Conditioning states',len(small.states),flush=True)
    nu={};parameters={};rows=[];wins={'first':0,'second':0,'sparse':0,'early_exact':0,'conditioned_complement':0};joint_wins=0
    conditioned_complement_prime=w.get('conditioned_complement_prime',0)
    assert isinstance(conditioned_complement_prime,int) and 0<=conditioned_complement_prime<K
    assert not conditioned_complement_prime or dynamic
    exact_early_prime=w.get('exact_early_prime',0);early_states={}
    atom_cap=w.get('recursive_atom_cap',0);atom_prime=w.get('recursive_atom_prime',0)
    assert isinstance(atom_cap,int) and atom_cap>=0 and isinstance(atom_prime,int)
    atoms=None
    if atom_cap:
        from recursive_atom_bound import RecursiveAtoms
        assert atom_cap>=2 and 0<atom_prime<K
        atoms=RecursiveAtoms(K,atom_cap,w.get('recursive_density_degree',1))
        for p in ps[:start]:atoms.update(p,0)
        wins['recursive_first']=0;wins['recursive_second']=0
    for idx,D in enumerate(ds,start):
        p=ps[idx];i=idx+1;d=I.rat(D,L)
        if coarse is not None and idx==prefix_index:
            # Both distributions have been propagated from the beginning.
            # No bins are reinterpreted or converted at the transition.
            full=coarse;coarse=None;rest=coarse_rest;coarse_rest=None
        if p<K:m1,m2=moments(K,p,C,T,nu,fs)
        else:m1,m2=C*I.rat(1,p-1),T*I.rat(1,(p-1)**2)
        candidates=[(m1.hi,'first')]
        if D:candidates.append(((m2/(FOUR*d*(ONE-d))).hi,'second'))
        baseline=full.excess_upper(D*(p-1)*S//L,C.hi,full.prefixes())
        correction=small.correction_lower(K,p,D,rest)
        joint_valid=idx>=r and p*ps[r]<K and (K-1)//p<ps[r]**2 and p*p*ps[r]>=K
        if joint_presence and joint_valid:
            from joint_presence_bound import joint_correction_lower
            joint=joint_correction_lower(small,K,p,D,ps[r:idx],parameters,ps[r],rest if statewise_presence else None)
            if joint>correction:correction=joint;joint_wins+=1
        if joint_product and joint_valid:
            from joint_product_bound import joint_product_correction_lower
            correction=max(correction,joint_product_correction_lower(small,K,p,D,ps[r:idx],parameters,ps[r]))
        if normalize_conditioning:
            correction=correction*S//min(S,small.probability_upper)
        assert baseline>=correction
        candidates.append((ceildiv((baseline-correction)*L,(p-1)*(L-D)),'sparse'))
        if p<=exact_early_prime:
            from early_exact_states import loss_upper
            value,count=loss_upper(K,p,D,ps[:idx],parameters,m1)
            candidates.append((value,'early_exact'));early_states[p]=count
        if p<=conditioned_complement_prime and idx>=r and (K-1)//p<ps[r]**2 and p*p*ps[r]>=K:
            from conditioned_complement_bound import loss_upper as conditioned_loss_upper
            value,_=conditioned_loss_upper(small,K,p,D,ps[r:idx],parameters,ps[r],state_cap,m1)
            candidates.append((value,'conditioned_complement'))
        if atoms is not None and p<=atom_prime:
            gains=atoms.moment_gains(p)
            am1=m1.hi-gains[1];am2=m2.hi-gains[2]
            assert am1>=0 and am2>=0
            candidates.append((am1,'recursive_first'))
            if D:candidates.append(((I(am2)/(FOUR*d*(ONE-d))).hi,'recursive_second'))
            atoms.update(p,D)
        loss,kind=min(candidates);wins[kind]+=1;mu-=loss
        assert mu>0,('nonpositive survivor bound',i,p)
        v=ONE/(ONE-d)
        C=C*(ONE+v*I.rat(1,p-1));T=T*(ONE+v*I.rat(3*p-1,(p-1)**2));nu[p]=v
        full.update(p,D)
        if coarse is not None:coarse.update(p,D)
        if p<K and idx>=r:
            rest.update(p,D)
            if coarse_rest is not None:coarse_rest.update(p,D)
        if dynamic:small.update(p,D)
        parameters[p]=D
        if i in (100,1000,10000,25000,50000,100000) or i==n:
            row=dict(i=i,p=p,mu_lower=fixed_text(mu),T_upper=fixed_text(T.hi),
                f_upper=fixed_text(ceildiv(T.hi*S,mu)),wins=wins.copy())
            rows.append(row);print('Verified stage',i,'prime',p,flush=True)
    if finite_only:
        result=dict(status='VERIFIED_FINITE_PREFIX',cutoff=K,terminal_index=n,
            terminal_prime=ps[n-1],mu_lower=fixed_text(mu),decimal_places=40,
            witness_sha256=hashlib.sha256(raw).hexdigest(),wins=wins,
            note='Finite losses only; no infinite-tail or minimum-modulus conclusion.')
        Path(output).write_text(json.dumps(result,indent=2)+'\n');return result
    fu=ceildiv(T.hi*S,mu);stop=stopping_bound(n);tail_info={}
    stopping_method=w.get('stopping_method','bbmst')
    if stopping_method in ('third_moment','fractional_moment'):
        assert ps[n-1]>=K
        if stopping_method=='third_moment':
            from third_moment_tail import third_moment,tail_bound
            tail,tail_info=tail_bound(ps[n-1],n,third_moment(ps,start,ds),w['tail_delta_numerator'],w.get('interpolate_tail',False),w.get('shifted_tail',False))
        else:
            from fractional_moment_tail import fractional_moment,fractional_tail_bound
            num=w['moment_numerator']
            den=w.get('moment_denominator',4)
            tail,tail_info=fractional_tail_bound(ps[n-1],n,fractional_moment(ps,start,ds,num,den),w['tail_delta_numerator'],num,den)
        if 'stratified_tail_orders' in w:
            from stratified_moment_tail import tail_bound as stratified_tail_bound
            stratified,details=stratified_tail_bound(full,ps,start,ds,w['stratified_tail_orders'],w['tail_delta_numerator'],C.hi,w.get('stratified_overflow_orders'))
            tail_info.update(details)
            tail_info['single_moment_tail_loss_upper']=fixed_text(tail)
            tail=min(tail,stratified)
            tail_info['tail_loss_upper']=fixed_text(tail)
            tail_info['stratified_tail_enabled']=True
        assert mu>tail,('moment tail failed',fixed_text(mu),fixed_text(tail))
        tail_info['infinite_survivor_lower']=fixed_text(mu-tail)
    else:
        assert stopping_method=='bbmst'
        assert n>=10 and fu<stop,('analytic stopping criterion failed',n,fixed_text(fu),fixed_text(stop))
    result=dict(status='VERIFIED',cutoff=K,minimum_modulus_at_most=K-1,
        method='finite conditioning correction and first-moment product bins',decimal_places=40,
        initial_mu_lower=fixed_text(initial_mu),terminal_index=n,terminal_prime=ps[n-1],
        f_upper=fixed_text(fu),stopping_lower=fixed_text(stop),bbmst_gap=fixed_text(stop-fu),
        certified_gap=fixed_text(stop-fu if stopping_method=='bbmst' else mu-tail),
        gap_kind='bbmst stopping threshold' if stopping_method=='bbmst' else 'infinite survivor mass',
        mu_lower=fixed_text(mu),bins=len(full.bins),cap=cap,resolution=resolution,
        prefix_resolution=prefix_resolution,prefix_index=prefix_index,product_tail_power=product_tail_power,
        support_mode=support_mode,
        normalize_conditioning=normalize_conditioning,
        dynamic=dynamic,joint_presence=joint_presence,joint_improved_stages=joint_wins,
        statewise_presence=statewise_presence,
        joint_product=joint_product,
        conditioned_complement_prime=conditioned_complement_prime,
        exact_early_prime=exact_early_prime,early_exact_states=early_states,
        omitted_probability_upper=fixed_text(full.error),
        small_primes=ps[:r],conditioning_states=len(small.states),state_cap=state_cap,
        retained_conditioning_probability_lower=fixed_text(small.probability_lower),
        retained_conditioning_probability_upper=fixed_text(min(S,small.probability_upper)),
        witness_sha256=hashlib.sha256(raw).hexdigest(),wins=wins,checkpoints=rows,
        stopping_method=stopping_method,**tail_info)
    Path(output).write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness');ap.add_argument('--output',required=True)
    ap.add_argument('--finite-only',action='store_true',help='Verify finite losses only; no modulus certificate.')
    args=ap.parse_args();t=time.time();print(json.dumps(verify(args.witness,args.output,args.finite_only),indent=2))
    print('seconds',time.time()-t)
