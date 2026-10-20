"""Adjoint optimization after a fixed, dynamically conditioned prefix.

The early parameters are fixed inputs, not silently differentiated through
the joint presence calculation. All subsequent parameters are optimized
together against the third-moment stopping margin.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from global_distortion_search import GlobalProblem,optimize
from sparse_conditioning_experiment import DynamicSparseStates
from bernoulli_conditioning_experiment import JointStates
from minimum_modulus_experiment import sieve,initial,moments,factors_table
from binned_distribution import Binned,MomentPrefixes
from reachable_support import reachable_upper_support
from third_moment_tail import tail_bound
from verify_minimum_modulus import ONE,S


class FrozenPrefixProblem(GlobalProblem):
    def __init__(self,witness,K=None,terminal=50000,freeze=100,knots=20,
                 state_cap=None,resolution=None,tail_delta_numerator=493000000000):
        self.K=K or witness['K'];self.original_start=witness['start'];self.small=witness['small']
        self.start=freeze;self.terminal=terminal;self.cap=witness['cap']
        self.state_cap=state_cap or witness['state_cap'];self.resolution=resolution or witness['resolution']
        self.support_mode='reachable';self.stopping_method='third_moment';self.normalize_conditioning=False
        if witness.get('stopping_method')=='fractional_moment':self.stopping_method='fractional_moment'
        self.tail_delta_numerator=tail_delta_numerator
        self.interpolate_tail=True
        self.statewise_presence=witness.get('statewise_presence',False)
        self.joint_product=witness.get('joint_product',False)
        self.shifted_tail=witness.get('shifted_tail',False)
        self.exact_early_prime=witness.get('exact_early_prime',0)
        self.conditioned_complement_prime=witness.get('conditioned_complement_prime',0)
        self.recursive_atom_cap=witness.get('recursive_atom_cap',0)
        self.recursive_atom_prime=witness.get('recursive_atom_prime',0)
        self.recursive_density_degree=witness.get('recursive_density_degree',1)
        from early_exact_experiment import EarlyExact
        early_exact=EarlyExact(self.K)
        ps=sieve(max(2000000,20*terminal));assert self.original_start<=self.small<freeze<terminal
        assert self.conditioned_complement_prime<ps[freeze]
        if self.stopping_method=='fractional_moment':
            from fractional_moment_tail import increment
            self.moment_numerator=witness['moment_numerator']
            self.moment_denominator=witness.get('moment_denominator',4)
            self.fractional_increments={q:increment(q,self.moment_numerator,self.moment_denominator).hi/S for q in ps[:terminal]}
        assert ps[freeze]*ps[self.small]>=self.K
        raw=list(witness['deltas']);assert len(raw)>=freeze-self.original_start
        raw.extend([400000000000]*max(0,terminal-self.original_start-len(raw)))
        self.prefix=raw[:freeze-self.original_start]
        self.primes=ps[freeze:terminal];self.control_count=len(self.primes)
        self.seed=np.clip(np.asarray(raw[freeze-self.original_start:terminal-self.original_start])/1e12,1e-6,.5-1e-6)
        self.logits=np.log(2*self.seed/(1-2*self.seed));self.basis=np.zeros((len(self.primes),knots))
        grid=np.linspace(math.log(self.primes[0]),math.log(self.primes[-1]),knots)
        for j in range(knots):self.basis[:,j]=np.interp(np.log(self.primes),grid,np.eye(knots)[j])
        self.bins=Binned(self.cap,self.resolution);self.values=self.bins.values.astype(float)
        self.upper=np.asarray(reachable_upper_support(tuple(map(int,self.values))),dtype=float)
        self.width=self.upper-self.values;self.invwidth=np.divide(1.,self.width,out=np.zeros_like(self.width),where=self.width>0)
        self.sm=DynamicSparseStates(self.K,ps,self.original_start,self.small,self.state_cap)
        P=np.zeros(len(self.values));P[1]=1.;M=P.copy();R=P.copy()
        for q in ps[:self.original_start]:P,M=self.bins.update_pair(P,M,q,1.)
        mu,C,T=initial(self.K,ps,self.original_start)
        U=math.prod(1+self.moment_increment(p) for p in ps[:self.original_start])
        fs=factors_table((self.K-1)//ps[self.original_start]+1);nus={};joint=None
        self.prefix_loss_rows=[]
        atoms=self.new_recursive_atoms(ps)
        for idx in range(self.original_start,freeze):
            p=ps[idx];d=raw[idx-self.original_start]/1e12;v=1/(1-d);A=(p-1)*d;t=min(A,self.cap)
            if idx==self.small:
                self.sm.prob=self.sm.normalizer*self.sm.coeff
                joint=JointStates(self.sm,self.K,p,p)
            tau,B,weights=self.sm.groups(self.K,p);rp=np.cumsum(R);rw=np.cumsum(R*self.values)
            j=np.searchsorted(self.values,A/tau,side='right')-1;k=np.searchsorted(self.values,(A+B)/tau,side='right')-1
            gain=tau*(rw[k]-rw[j])-A*(rp[k]-rp[j])+B*(rp[-1]-rp[k])
            correction=float(np.sum(weights*np.maximum(0.,gain)))
            if joint is not None and p*ps[self.small]<self.K and (self.K-1)//p<ps[self.small]**2 and p*p*ps[self.small]>=self.K:
                correction=max(correction,joint.correction(p,d,ps[self.small:idx],[nus[q] for q in ps[self.small:idx]],R if self.statewise_presence else None,self.values))
                if self.joint_product:
                    from joint_product_experiment import correction as product_correction
                    correction=max(correction,product_correction(joint,p,d,ps[self.small:idx],[nus[q] for q in ps[self.small:idx]]))
            mp=MomentPrefixes(self.bins.values,P,M,self.resolution,'reachable')
            loss=(C-t+float(mp.complement(np.asarray([t]))[0])-correction)/((p-1)*(1-d))
            if p<self.K:
                m1=moments(self.K,p,C,T,nus,fs,first=True)[0];loss=min(loss,m1)
                if p<=self.exact_early_prime:loss=min(loss,early_exact.loss(p,d,ps[:idx],nus,m1))
                if p<=self.conditioned_complement_prime and idx>=self.small and (self.K-1)//p<ps[self.small]**2 and p*p*ps[self.small]>=self.K:
                    from joint_product_experiment import conditioned_loss
                    loss=min(loss,conditioned_loss(joint,p,d,ps[self.small:idx],[nus[q] for q in ps[self.small:idx]],ps[:self.small],[nus.get(q,1.) for q in ps[:self.small]],self.state_cap,m1))
            if atoms is not None and p<=self.recursive_atom_prime:
                am1,am2=moments(self.K,p,C,T,nus,fs)
                loss=min(loss,atoms.loss(p,d,am1,am2))
                atoms.update(p,d)
            self.prefix_loss_rows.append(dict(prime=p,loss=loss))
            mu-=loss;P,M=self.bins.update_pair(P,M,p,v)
            if idx>=self.small:R=self.bins.update(R,p,v)
            C*=1+v/(p-1);T*=1+v*(3*p-1)/(p-1)**2;U*=1+v*self.moment_increment(p)
            nus[p]=v;self.sm.update(p,v)
        self.initial=(P,M,R);self.mu0=mu;self.C0=C;self.T0=T;self.U0=U
        self.search_tail_weight=witness.get('search_tail_weight',1.)
        assert 0<=self.search_tail_weight<=1
        self.tail_coefficient=self.search_tail_weight*self.moment_tail_coefficient(self.primes[-1],terminal)
        self.search_stratified_tail=witness.get('search_stratified_tail',False)
        self.stratified_tail_orders=witness.get('stratified_tail_orders',[])
        self.stratified_overflow_orders=witness.get('stratified_overflow_orders')
        self.stratified_objectives={}
        if self.search_stratified_tail:assert self.search_tail_weight==1 and self.stopping_method=='fractional_moment'

    def new_recursive_atoms(self,ps):
        if not self.recursive_atom_cap:return None
        from recursive_atom_experiment import FloatRecursiveAtoms
        assert 2<=self.recursive_atom_cap and 0<self.recursive_atom_prime<ps[self.start]
        atoms=FloatRecursiveAtoms(self.K,self.recursive_atom_cap,self.recursive_density_degree)
        for p in ps[:self.original_start]:atoms.update(p,0.)
        return atoms

    def evaluate(self,x,gradient=False,save=False):
        result=super().evaluate(x,gradient,save);r=result[-1]
        r.update(start=self.original_start,dynamic=True,joint_presence=True,
            frozen_index=self.start,tail_delta_numerator=self.tail_delta_numerator,interpolate_tail=True)
        r['statewise_presence']=self.statewise_presence
        r['joint_product']=self.joint_product
        r['shifted_tail']=self.shifted_tail
        r['exact_early_prime']=self.exact_early_prime
        r['conditioned_complement_prime']=self.conditioned_complement_prime
        if self.recursive_atom_cap:
            r.update(recursive_atom_cap=self.recursive_atom_cap,recursive_atom_prime=self.recursive_atom_prime,recursive_density_degree=self.recursive_density_degree)
        if save:r['deltas']=self.prefix+r['deltas']
        return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness');ap.add_argument('--K',type=int)
    ap.add_argument('--terminal',type=int,default=50000);ap.add_argument('--freeze',type=int,default=100)
    ap.add_argument('--knots',type=int,default=20);ap.add_argument('--iterations',type=int,default=12)
    ap.add_argument('--state-cap',type=int);ap.add_argument('--resolution',type=int)
    ap.add_argument('--initial-scale',type=float,default=100.)
    ap.add_argument('--max-step',type=float,default=1.5)
    ap.add_argument('--prime-preconditioning',action='store_true')
    ap.add_argument('--tail-delta-numerator',type=int,default=493000000000);ap.add_argument('--output',required=True)
    args=vars(ap.parse_args());w=json.loads(Path(args.pop('witness')).read_text())
    out=args.pop('output');it=args.pop('iterations');t=time.time()
    scale=args.pop('initial_scale');step=args.pop('max_step')
    precondition=args.pop('prime_preconditioning')
    problem=FrozenPrefixProblem(w,**args);problem.prime_preconditioning=precondition
    optimize(problem,it,out,initial_scale=scale,max_step=step);print('seconds',time.time()-t)
