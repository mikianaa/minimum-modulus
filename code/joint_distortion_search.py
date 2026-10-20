"""Joint early/late terminal optimization, for discovery only.

The late adjoint includes derivatives with respect to the conditioned-state
weights. A finite-difference pullback through the short early calculation
then differentiates its distributions, probabilities, losses and moments.
All saved parameters still require the independent integer verifier.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from frozen_prefix_search import FrozenPrefixProblem
from global_distortion_search import GlobalProblem
from bernoulli_conditioning_experiment import JointStates
from minimum_modulus_experiment import initial,moments,factors_table,sieve
from binned_distribution import MomentPrefixes


class JointProblem(FrozenPrefixProblem):
    def __init__(self,witness,K=None,terminal=50000,start=2,freeze=100,
                 early_knots=6,knots=16,state_cap=None,resolution=None,
                 tail_delta_numerator=500000000000):
        witness=dict(witness);oldstart=witness['start']
        self.prefix_difference=witness.get('prefix_difference','central')
        assert self.prefix_difference in ('central','forward')
        assert 0<=start<=oldstart
        witness['deltas']=[0]*(oldstart-start)+list(witness['deltas']);witness['start']=start
        super().__init__(witness,K,terminal,freeze,knots,state_cap,resolution,tail_delta_numerator)
        self.ps=sieve(max(2000000,20*terminal));self.track_conditioning_weights=True
        self.fixed_early_primes=witness.get('fixed_early_primes',[])
        self.fixed_early_indices=[self.ps.index(q)-start for q in self.fixed_early_primes]
        self.statewise_presence=True
        self.early_seed=np.asarray(self.prefix)/1e12
        self.early_scales=np.minimum(.1,np.maximum(.002,self.early_seed))
        control_small=max(self.small,witness.get('individual_prefix_count',16))
        self.control_small=control_small
        individual=control_small-start;self.individual=individual
        self.early_basis=np.zeros((freeze-start,individual+early_knots))
        self.early_basis[:individual,:individual]=np.eye(individual)
        grid=np.linspace(math.log(self.ps[control_small]),math.log(self.ps[freeze-1]),early_knots)
        for j in range(early_knots):
            self.early_basis[individual:,individual+j]=np.interp(np.log(self.ps[control_small:freeze]),grid,np.eye(early_knots)[j])
        self.early_count=self.early_basis.shape[1];self.dimension=self.early_count+knots
        self.prefix_initial=initial(self.K,self.ps,start)
        self.prefix_U=math.prod(1+self.moment_increment(q) for q in self.ps[:start])
        self.fs=factors_table((self.K-1)//self.ps[start]+1)
        p0=np.zeros(len(self.values));p0[1]=1.;m0=p0.copy()
        for q in self.ps[:start]:p0,m0=self.bins.update_pair(p0,m0,q,1.)
        self.prefix_PM=(p0,m0)
        if self.exact_early_prime>=self.ps[self.small-1]:
            from packed_conditioning import PackedConditioning
            self.sm=PackedConditioning(self.sm,self.ps[:self.small])
        self.prefix_joint=JointStates(self.sm,self.K,self.ps[self.small],self.ps[self.small])
        from early_exact_experiment import EarlyExact
        self.early_exact=EarlyExact(self.K)

    def project(self,x):
        x=np.asarray(x).copy();n=self.individual
        x[:n]=np.clip(x[:n],-self.early_seed[:n]/self.early_scales[:n],(.5-self.early_seed[:n])/self.early_scales[:n])
        if self.original_start==0:
            x[0]=min(x[0],(.25-self.early_seed[0])/self.early_scales[0])
        for j in self.fixed_early_indices:
            assert 0<=j<n
            x[j]=0.
        return x

    def set_prefix(self,x):
        ds=np.clip(self.early_seed+self.early_scales*(self.early_basis@x),0.,.5)
        if self.original_start==0:ds[0]=min(ds[0],.25)
        sm=self.sm
        if hasattr(sm,'reset'):sm.reset()
        else:sm.coeff=1./sm.ns;sm.normalizer=1.;sm.weight_cache.clear()
        for q in self.ps[:self.original_start]:sm.update(q,1.)
        P,M=(a.copy() for a in self.prefix_PM);R=np.zeros(len(P));R[1]=1.
        mu,C,T=self.prefix_initial;U=self.prefix_U;nus={}
        atoms=self.new_recursive_atoms(self.ps)
        self.prefix_joint.prepared_p=None
        for idx,d in enumerate(ds,self.original_start):
            p=self.ps[idx];v=1/(1-d);A=(p-1)*d;t=min(A,self.cap)
            if idx==self.small:sm.prob=sm.normalizer*sm.coeff
            use_complement=(p<=self.conditioned_complement_prime and idx>=self.small and (self.K-1)//p<self.ps[self.small]**2 and p*p*self.ps[self.small]>=self.K)
            if p<=self.exact_early_prime or use_complement:
                m1=moments(self.K,p,C,T,nus,self.fs,first=True)[0]
                if p<=self.exact_early_prime:
                    loss=min(m1,self.early_exact.loss(p,d,self.ps[:idx],nus,m1))
                else:
                    from joint_product_experiment import conditioned_loss
                    loss=min(m1,conditioned_loss(self.prefix_joint,p,d,self.ps[self.small:idx],[nus[q] for q in self.ps[self.small:idx]],self.ps[:self.small],[nus.get(q,1.) for q in self.ps[:self.small]],self.state_cap,m1))
                if atoms is not None and p<=self.recursive_atom_prime:
                    am1,am2=moments(self.K,p,C,T,nus,self.fs)
                    loss=min(loss,atoms.loss(p,d,am1,am2));atoms.update(p,d)
                mu-=loss;P,M=self.bins.update_pair(P,M,p,v)
                if idx>=self.small:R=self.bins.update(R,p,v)
                C*=1+v/(p-1);T*=1+v*(3*p-1)/(p-1)**2;U*=1+v*self.moment_increment(p)
                nus[p]=v;sm.update(p,v)
                continue
            tau,B,weights=sm.groups(self.K,p);rp=np.cumsum(R);rw=np.cumsum(R*self.values)
            j=np.searchsorted(self.values,A/tau,side='right')-1
            k=np.searchsorted(self.values,(A+B)/tau,side='right')-1
            gain=tau*(rw[k]-rw[j])-A*(rp[k]-rp[j])+B*(rp[-1]-rp[k])
            correction=float(weights@np.maximum(0.,gain))
            if idx>=self.small and p*self.ps[self.small]<self.K and (self.K-1)//p<self.ps[self.small]**2 and p*p*self.ps[self.small]>=self.K:
                correction=max(correction,self.prefix_joint.correction(p,d,self.ps[self.small:idx],[nus[q] for q in self.ps[self.small:idx]],R if self.statewise_presence else None,self.values))
                if self.joint_product:
                    from joint_product_experiment import correction as product_correction
                    correction=max(correction,product_correction(self.prefix_joint,p,d,self.ps[self.small:idx],[nus[q] for q in self.ps[self.small:idx]]))
            mp=MomentPrefixes(self.bins.values,P,M,self.resolution,'reachable')
            loss=(C-t+float(mp.complement(np.asarray([t]))[0])-correction)/((p-1)*(1-d))
            if p<self.K:
                m1=moments(self.K,p,C,T,nus,self.fs,first=True)[0];loss=min(loss,m1)
                if p<=self.exact_early_prime:loss=min(loss,self.early_exact.loss(p,d,self.ps[:idx],nus,m1))
                if p<=self.conditioned_complement_prime and idx>=self.small and (self.K-1)//p<self.ps[self.small]**2 and p*p*self.ps[self.small]>=self.K:
                    from joint_product_experiment import conditioned_loss
                    loss=min(loss,conditioned_loss(self.prefix_joint,p,d,self.ps[self.small:idx],[nus[q] for q in self.ps[self.small:idx]],self.ps[:self.small],[nus.get(q,1.) for q in self.ps[:self.small]],self.state_cap,m1))
            if atoms is not None and p<=self.recursive_atom_prime:
                am1,am2=moments(self.K,p,C,T,nus,self.fs)
                loss=min(loss,atoms.loss(p,d,am1,am2));atoms.update(p,d)
            mu-=loss;P,M=self.bins.update_pair(P,M,p,v)
            if idx>=self.small:R=self.bins.update(R,p,v)
            C*=1+v/(p-1);T*=1+v*(3*p-1)/(p-1)**2;U*=1+v*self.moment_increment(p)
            nus[p]=v;sm.update(p,v)
        self.initial=(P,M,R);self.mu0=mu;self.C0=C;self.T0=T;self.U0=U
        self.prefix=[round(d*1e12) for d in ds]
        sm.prob=sm.normalizer*sm.coeff

    def evaluate(self,x,gradient=False,save=False):
        x=np.asarray(x);ne=self.early_count;self.set_prefix(x[:ne])
        result=GlobalProblem.evaluate(self,x[ne:],gradient,save)
        record=result[-1]
        record.update(start=self.original_start,dynamic=True,joint_presence=True,
            joint_early_late=True,early_controls=ne,tail_delta_numerator=self.tail_delta_numerator,
            interpolate_tail=True,x=list(map(float,x)))
        record['statewise_presence']=self.statewise_presence
        record['joint_product']=self.joint_product
        record['shifted_tail']=self.shifted_tail
        record['exact_early_prime']=self.exact_early_prime
        record['conditioned_complement_prime']=self.conditioned_complement_prime
        record['prefix_complement_only_search']=bool(self.conditioned_complement_prime)
        record['prefix_difference']=self.prefix_difference
        record['individual_prefix_count']=self.control_small
        record['fixed_early_primes']=self.fixed_early_primes
        if self.recursive_atom_cap:
            record.update(recursive_atom_cap=self.recursive_atom_cap,recursive_atom_prime=self.recursive_atom_prime,recursive_density_degree=self.recursive_density_degree)
        if save:record['deltas']=self.prefix+record['deltas']
        if not gradient:return result[0],record
        AP,AM,AR,AC,AU,potential=self.prefix_sensitivity
        self.set_prefix(x[:ne])
        P0,M0,R0=(a.copy() for a in self.initial)
        mu0,C0,U0,prob0=self.mu0,self.C0,self.U0,self.sm.prob.copy()
        def surrogate(early):
            self.set_prefix(early);P,M,R=self.initial
            return -(self.mu0-mu0)+AP@(P-P0)+AM@(M-M0)+AR@(R-R0)+AC*(self.C0-C0)+AU*(self.U0-U0)+potential@(self.sm.prob-prob0)
        h=1e-4 if self.prefix_difference=='forward' else 2e-5;g=np.empty(ne)
        if self.prefix_difference=='forward':base=surrogate(x[:ne])
        for j in range(ne):
            if j in self.fixed_early_indices:g[j]=0.;continue
            e=np.zeros(ne);e[j]=h
            if self.prefix_difference=='forward':g[j]=(surrogate(x[:ne]+e)-base)/h
            else:g[j]=(surrogate(x[:ne]+e)-surrogate(x[:ne]-e))/(2*h)
        return result[0],np.r_[g,result[1]],record


def optimize_joint(problem,iterations,output,initial_scale=1000.,max_step=1.5):
    n=problem.dimension;x=np.zeros(n);diagonal=np.full(n,initial_scale)
    if getattr(problem,'prime_preconditioning',False):
        count=n-problem.early_count
        diagonal[problem.early_count:]*=np.geomspace(problem.primes[0],problem.primes[-1],count)/problem.primes[0]
    H0=np.diag(diagonal);H=H0.copy();history=[]
    f,g,r=problem.evaluate(x,True)
    for it in range(iterations+1):
        _,record=problem.evaluate(x,save=True)
        projected_g=x-problem.project(x-g)
        history.append(dict(iteration=it,objective=f,mu=r['mu'],gradient_norm=float(np.linalg.norm(g)),projected_gradient_norm=float(np.linalg.norm(projected_g)),gradient=g.tolist(),x=x.tolist()))
        record.update(history=history,optimizer_initial_scale=initial_scale,optimizer_max_step=max_step)
        Path(output).write_text(json.dumps(record)+'\n')
        print(it,'margin',-f,'mu',r['mu'],'projected gradient',np.linalg.norm(projected_g),flush=True)
        if it==iterations:break
        direction=-H@projected_g
        # Remove outward directions at active bounds before the descent check.
        direction=problem.project(x+direction)-x
        if direction@g>=0:
            H=H0.copy()
            direction=problem.project(x-H0@projected_g)-x
        direction*=min(1.,max_step/max(np.linalg.norm(direction),1e-30))
        alpha=1.;accepted=False
        for _ in range(10):
            xx=problem.project(x+alpha*direction);ff,_=problem.evaluate(xx)
            if ff<f+1e-4*(g@(xx-x)):accepted=True;break
            alpha*=.5
        if not accepted:break
        ff,gg,rr=problem.evaluate(xx,True);s=xx-x;y=gg-g;sy=s@y
        if sy>1e-12:
            rho=1/sy;V=np.eye(n)-rho*np.outer(s,y);H=V@H@V.T+rho*np.outer(s,s)
        else:H=H0.copy()
        x,f,g,r=xx,ff,gg,rr
    return record


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness');ap.add_argument('--K',type=int)
    for name,default in [('terminal',50000),('start',2),('freeze',100),('early_knots',6),('knots',16),('iterations',12)]:
        ap.add_argument('--'+name.replace('_','-'),type=int,default=default)
    ap.add_argument('--state-cap',type=int);ap.add_argument('--resolution',type=int)
    ap.add_argument('--initial-scale',type=float,default=1000.)
    ap.add_argument('--max-step',type=float,default=1.5)
    ap.add_argument('--prime-preconditioning',action='store_true')
    ap.add_argument('--output',required=True);args=vars(ap.parse_args())
    witness=json.loads(Path(args.pop('witness')).read_text());out=args.pop('output');it=args.pop('iterations')
    scale=args.pop('initial_scale');step=args.pop('max_step')
    precondition=args.pop('prime_preconditioning')
    t=time.time();problem=JointProblem(witness,**args);problem.prime_preconditioning=precondition
    optimize_joint(problem,it,out,scale,step)
    print('seconds',time.time()-t)
