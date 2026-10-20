"""Experimental adjoint optimization of the terminal stopping margin.

No proof relies on these floating-point derivatives or on optimization
convergence. The resulting rational parameters require the integer verifier.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from sparse_conditioning_experiment import SparseStates
from minimum_modulus_experiment import sieve,initial
from binned_distribution import Binned,upper_support


class GlobalProblem:
    def stratified_objective(self,p,i):
        if i not in self.stratified_objectives:
            from stratified_tail_search import StratifiedObjective
            self.stratified_objectives[i]=StratifiedObjective(self.values,self.upper,p,i,self.tail_delta_numerator,self.stratified_tail_orders,(self.moment_numerator,self.moment_denominator))
        return self.stratified_objectives[i]

    def moment_increment(self,p):
        if hasattr(self,'fractional_increments'):return self.fractional_increments[p]
        return (7*p*p-2*p+1)/(p-1)**3

    def moment_tail_coefficient(self,p,index):
        from verify_minimum_modulus import ONE,S
        if self.stopping_method=='fractional_moment':
            from fractional_moment_tail import fractional_tail_bound
            return fractional_tail_bound(p,index,ONE,self.tail_delta_numerator,self.moment_numerator,getattr(self,'moment_denominator',4))[0]/S
        from third_moment_tail import tail_bound
        return tail_bound(p,index,ONE,getattr(self,'tail_delta_numerator',450000000000),getattr(self,'interpolate_tail',False),getattr(self,'shifted_tail',False))[0]/S

    def __init__(self,witness,K=None,terminal=10000,control_prime=20000,knots=12,
                 resolution=None,state_cap=None,support_mode=None,stopping_method=None,
                 normalize_conditioning=False):
        self.K=K or witness['K'];self.start=witness['start'];self.small=witness['small']
        self.cap=witness['cap'];self.resolution=resolution or witness['resolution']
        self.state_cap=state_cap or witness['state_cap'];self.terminal=terminal
        self.support_mode=support_mode or witness.get('support_mode','uniform')
        self.stopping_method=stopping_method or witness.get('stopping_method','bbmst')
        self.normalize_conditioning=normalize_conditioning or witness.get('normalize_conditioning',False)
        ps=sieve(max(2000000,20*terminal))
        self.primes=ps[self.start:terminal]
        assert len(witness['deltas'])>=len(self.primes)
        self.seed=np.asarray(witness['deltas'][:len(self.primes)],dtype=float)/10**12
        self.seed=np.clip(self.seed,1e-6,.5-1e-6)
        self.logits=np.log(2*self.seed/(1-2*self.seed))
        self.control_count=sum(p<=control_prime for p in self.primes)
        self.basis=np.zeros((self.control_count,knots))
        grid=np.linspace(math.log(self.primes[0]),math.log(self.primes[self.control_count-1]),knots)
        for j in range(knots):self.basis[:,j]=np.interp(np.log(self.primes[:self.control_count]),grid,np.eye(knots)[j])
        self.sm=SparseStates(self.K,self.primes[0],ps[:self.small],self.state_cap,self.normalize_conditioning)
        self.bins=Binned(self.cap,self.resolution);self.values=self.bins.values.astype(float)
        self.upper=upper_support(self.cap,self.resolution).astype(float)
        if self.support_mode=='reachable':
            from reachable_support import reachable_upper_support
            self.upper=np.asarray(reachable_upper_support(tuple(map(int,self.values))),dtype=float)
        self.width=self.upper-self.values;self.invwidth=np.divide(1.,self.width,out=np.zeros_like(self.width),where=self.width>0)
        self.mu0,self.C0,self.T0=initial(self.K,ps,self.start)
        self.U0=math.prod(1+(7*p*p-2*p+1)/(p-1)**3 for p in ps[:self.start])
        if self.stopping_method=='third_moment':
            from third_moment_tail import tail_bound
            from verify_minimum_modulus import ONE,S
            self.tail_coefficient=tail_bound(self.primes[-1],terminal,ONE)[0]/S
        P=np.zeros(len(self.values));P[1]=1.;M=P.copy();R=P.copy()
        for j,q in enumerate(ps[:self.start]):
            P,M=self.bins.update_pair(P,M,q,1.)
            if j>=self.small:R=self.bins.update(R,q,1.)
        self.initial=(P,M,R)

    def controls(self,x):
        ds=self.seed.copy();z=self.logits[:self.control_count]+self.basis@x
        ds[:self.control_count]=.5/(1+np.exp(-np.clip(z,-35,35)))
        return ds

    def coefficients(self,p,d,derivative=False):
        A=(p-1)*d;t=min(A,self.cap)
        inside=self.upper<=t;cross=(self.values<t)&~inside
        ratio=np.where(cross,(t-self.values)*self.invwidth,0.)
        cp=np.where(inside,t,self.upper*ratio);cm=np.where(inside,-1.,-ratio)
        h=np.zeros(len(cp));indicator=h.copy()
        if p<self.K:
            tau,B,weights=self.sm.groups(self.K,p)
            j=np.searchsorted(self.values,A/tau,side='right')
            k=np.searchsorted(self.values,(A+B)/tau,side='right')
            size=len(cp)+1
            slope=np.bincount(j,weights=weights*tau,minlength=size)-np.bincount(k,weights=weights*tau,minlength=size)
            intercept=np.bincount(j,weights=-weights*A,minlength=size)+np.bincount(k,weights=weights*(A+B),minlength=size)
            h=self.values*np.cumsum(slope)[:-1]+np.cumsum(intercept)[:-1]
            # Residual roundoff in a mathematically nonnegative coefficient.
            h=np.maximum(0.,h)
            if derivative:
                indicator=np.cumsum(np.bincount(j,weights=weights,minlength=size)-np.bincount(k,weights=weights,minlength=size))[:-1]
        if not derivative:return cp,cm,h
        dcp=np.where(inside,1.,np.where(cross,self.upper*self.invwidth,0.))
        dcm=np.where(cross,-self.invwidth,0.)
        if A>=self.cap:dcp*=0;dcm*=0
        return cp,cm,h,dcp,dcm,indicator,float(A<self.cap)

    def transpose(self,adj,q,nu,moment=False):
        out=(1-nu/q)*adj;weight=nu*(q-1)/q**2;h=2
        while weight>1e-30:
            keep,indices=self.bins.maps[h]
            out[keep]+=weight*(h if moment else 1)*adj[indices]
            weight/=q;h+=1
        return out

    def evaluate(self,x,gradient=False,save=False):
        ds=self.controls(np.asarray(x));P,M,R=(a.copy() for a in self.initial)
        mu=self.mu0;C=self.C0;T=self.T0;U=self.U0;N=self.control_count;B=len(P)
        # Checkpoint the forward state, then reconstruct short blocks during
        # the reverse pass. Memory no longer grows as (primes x bins).
        checkpoints={};block_size=256
        min_mu=mu;terminal_checkpoints=[]
        weight_potentials={}
        for j,(p,d) in enumerate(zip(self.primes,ds)):
            if gradient and j<N and j%block_size==0:
                checkpoints[j]=(P.copy(),M.copy(),R.copy(),C)
            v=1/(1-d);a=(3*p-1)/(p-1)**2
            cp,cm,h=self.coefficients(p,d)
            numerator=C-min((p-1)*d,self.cap)+np.sum(cp*P+cm*M-h*R)
            loss=v/(p-1)*numerator
            if gradient and getattr(self,'track_conditioning_weights',False) and p<self.K:
                tau,bb,weights=self.sm.groups(self.K,p)
                a0=(p-1)*d
                rp=np.cumsum(R);rw=np.cumsum(R*self.values)
                jj=np.searchsorted(self.values,a0/tau,side='right')-1
                kk=np.searchsorted(self.values,(a0+bb)/tau,side='right')-1
                gain=tau*(rw[kk]-rw[jj])-a0*(rp[kk]-rp[jj])+bb*(rp[-1]-rp[kk])
                powers=[];pw=p
                while pw<self.K:powers.append(pw);pw*=p
                key=(min(p,self.sm.last+1),tuple((self.K-1)//pw+1 for pw in powers))
                if key not in weight_potentials:weight_potentials[key]=np.zeros_like(gain)
                weight_potentials[key]-=v/(p-1)*gain
            mu-=loss;min_mu=min(min_mu,mu)
            Pnew,Mnew=self.bins.update_pair(P,M,p,v)
            Rnew=self.bins.update(R,p,v) if p<self.K else R
            P,M,R=Pnew,Mnew,Rnew;C*=1+v/(p-1);T*=1+a*v
            U*=1+v*self.moment_increment(p)
            i=self.start+j+1
            if self.stopping_method in ('third_moment','fractional_moment') and i in (50000,75000,100000,125000,150000,200000):
                coefficient=self.moment_tail_coefficient(p,i)
                checkpoint_tail=U*coefficient
                excluded=[]
                if getattr(self,'search_stratified_tail',False):
                    stratified=self.stratified_objective(p,i)
                    checkpoint_tail=stratified.evaluate(P,M,U)[0]
                    excluded=stratified.excluded_orders
                terminal_checkpoints.append(dict(i=i,p=p,mu=float(mu),moment_value=float(U),tail=float(checkpoint_tail),margin=float(mu-checkpoint_tail),excluded_stratified_orders=excluded))
        stop=self.terminal*(math.log(self.terminal)+math.log(math.log(self.terminal))-3)**2
        terminal_cost=U*self.tail_coefficient if self.stopping_method in ('third_moment','fractional_moment') else T/stop
        moment_terminal_cost=terminal_cost
        terminal_AP=np.zeros(B);terminal_AM=np.zeros(B)
        if getattr(self,'search_stratified_tail',False):
            stratified=self.stratified_objective(self.primes[-1],self.terminal)
            assert not stratified.excluded_orders, 'A requested terminal order is inadmissible.'
            terminal_cost,terminal_AP,terminal_AM,moment_terminal_cost=stratified.evaluate(P,M,U)
        objective=-mu+terminal_cost
        self.terminal_distribution=(P,M)
        result=dict(status='EXPERIMENTAL',K=self.K,start=self.start,small=self.small,state_cap=self.state_cap,
            cap=self.cap,resolution=self.resolution,terminal_index=self.terminal,
            support_mode=self.support_mode,
            normalize_conditioning=self.normalize_conditioning,
            stopping_method=self.stopping_method,third_moment=U,terminal_cost=terminal_cost,
            terminal_prime=self.primes[-1],mu=mu,T=T,objective=objective,
            f=T/mu if mu>0 else None,stop=stop,
            x=list(map(float,x)),minimum_mu=min_mu)
        result['terminal_checkpoints']=terminal_checkpoints
        result['search_tail_weight']=getattr(self,'search_tail_weight',1.)
        result['prime_preconditioning']=getattr(self,'prime_preconditioning',False)
        result['search_stratified_tail']=getattr(self,'search_stratified_tail',False)
        if result['search_stratified_tail']:result['stratified_tail_orders']=self.stratified_tail_orders
        if getattr(self,'stratified_overflow_orders',None) is not None:result['stratified_overflow_orders']=self.stratified_overflow_orders
        if result['search_tail_weight']!=1:
            result['objective_note']='Experimental weighted tail proxy; actual tail must be independently evaluated and certified.'
        if save:result['deltas']=[round(d*10**12) for d in ds]
        if self.stopping_method=='third_moment':result['tail_delta_numerator']=450000000000
        if self.stopping_method=='fractional_moment':
            result['fractional_moment']=result.pop('third_moment')
            result['moment_numerator']=self.moment_numerator
            result['moment_denominator']=getattr(self,'moment_denominator',4)
        if not gradient:return objective,result
        AP=terminal_AP.copy();AM=terminal_AM.copy();AR=np.zeros(B);AC=0.;gd=np.zeros(N);loaded=-1
        for j in range(len(ds)-1,-1,-1):
            p=self.primes[j];d=ds[j];v=1/(1-d);a=(3*p-1)/(p-1)**2
            if j<N:
                lo=j//block_size*block_size
                if lo!=loaded:
                    loaded=lo;hi=min(N,lo+block_size)
                    P,M,R,C=checkpoints[lo]
                    tangents=np.zeros((hi-lo,3,B));ctangents=np.zeros(hi-lo);direct=np.zeros(hi-lo)
                    for jj in range(lo,hi):
                        q=self.primes[jj];dd=ds[jj];vv=1/(1-dd);local=jj-lo
                        cp,cm,h,dcp,dcm,ind,unclipped=self.coefficients(q,dd,True)
                        loss=vv/(q-1)*(C-min((q-1)*dd,self.cap)+np.sum(cp*P+cm*M-h*R))
                        direct[local]=vv*loss+vv*(-unclipped+np.sum(dcp*P+dcm*M+ind*R))
                        ctangents[local]=C/(q-1)
                        Pnew,Mnew=self.bins.update_pair(P,M,q,vv)
                        Rnew=self.bins.update(R,q,vv) if q<self.K else R
                        tangents[local,0]=(Pnew-P)/vv;tangents[local,1]=(Mnew-M)/vv;tangents[local,2]=(Rnew-R)/vv
                        P,M,R=Pnew,Mnew,Rnew;C*=1+vv/(q-1)
                local=j-lo
                gd[j]=direct[local]+v*v*(np.sum(AP*tangents[local,0]+AM*tangents[local,1]+AR*tangents[local,2])+AC*ctangents[local])
                aa=self.moment_increment(p) if self.stopping_method in ('third_moment','fractional_moment') else a
                gd[j]+=moment_terminal_cost*aa*v*v/(1+aa*v)
            cp,cm,h=self.coefficients(p,d)
            AP=self.transpose(AP,p,v)+v/(p-1)*cp
            AM=self.transpose(AM,p,v,True)+v/(p-1)*cm
            if p<self.K:AR=self.transpose(AR,p,v)-v/(p-1)*h
            AC=AC*(1+v/(p-1))+v/(p-1)
        grad=self.basis.T@(gd*ds[:N]*(1-2*ds[:N]))
        if getattr(self,'track_conditioning_weights',False):
            potential=np.zeros(len(self.sm.ns))
            for key,gg in weight_potentials.items():
                _,inv=self.sm.cache[key];potential+=gg[inv]
            self.prefix_sensitivity=(AP,AM,AR,AC,moment_terminal_cost/self.U0,potential)
        return objective,grad,result


def optimize(problem,iterations,output,x=None,initial_scale=100,max_step=1.5):
    n=problem.basis.shape[1];x=np.zeros(n) if x is None else np.asarray(x,dtype=float)
    diagonal=np.full(n,initial_scale)
    if getattr(problem,'prime_preconditioning',False):
        diagonal*=np.geomspace(problem.primes[0],problem.primes[problem.control_count-1],n)/problem.primes[0]
    H0=np.diag(diagonal);H=H0.copy();history=[];f,g,r=problem.evaluate(x,True)
    for it in range(iterations+1):
        r.update(iteration=it,gradient_norm=float(np.linalg.norm(g)))
        history.append(r);print(it,r['mu'],r['f'],r['gradient_norm'],flush=True)
        _,record=problem.evaluate(x,save=True);record['history']=history
        Path(output).write_text(json.dumps(record))
        if it==iterations:break
        direction=-H@g
        if direction@g>=0:H=H0.copy();direction=-H@g
        direction*=min(1.,max_step/max(np.linalg.norm(direction),1e-30))
        alpha=1.;accepted=False
        for _ in range(10):
            ff,_=problem.evaluate(x+alpha*direction)
            if ff<f+1e-4*alpha*(g@direction):accepted=True;break
            alpha*=.5
        if not accepted:break
        xx=x+alpha*direction;ff,gg,rr=problem.evaluate(xx,True)
        s=xx-x;y=gg-g;sy=s@y
        if sy>1e-12:
            rho=1/sy;V=np.eye(n)-rho*np.outer(s,y);H=V@H@V.T+rho*np.outer(s,s)
        else:H=H0.copy()
        x,f,g,r=xx,ff,gg,rr
    return record


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness');ap.add_argument('--K',type=int)
    ap.add_argument('--terminal',type=int,default=10000);ap.add_argument('--control-prime',type=int,default=20000)
    ap.add_argument('--knots',type=int,default=12);ap.add_argument('--iterations',type=int,default=15)
    ap.add_argument('--resolution',type=int);ap.add_argument('--state-cap',type=int)
    ap.add_argument('--support-mode',choices=['uniform','reachable'])
    ap.add_argument('--stopping-method',choices=['bbmst','third_moment'])
    ap.add_argument('--normalize-conditioning',action='store_true')
    ap.add_argument('--output',required=True);ap.add_argument('--initial')
    args=vars(ap.parse_args());w=json.loads(Path(args.pop('witness')).read_text())
    out=args.pop('output');iterations=args.pop('iterations');init=args.pop('initial')
    x=json.loads(Path(init).read_text())['x'] if init else None
    t=time.time();problem=GlobalProblem(w,**args);optimize(problem,iterations,out,x)
    print('seconds',time.time()-t)
