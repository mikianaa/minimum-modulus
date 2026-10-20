"""Experimental mixtures of actual moment bounds and comparison stop-loss.

The pointwise inequalities are valid, but this evaluator uses floating
comparison bounds and cannot establish a modulus certificate.
"""
import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction as F
from minimum_modulus_experiment import sieve
from sparse_conditioning_experiment import DynamicSparseStates
from bernoulli_conditioning_experiment import JointStates
from early_exact_experiment import EarlyExact
from joint_product_experiment import conditioned_loss


def check_inequalities():
    count=0
    for di in (1,3,7,10):
        d=F(di,20)
        for j in (1,3,7,9):
            c=F(j,10);t=d/(1-c)
            u=d*(1+F(j,10))
            for xi in range(41):
                x=F(xi,20);loss=max(0,x-d)
                assert loss<=c*x+(1-c)*max(0,x-t)
                assert loss<=(u-d)*x*x/(u*u)+(2*d/u-1)*max(0,x-u)
                count+=2
    return count


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('profile');ap.add_argument('--output',required=True)
    a=ap.parse_args();raw=Path(a.profile).read_bytes();r=json.loads(raw)
    source=Path(r['source']).read_bytes();assert hashlib.sha256(source).hexdigest()==r['source_sha256']
    w=json.loads(source);K=r['cutoff'];ps=sieve(K);start=w['start'];small=w['small']
    sm=DynamicSparseStates(K,ps,start,small,w['state_cap'])
    joint=JointStates(sm,K,ps[small],ps[small]);early=EarlyExact(K);nus={};rows=[]
    wanted={31,53,79,101,151,251,557,829,859};profiles={v['prime']:v for v in r['rows']}
    checks=check_inequalities()
    for idx,p in enumerate(ps[start:],start):
        if p>max(wanted):break
        d=w['deltas'][idx-start]/1e12
        if p in wanted:
            row=profiles[p];old1=float(row['first_upper'])
            m1=old1-float(row['first_gain_lower'])
            m2=float(row['second_upper'])-float(row['second_gain_lower'])
            sm.prob=sm.normalizer*sm.coeff
            def stop(t):
                assert 0<=t<1
                if p<=w['exact_early_prime']:
                    return (1-t)*early.loss(p,t,ps[:idx],nus,old1)
                return (1-t)*conditioned_loss(joint,p,t,ps[small:idx],
                      [nus[q] for q in ps[small:idx]],ps[:small],
                      [nus.get(q,1.) for q in ps[:small]],w['state_cap'],old1)
            candidates=[dict(kind='comparison',parameter=d,loss=stop(d)/(1-d))]
            for c in (.01,.025,.05,.1,.2,.3,.4):
                t=d/(1-c)
                if t<1:
                    candidates.append(dict(kind='first_mixture',parameter=c,
                                           loss=(c*m1+(1-c)*stop(t))/(1-d)))
            for fraction in (.01,.025,.05,.1,.2,.3,.4,.6,.8,1.):
                t=d*(1+fraction)
                if t<1:
                    candidates.append(dict(kind='second_mixture',parameter=t,
                         loss=((t-d)*m2/(t*t)+(2*d/t-1)*stop(t))/(1-d)))
            best=min(candidates,key=lambda z:z['loss'])
            previous=row['existing_loss_experimental']
            rows.append(dict(prime=p,delta=d,existing_loss=previous,best=best,
                             improvement=max(0.,previous-best['loss']),candidates=candidates))
            print(p,best,flush=True)
        nus[p]=1/(1-d);sm.update(p,nus[p])
    result=dict(status='EXPERIMENTAL',profile=a.profile,profile_sha256=hashlib.sha256(raw).hexdigest(),
                pointwise_rational_checks=checks,rows=rows,
                total_sampled_improvement=sum(row['improvement'] for row in rows),
                note='Finite threshold grids at sampled primes; no modulus conclusion.')
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
