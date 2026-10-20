"""Certified moment improvements, with an explicitly experimental loss comparison."""
import argparse, hashlib, json, time
from pathlib import Path
from recursive_atom_bound import RecursiveAtoms
from verify_minimum_modulus import I, ONE, S, initial, moments, primes_to, factorizations, fixed_text
from verify_stoploss import L


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('source')
    ap.add_argument('--K',type=int);ap.add_argument('--cap',type=int,default=512)
    ap.add_argument('--last-prime',type=int,default=859);ap.add_argument('--output',required=True)
    ap.add_argument('--density-degree',type=int,choices=(1,2),default=1)
    a=ap.parse_args();raw=Path(a.source).read_bytes();w=json.loads(raw)
    K=a.K or w['K'];ps=primes_to(max(K,a.last_prime+100))
    atoms=RecursiveAtoms(K,a.cap,a.density_degree);_,C,T=initial(K,ps,w['start'])
    for p in ps[:w['start']]:atoms.update(p,0)
    fs=factorizations(K);nu={};rows=[];t=time.time()
    for j,D in enumerate(w['deltas'],w['start']):
        p=ps[j]
        if p>a.last_prime:break
        m1,m2=moments(K,p,C,T,nu,fs);g=atoms.moment_gains(p)
        new1=m1.hi-g[1];new2=m2.hi-g[2];assert new1>=0 and new2>=0
        loss=min(new1,((I(new2)/(I.rat(4*D,L)*I.rat(L-D,L))).hi if D else new1))
        rows.append(dict(prime=p,first_upper=fixed_text(m1.hi),second_upper=fixed_text(m2.hi),
                         first_gain_lower=fixed_text(g[1]),second_gain_lower=fixed_text(g[2]),
                         refined_moment_loss_upper=fixed_text(loss)))
        atoms.update(p,D);v=I.rat(L,L-D);nu[p]=v
        C=C*(ONE+v*I.rat(1,p-1));T=T*(ONE+v*I.rat(3*p-1,(p-1)**2))
        if j%25==0:print('atoms',p,atoms.count,'seconds',time.time()-t,flush=True)
    result=dict(status='VERIFIED_MOMENT_BOUNDS',cutoff=K,cap=a.cap,last_prime=a.last_prime,
                source=a.source,source_sha256=hashlib.sha256(raw).hexdigest(),decimal_places=40,
                refined_atoms=atoms.count,density_degree=a.density_degree,rows=rows,
                note='Actual-measure moment bounds only; not a minimum-modulus certificate.')
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print('saved',a.output,'seconds',time.time()-t,flush=True)
