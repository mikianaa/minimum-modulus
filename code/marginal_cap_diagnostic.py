"""Certified marginal atom caps; these are not sequential conditional caps.

The BBMST fibre density is bounded by g(alpha)=1/(1-min(alpha,delta))
and g(alpha)<=1+nu*alpha. Thus every residue modulo p**a, a>=1,
has unconditional probability at most min(nu,1+nu*M1)/p**a.
Multiplying these marginal bounds across primes is NOT justified.
"""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from verify_minimum_modulus import I,ONE,S,initial,moments,factorizations,primes_to,fixed_text
from verify_stoploss import L


if __name__ == '__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('witness');ap.add_argument('--last-prime',type=int,default=859)
    ap.add_argument('--output',required=True);a=ap.parse_args()
    raw=Path(a.witness).read_bytes();w=json.loads(raw);K=w['K'];start=w['start']
    ps=primes_to(max(K,a.last_prime+100));_,C,T=initial(K,ps,start)
    fs=factorizations((K-1)//ps[start]+1);nu={};rows=[]
    for j,D in enumerate(w['deltas'],start):
        p=ps[j]
        if p>a.last_prime:break
        m1,_=moments(K,p,C,T,nu,fs)
        v=I.rat(L,L-D);gamma=min(v.hi,(ONE+v*m1).hi)
        rows.append(dict(prime=p,delta_numerator=D,first_moment_upper=fixed_text(m1.hi),
                         sequential_cap_upper=fixed_text(v.hi),marginal_cap_upper=fixed_text(gamma)))
        C=C*(ONE+v*I.rat(1,p-1));T=T*(ONE+v*I.rat(3*p-1,(p-1)**2));nu[p]=v
    ranked=sorted(rows,key=lambda r:Fraction(r['marginal_cap_upper'])/Fraction(r['sequential_cap_upper']))
    result=dict(status='VERIFIED_LOCAL_INEQUALITIES',cutoff=K,source=a.witness,
                source_sha256=hashlib.sha256(raw).hexdigest(),decimal_places=40,
                warning='Marginal bounds cannot be multiplied across primes. This is not a modulus certificate and is not used by the current verifier.',
                strongest_relative_improvements=ranked[:10],rows=rows)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['strongest_relative_improvements'][:3],indent=2))
