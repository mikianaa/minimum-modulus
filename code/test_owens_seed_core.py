"""Exhaust the early seed's claimed holes, independently of region indexing."""
import json
from math import prod
from pathlib import Path

from graded_arrow_blocks import np
from construction_kernel import verify
from owens_seed_core import seed


def main():
    checks=[]
    for minimum in (42,50):
        data=seed(minimum,True)
        verify(data)
        ps,ds=data['primes'],data['depths']
        dims=[p**(E-1)*(p-1) for p,E in zip(ps,ds)]
        states=np.indices(dims,sparse=True)
        values=[s//(p-1)+(s%(p-1)+1)*p**(E-1) for s,p,E in zip(states,ps,ds)]
        covered=np.zeros(dims,dtype=bool)
        for row in data['rows']:
            hit=True
            for s,x,p,E,e,r in zip(states,values,ps,ds,row['exponents'],row['residues']):
                hit=hit & (s==r if e==E else (x-r)%p**e==0)
            covered|=hit
        missing=~covered
        holes=[('8 mod 16',values[0]%16==8),
               ('4 mod 8 outside 1 mod 5',(values[0]%8==4)&(values[2]%5!=1)),
               ('16 mod 32 outside 4 mod 5',(values[0]%32==16)&(values[2]%5!=4))]
        counts={name:int((missing&mask).sum()) for name,mask in holes}
        if minimum==42:
            assert all(n==0 for n in counts.values())
        checks.append(dict(minimum=minimum,rows=len(data['rows']),kernel_states=prod(dims),
                           uncovered_states=int(missing.sum()),hole_gaps=counts))
    report=dict(status='PASS',checks=checks,
                scope='Three local seed claims only; neither the full Owens thesis nor a minimum-50 cover is certified.')
    Path('data/owens-seed-branch-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
