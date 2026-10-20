"""Owens-inspired literal early branches, with all moduli below 50 removed.

This is an explicit seed, not a transcription or verification of the
entire thesis. The optional prime-7 stage has explicit audited placements;
later stages are not imported. Unused resources are placed by exact marginals.
"""
from itertools import product
import json
from math import prod
from pathlib import Path

from exact_core_array import build


def seed(minimum=50,include_seven=False):
    ps,ds=[2,3,5,7],[6,4,3,2]
    families=[]
    F=lambda e,r:('fixed',e,r,0)
    T=lambda e,r,d=1:('tail',e,r,d)
    def add(**kw):
        families.append({int(p[1:]):v for p,v in kw.items()})
    # Pure binary tail; the two odd binary classes in the ternary arrow.
    add(p2=T(6,0))
    add(p2=F(1,1),p3=T(1,0,1))
    add(p2=T(2,1),p3=T(1,0,2))
    add(p3=T(4,21,1))
    holes=[F(2,2),F(3,4),F(4,8),F(5,16)]
    for binary in holes[2:]:
        add(p2=binary,p5=F(1,1))
    for binary in holes[:3]:
        add(p2=binary,p3=T(2,0,1),p5=F(1,2))
    add(p2=T(5,8),p3=T(2,0,2),p5=F(1,2))
    add(p2=T(6,16),p5=F(1,2))
    add(p2=T(6,16),p3=F(1,1),p5=F(1,3))
    for binary in holes:
        add(p2=binary,p3=F(1,2),p5=F(1,3))
    add(p3=T(2,0,1),p5=F(1,3))
    add(p2=F(1,0),p3=T(2,0,2),p5=F(1,3))
    add(p2=F(1,0),p5=F(2,5))
    for binary in holes:
        add(p2=binary,p5=F(2,10))
    add(p3=T(1,0,1),p5=F(2,15))
    add(p2=F(1,0),p3=T(1,0,2),p5=F(2,15))
    add(p2=T(5,8),p3=T(1,0,1),p5=F(2,20))
    for binary in holes[:3]:
        add(p2=binary,p3=T(1,0,2),p5=F(2,20))
    add(p2=T(6,16),p5=F(2,20))
    add(p5=T(3,0,1))
    add(p2=F(1,0),p5=T(3,0,2))
    add(p3=T(1,0,1),p5=T(3,0,3))
    add(p2=F(1,0),p3=T(1,0,2),p5=T(3,0,3))
    add(p2=T(2,0),p5=T(3,0,4))
    if include_seven:
        # Explicit placements for the spare 125-tail/3-tail families.
        for digit,binary in [(1,F(2,0)),(2,F(3,4)),(3,T(4,4))]:
            add(p2=binary,p3=T(1,0,1),p5=T(3,20,digit))
        def seven(d,**kw):
            add(p7=T(1,0,d),**kw)
        for binary in [F(3,4),F(4,8)]:
            seven(1,p2=binary)
        seven(2,p2=F(3,4),p3=T(1,0,1))
        seven(2,p2=T(4,4),p3=T(1,0,2))
        seven(2,p2=T(5,8))
        seven(3,p2=F(1,0),p3=F(1,1))
        seven(3,p2=F(2,0),p3=F(1,2))
        seven(3,p3=T(2,0,1))
        seven(3,p2=F(1,0),p3=T(2,0,2))
        seven(4,p2=F(2,0),p3=F(1,1),p5=F(1,2))
        for binary in [F(3,4),F(4,8)]:
            seven(4,p2=binary,p3=F(1,2),p5=F(1,2))
        seven(4,p2=F(2,0),p3=T(2,0,2),p5=F(1,2))
        seven(4,p3=F(1,1),p5=F(1,3))
        seven(4,p2=F(1,0),p5=F(1,4))
        for binary in [F(3,4),F(4,8)]:
            seven(5,p2=binary,p5=F(1,2))
        seven(5,p3=T(2,1,1),p5=F(1,3))
        seven(5,p2=F(1,0),p3=T(2,1,2),p5=F(1,3))
        seven(5,p5=T(2,4,1))
        seven(5,p2=F(1,0),p5=T(2,4,2))
        seven(5,p3=T(1,0,1),p5=T(2,4,3))
        seven(5,p2=F(1,0),p3=T(1,0,2),p5=T(2,4,3))
        seven(5,p2=F(2,0),p5=T(2,4,4))
        seven(6,p2=T(5,8),p5=F(1,2))
        seven(6,p2=F(3,4),p3=T(2,1,1),p5=F(1,2))
        seven(6,p2=T(4,4),p3=T(2,0,2),p5=F(1,2))
        seven(6,p2=F(3,4),p5=F(2,7))
        seven(6,p2=T(4,4),p5=F(2,12))
        for residue,binary in [(17,F(2,0)),(22,F(3,4)),(2,T(4,4))]:
            seven(6,p2=binary,p3=T(2,1,2),p5=F(2,residue))
            seven(6,p2=binary,p3=F(1,2),p5=F(2,residue))
        seven(6,p2=F(1,0),p3=F(1,1),p5=F(1,3))
        seven(6,p2=F(2,0),p5=F(1,4))
        for digit,binary in [(4,F(2,0)),(5,F(3,4)),(6,T(4,4))]:
            seven(digit,p2=binary,p3=T(1,0,1),p5=T(3,20,4))
        seven(3,p2=F(2,1),p3=F(2,6))
        seven(5,p2=T(3,2),p5=T(3,4,4))
    rows=[]
    seen=set()
    for family in families:
        choices=[]
        for p,E in zip(ps,ds):
            if p not in family:
                choices.append([(0,0)])
                continue
            kind,e,a,d=family[p]
            if kind=='fixed':
                assert e<E
                choices.append([(e,a)])
            else:
                assert e<=E and 0<=a<p**(e-1) and 1<=d<p
                choices.append([(k,a+d*p**(k-1)) for k in range(e,E)]
                               +[(E,a*(p-1)+d-1)])
        for row in product(*choices):
            es=tuple(e for e,r in row)
            if prod(p**e for p,e in zip(ps,es))<minimum:
                continue
            assert es not in seen,('repeated resource',es)
            seen.add(es)
            rows.append(dict(exponents=list(es),residues=[r for e,r in row]))
    return dict(minimum=minimum,primes=ps,depths=ds,rows=rows,
                status='EXPLICIT_PARTIAL_SEED',families=len(families),
                provenance='Owens-inspired explicit branch placements; all rows and resource distinctions are checked independently. No claim that this reproduces or verifies the full thesis.',
                include_seven=include_seven)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--include-seven',action='store_true')
    args=parser.parse_args()
    data=seed(include_seven=args.include_seven)
    name='owens-seven' if args.include_seven else 'owens-inspired'
    Path(f'data/construction-50-{name}-seed.json').write_text(json.dumps(data,indent=2)+'\n')
    cap=prod(p**E for p,E in zip(data['primes'],data['depths'])) if args.include_seven else 200000
    result=build(50,data['primes'],data['depths'],cap=cap,seed_rows=data['rows'])
    result['seed_file']=f'data/construction-50-{name}-seed.json'
    Path(f'data/construction-50-{name}-core.json').write_text(json.dumps(result,indent=2)+'\n')
    print(dict(seed_rows=len(data['rows']),status=result['status'],uncovered=result['exact_uncovered_states']))
