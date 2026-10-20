"""Direct-period checks of the squarefree model and its verifier adapter."""
from fractions import Fraction as F
from itertools import product
import json
from math import prod
from pathlib import Path
import random
from construction_kernel import verify
from search_squarefree_construction import profile
from structural_prime_patches import crt_prefix

ROOT = Path(__file__).resolve().parents[1]

def main():
    rng = random.Random(30)
    checks = 0
    for ps in [[2,3,5],[2,3,5,7]]:
        options = [es for es in product([0,1],repeat=len(ps)) if any(es)]
        for _ in range(100):
            chosen = rng.sample(options,rng.randrange(len(options)+1))
            rows,classes = [],[]
            for es in chosen:
                rs = [rng.randrange(p) if e else 0 for p,e in zip(ps,es)]
                rows.append({'exponents':list(es),'residues':rs})
                classes.append(crt_prefix(es,rs,ps))
            cert = {'minimum':2,'primes':ps,'depths':[2]*len(ps),'rows':rows}
            truth = all(any(x % n == a for a,n in classes) for x in range(prod(ps)))
            checked = verify(cert)
            assert (checked['status'] == 'VERIFIED') == truth
            if not truth:
                r = checked['uncovered_progression']['residue'] % prod(ps)
                assert all(r % n != a for a,n in classes)
            checks += 1
        total = sum((F(1,prod(p**e for p,e in zip(ps,es))) for es in options
                     if prod(p**e for p,e in zip(ps,es)) >= 30),F())
        assert total == F(profile(30,max(ps))[-1]['full_squarefree_capacity'])
    c = json.loads((ROOT/'data/construction-squarefree-control-2.json').read_text())
    assert prod(c['primes']) == 210
    ns = [r['modulus'] for r in c['congruences']]
    assert len(ns) == len(set(ns)) and min(ns) == 2
    assert all(210 % n == 0 for n in ns)
    assert all(any(x % r['modulus'] == r['residue'] for r in c['congruences'])
               for x in range(210))
    out = {'status':'PASS','random_direct_period_comparisons':checks,
           'squarefree_control_period_exhaustively_checked':210,
           'squarefree_control_classes':len(ns)}
    (ROOT/'data/construction-squarefree-tests.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__ == '__main__':
    main()
