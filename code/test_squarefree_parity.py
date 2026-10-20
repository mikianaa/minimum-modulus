"""Independent finite checks for normalization, pair moves, and anchor bounds."""
from collections import Counter
from itertools import product,combinations
import json
from math import prod
from pathlib import Path
import random

from squarefree_parity_repair import top_two,anchor_audit

ROOT = Path(__file__).resolve().parents[1]


def covered(rows,point):
    return any(all(point[i] == a for i,a in row.items()) for row in rows.values())


def main():
    rng = random.Random(303)
    pair_checks = 0
    for d in [3,5,7,11,15]:
        for _ in range(100):
            g = [[rng.randrange(10) for _ in range(d)] for _ in [0,1]]
            H = [a+b for a,b in zip(*g)]
            G = [max(a,b) for a,b in zip(*g)]
            reduced = max(H[a]+G[b] for a in top_two(H) for b in top_two(G) if a != b)
            literal = max(H[a]+(g[e][b] if a != b else 0)
                          for a,b,e in product(range(d),range(d),range(2)))
            assert reduced == literal
            pair_checks += 1
    ps = [2,3,5,7]
    points = list(product(*(range(p) for p in ps)))
    supports = [S for k in range(1,5) for S in combinations(range(4),k) if S != (0,)]
    normalization_checks = 0
    for _ in range(100):
        rows = {S:{i:rng.randrange(ps[i]) for i in S}
                for S in rng.sample(supports,rng.randrange(len(supports)+1))}
        before = sum(covered(rows,x) for x in points)
        for S in supports:
            rows.setdefault(S,{i:0 for i in S})
        for i in range(1,4):
            a,b = rows[(i,)][i],rows[(0,i)][i]
            if a == b:
                b = (a+1) % ps[i]
                rows[(0,i)][i] = b
            order = [a,b]+[c for c in range(ps[i]) if c not in (a,b)]
            permutation = {v:k for k,v in enumerate(order)}
            for row in rows.values():
                if i in row:
                    row[i] = permutation[row[i]]
        flip = rows[(0,1)][0]
        for row in rows.values():
            if 0 in row:
                row[0] ^= flip
            odd_coordinates = [i for i in row if i]
            if len(odd_coordinates) >= 2 and any(row[i] == 0 for i in odd_coordinates):
                for i in odd_coordinates:
                    if row[i] == 0:
                        row[i] = 1
        assert all(rows[(i,)][i] == 0 and rows[(0,i)][i] == 1 for i in range(1,4))
        assert rows[(0,1)][0] == 0
        assert sum(covered(rows,x) for x in points) >= before
        normalization_checks += 1
    bound_checks = 0
    for ps in [[3],[3,5],[3,5,7],[3,5,7,11],[3,5,7,11,13]]:
        audit = anchor_audit(ps)
        supports = [S for k in range(2,len(ps)+1) for S in combinations(range(len(ps)),k)]
        for record in audit['parity_patterns']:
            colors = record['anchor_parities']
            residual = [x for x in range(2*prod(ps))
                        if all(x % p and not (x % p == 1 and x % 2 == b) for p,b in zip(ps,colors))]
            actual = 0
            for S in supports:
                d = prod(ps[i] for i in S)
                actual += max(Counter(x % d for x in residual).values())
                actual += max(Counter(x % (2*d) for x in residual).values())
            assert len(residual) == record['residual_before_composites']
            assert actual == record['composite_capacity_upper_bound']
            bound_checks += 1
    out = {'status':'PASS','joint_pair_exhaustive_comparisons':pair_checks,
           'normalization_direct_period_comparisons':normalization_checks,
           'anchor_capacity_direct_period_checks':bound_checks}
    (ROOT/'data/squarefree-parity-tests.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
