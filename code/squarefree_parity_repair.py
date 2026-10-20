"""Exact joint relocation of the d and 2d classes, with a frozen period.

After removing a pair, H(a)=g0(a)+g1(a) is the gain of the odd class.
The even class at a different odd residue b gains max(g0(b),g1(b)).
Two distinct residue cylinders modulo d are disjoint. The optimum distinct
pair is among the top two H entries and the top two even-gain entries.
Only strict improvements are accepted. Every integer in a period is used.
"""
import argparse
from copy import deepcopy
from fractions import Fraction
from itertools import combinations,product
import json
from math import prod
from pathlib import Path
from collections import Counter
from heapq import nlargest

from construction_kernel import verify,prime
from owens_branch_search import primes_through

ROOT = Path(__file__).resolve().parents[1]


def anchor_audit(ps):
    """Universal union bound after prime/twice-prime normalization."""
    assert ps[0] == 3 and all(p > 2 for p in ps)
    supports = [S for k in range(2,len(ps)+1) for S in combinations(range(len(ps)),k)]
    records = []
    for tail in product([0,1],repeat=len(ps)-1):
        colors = (0,)+tail
        dims = [[p-1-int(c == b) for p,c in zip(ps,colors)] for b in [0,1]]
        sizes = [prod(v) for v in dims]
        capacity = 0
        for S in supports:
            volumes = [sizes[b]//prod(dims[b][i] for i in S) for b in [0,1]]
            capacity += sum(volumes)+max(volumes)
        records.append({'anchor_parities':list(colors),'residual_before_composites':sum(sizes),
                        'composite_capacity_upper_bound':capacity,
                        'uncovered_lower_bound':max(0,sum(sizes)-capacity)})
    lower = min(r['uncovered_lower_bound'] for r in records)
    return {'odd_primes':ps,'period':2*prod(ps),'parity_patterns':records,
            'universal_uncovered_lower_bound':lower,
            'status':'COVER_EXCLUDED' if lower else 'NOT_EXCLUDED_BY_THIS_BOUND'}


def top_two(values):
    return nlargest(2,range(len(values)),key=lambda a:(values[a],-a))


def extend(data,q):
    """Add a fresh prime layer, then permit all old classes to move in repair."""
    assert prime(q) and q > max(data['primes'])
    ps = data['primes']
    L = data['period']*q
    classes = deepcopy(data['congruences'])
    counts = [0]*L
    for c in classes:
        for x in range(c['residue'],L,c['modulus']):
            counts[x] += 1
    before = counts.count(0)
    divisors = sorted(prod(p for p,e in zip(ps,es) if e) for es in product([0,1],repeat=len(ps)))
    for d in divisors:
        n = q*d
        if d == 1:
            a = 0
        elif d == 2:
            a = max([1,q+1],key=lambda a:sum(counts[x] == 0 for x in range(a,L,n)))
        else:
            hist = Counter(x % n for x,c in enumerate(counts) if c == 0)
            a = min(hist,key=lambda a:(-hist[a],a)) if hist else 0
        classes.append({'residue':a,'modulus':n})
        for x in range(a,L,n):
            counts[x] += 1
    return {'minimum':3,'primes':ps+[q],'period':L,'congruences':classes,
            'extension_initialization':{'new_prime':q,'uncovered_before':before,
                                        'uncovered_after':counts.count(0)}}


def repair(data,max_sweeps=30):
    ps,L = data['primes'],data['period']
    assert ps[0] == 2 and len(ps) == len(set(ps)) and all(prime(p) for p in ps)
    assert L == prod(ps)
    classes = deepcopy(data['congruences'])
    by_modulus = {c['modulus']:c for c in classes}
    assert len(by_modulus) == len(classes) and data['minimum'] == 3
    counts = [0]*L
    def change(c,delta):
        for x in range(c['residue'],L,c['modulus']):
            counts[x] += delta
    for c in classes:
        change(c,1)
    initial = counts.count(0)
    trace = []
    odd_composites = sorted(n for n in by_modulus if n % 2 and n not in ps)
    for sweep in range(max_sweeps):
        start_holes = counts.count(0)
        moves = 0
        for p in ps[1:]:
            c = by_modulus[2*p]
            change(c,-1)
            gains = [sum(counts[x] == 0 for x in range(1 if b else p+1,L,2*p)) for b in [0,1]]
            old_parity = c['residue'] % 2
            b = max(range(2),key=lambda b:gains[b])
            if gains[b] > gains[old_parity]:
                c['residue'] = 1 if b else p+1
                moves += 1
            change(c,1)
        for d in odd_composites:
            odd,even = by_modulus[d],by_modulus[2*d]
            before = counts.count(0)
            change(odd,-1)
            change(even,-1)
            g = [[0]*d for _ in range(2)]
            for x,c in enumerate(counts):
                if c == 0:
                    g[x % 2][x % d] += 1
            H = [a+b for a,b in zip(*g)]
            G = [max(a,b) for a,b in zip(*g)]
            candidates = [(H[a]+G[b],a,b) for a in top_two(H) for b in top_two(G) if a != b]
            gain,a,b = max(candidates,key=lambda t:(t[0],-t[1],-t[2]))
            removed_holes = sum(H)
            if gain > removed_holes-before:
                odd['residue'] = a
                parity = 0 if g[0][b] >= g[1][b] else 1
                even['residue'] = b if b % 2 == parity else b+d
                moves += 1
            change(odd,1)
            change(even,1)
            assert counts.count(0) <= before
        remaining = counts.count(0)
        trace.append({'sweep':sweep+1,'uncovered_before':start_holes,
                      'uncovered_after':remaining,'strict_moves':moves})
        print(json.dumps(trace[-1]),flush=True)
        if not moves or not remaining:
            break
    result = {'minimum':3,'primes':ps,'depths':[2]*len(ps),'period':L,
              'method':'exact_joint_odd_even_resource_relocation',
              'initial_uncovered':initial,'trace':trace,'congruences':classes,'rows':[]}
    for c in classes:
        n,a = c['modulus'],c['residue']
        es = [int(n % p == 0) for p in ps]
        assert prod(p for p,e in zip(ps,es) if e) == n
        result['rows'].append({'exponents':es,'residues':[a % p if e else 0 for p,e in zip(ps,es)]})
    missed = [x for x in range(L) if all(x % c['modulus'] != c['residue'] for c in classes)]
    assert len(missed) == counts.count(0)
    result['direct_period_check'] = {'status':'UNCOVERED' if missed else 'VERIFIED',
        'period_checked':L,'uncovered_count':len(missed),
        'uncovered_density':str(Fraction(len(missed),L)),
        'first_uncovered_residues':missed[:100],'all_moduli_distinct_squarefree':True}
    result['verification'] = verify(result,node_limit=1000000)
    result['status'] = result['direct_period_check']['status']
    assert result['verification']['status'] == result['status']
    if missed:
        result['uncovered_progression'] = {'residue':missed[0],'modulus':L}
        # For any new q, the possible moduli q*d have d|L. Their maximum
        # total incidence on the residual must be at least q*|R| to repair it.
        divisors = sorted(prod(p for p,e in zip(ps,es) if e) for es in product([0,1],repeat=len(ps)))
        terms = [{'divisor':d,'maximum_residual_intersection':max(Counter(x % d for x in missed).values())}
                 for d in divisors]
        bound = sum(t['maximum_residual_intersection'] for t in terms)
        result['one_fresh_prime_extension_bound'] = {
            'incidence_capacity':bound,'residual_size':len(missed),
            'necessary_q_upper_bound':str(Fraction(bound,len(missed))),
            'divisor_terms':terms,
            'scope':'frozen old classes; only new moduli q*d with d dividing the old period'}
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidate',nargs='?')
    parser.add_argument('--anchor-audit',action='store_true')
    parser.add_argument('--sweeps',type=int,default=30)
    parser.add_argument('--extend-prime',type=int)
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    if args.anchor_audit:
        result = {'method':'exact_normalized_anchor_capacity',
                  'audits':[anchor_audit([p for p in primes_through(y) if p > 2]) for y in [3,5,7,11,13]]}
    else:
        data = json.loads((ROOT/args.candidate).read_text())
        if args.extend_prime:
            data = extend(data,args.extend_prime)
        result = repair(data,args.sweeps)
        if 'extension_initialization' in data:
            result['extension_initialization'] = data['extension_initialization']
        result['seed_file'] = args.candidate
    (ROOT/args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','congruences','audits','one_fresh_prime_extension_bound')},indent=2))
