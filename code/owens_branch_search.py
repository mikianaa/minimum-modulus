"""Restricted, reproducible repair search for Owens's modulus-42 branch.

No result here is a covering of Z. A PARTIAL_BRANCH_CERTIFICATE covers only
its stated progression, using fresh largest-prime families. See the proof in
src/construction-search.tex. Counts are sufficient constructions, not upper
bounds, except for the separately labelled rational relaxation.
"""
from fractions import Fraction as F
from itertools import permutations
from math import gcd, prod
from pathlib import Path
import argparse
import json
import random

ROOT = Path(__file__).resolve().parents[1]


def primes_through(n):
    sieve = bytearray(b'\1') * (n + 1)
    sieve[:2] = b'\0\0'
    for p in range(2, int(n**.5) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = b'\0' * ((n-p*p)//p + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def reserved_235(a, b, c):
    """The early 2,3,5 stages reserve every smooth modulus >=42."""
    if c == 0:
        return ((b == 0 and a >= 6) or (a == 0 and b >= 4)
                or (a >= 1 and b >= 1 and 2**a * 3**b >= 42))
    if c == 1:
        return (b == 0 and a >= 4) or (b == 1 and a >= 2) or b >= 2
    if c == 2:
        return (b == 0 and a >= 1) or b >= 1
    return True


def reserved_7(a, b, c, reserve_36=True):
    """Only unambiguous families from 3.4; deliberately under-reserves.

    In particular, omit the final 125-arrow additions and the extra
    125-arrow*8-arrow. These omissions enlarge the putative spare pool.
    """
    if c == 0:
        return (a >= 3 or (b == 1 and a in (1, 2))
                or (b >= 2 and a <= 1)
                or (reserve_36 and a == 2 and b == 2))
    if c == 1:
        return ((b == 0 and a >= 1) or (b == 1 and a <= 4) or b >= 2)
    if c == 2:
        return True
    return (b == 0 and a <= 2) or (b >= 1 and a <= 1)


def spare_category(a, b, c, d):
    if d == 0 or 2**a * 3**b * 5**c * 7**d < 43:
        return None
    if c == 0 and b == 0 and a in (0, 1) and d >= 2:
        return 0
    if c == 0 and b == 1 and a == 0 and d >= 2:
        return 0
    if c == 0 and b == 0 and a == 2 and d >= 2:
        return 1
    if c == 1 and b == 0 and a == 0 and d >= 2:
        return 2
    if c == 0 and a == 2 and b >= 3:
        return 3
    if c == 1 and b == 1 and a >= 5:
        return 4
    if c >= 3 and b == 0 and a >= 3:
        return 5
    if c >= 3 and b >= 1 and a >= 2:
        return 6
    return None


def density_audit():
    # Each term is sum gcd(n,42)/n over the indicated disjoint family.
    terms = [
        ('7^d, 2*7^d, 3*7^d; d>=2', F(1, 2)),
        ('4*7^d; d>=2', F(1, 12)),
        ('5*7^d; d>=2', F(1, 30)),
        ('4*3^b*7^d; b>=3,d>=1', F(7, 72)),
        ('2^a*3*5*7^d; a>=5,d>=1', F(7, 240)),
        ('2^a*5^c*7^d; a>=3,c>=3,d>=1', F(7, 1200)),
        ('2^a*3^b*5^c*7^d; a>=2,b>=1,c>=3,d>=1', F(7, 400)),
    ]
    checks = 0
    partial = [F(0) for _ in terms]
    for a in range(13):
        for b in range(10):
            for c in range(8):
                n0 = 2**a * 3**b * 5**c
                assert n0 < 42 or reserved_235(a, b, c)
                for d in range(1, 8):
                    n = n0 * 7**d
                    category = spare_category(a, b, c, d)
                    assert (category is not None) == (n >= 43 and not reserved_7(a, b, c))
                    if category is not None:
                        partial[category] += F(gcd(n, 42), n)
                    checks += 1
    assert all(s <= t[1] for s, t in zip(partial, terms))
    total = sum((x[1] for x in terms), F(0))
    assert total == F(2759, 3600)
    return {
        'status': 'RESTRICTED_WHOLE_BRANCH_DENSITY_OBSTRUCTION',
        'scope': 'Reserve the indicated symbolic infinite families; cover the whole 10 mod 42 class using only spare 7-smooth moduli. Does not exclude overlap-assisted or rearranged repairs, or release of powers after finite truncation.',
        'terms': [{'family': name, 'conditional_mass': str(v)} for name, v in terms],
        'total': str(total),
        'without_reserving_36_times_7_powers': str(total + F(7, 36)),
        'finite_inventory_cross_checks': checks,
    }


def evaluate(order, fixed, trace=False):
    count = 1
    result = []
    for p in order:
        e = fixed.get(p, 0)
        nxt = (e + 1) * count + count // (p - 1)
        if trace:
            result.append({'prime': p, 'fixed_exponent': e, 'before': count,
                           'tail_groups': count // (p - 1), 'after': nxt})
        count = nxt
    return result if trace else count


def optimize_order(ps, fixed, rng, iterations):
    # Always include the explicit natural construction as a reproducible seed.
    seeds = [ps[:], sorted(ps, reverse=True)]
    special = [p for p in ps if fixed.get(p, 0)]
    for head in permutations(special):
        seeds.append(list(head) + [p for p in ps if p not in head])
    best = max(seeds, key=lambda x: evaluate(x, fixed))[:]
    score = evaluate(best, fixed)
    current = best[:]
    cur_score = score
    for k in range(iterations):
        if k % 500 == 0:
            current = best[:]
            cur_score = score
        proposal = current[:]
        i, j = rng.sample(range(len(ps)), 2)
        if k % 2:
            proposal[i], proposal[j] = proposal[j], proposal[i]
        else:
            v = proposal.pop(i)
            proposal.insert(j, v)
        v = evaluate(proposal, fixed)
        if v >= cur_score or rng.random() < .015:
            current, cur_score = proposal, v
        if v > score:
            best, score = proposal[:], v
    return best, score


def relaxed_capacity(ps, fixed):
    return prod(F(fixed.get(p, 0) + 1) + F(1, p - 1) for p in ps)


def verify_partial(cert):
    """Check an inductively specified construction, without expanding its DAG.

    At each stage keep old sets, add e fixed-power copies of every set, and
    group disjoint old sets in batches of p-1 for the p^(e+1) arrow copies.
    Every stage uses a previously absent prime. This is what makes the
    resource-disjointness and coverage induction valid.
    """
    q = cert['fresh_prime']
    assert q in primes_through(q) and q >= 97
    order = cert['order']
    fixed = {int(p): e for p, e in cert['fixed_exponents'].items()}
    assert len(set(order)) == len(order)
    assert all(p in primes_through(q - 1) for p in order)
    assert all(isinstance(e, int) and e >= 0 for e in fixed.values())
    assert set(fixed) <= set(order)
    a, m = cert['progression']['residue'], cert['progression']['modulus']
    assert m == prod(p**e for p, e in fixed.items())
    assert 0 <= a < m and m % 42 == 0 and a % 42 == 10
    trace = evaluate(order, fixed, True)
    assert trace == cert['trace']
    assert trace[-1]['after'] >= q - 1
    assert cert['root_children'] == list(range(q - 1))
    assert cert['relative_density_in_42_branch'] == str(F(42, m))
    return {'status': 'VERIFIED_PARTIAL_BRANCH_SCHEMA', 'minimum_modulus': q,
            'progression': cert['progression'],
            'available_disjoint_sets': trace[-1]['after'], 'root_inputs': q - 1,
            'is_cover_of_all_integers': False}


def make_partial(order, fixed, residue, q):
    m = prod(p**e for p, e in fixed.items())
    cert = {
        'status': 'PARTIAL_BRANCH_CERTIFICATE',
        'fresh_prime': q, 'order': order,
        'fixed_exponents': {str(p): e for p, e in fixed.items()},
        'progression': {'residue': residue % m, 'modulus': m},
        'relative_density_in_42_branch': str(F(42, m)),
        'trace': evaluate(order, fixed, True), 'root_children': list(range(q - 1)),
        'recipe': 'Start with the unit set. For each prime p with fixed exponent e, retain the old list, append the p^j copies of each old set (j=1..e, residue matching the progression), then append arrow groups of consecutive old sets in batches of p-1, starting at exponent e+1 with the progression prefix. Discard an incomplete final batch. At the root use q-arrow on the first q-1 final sets. Finitize with auxiliary primes distinct from all construction primes and from each other.',
    }
    cert['verification'] = verify_partial(cert)
    return cert


def partial_union_miss(partials):
    # This misses the new partial covers, not necessarily the retained cover.
    # Every certificate starts with primes 2,5: child 3 is the fixed 8
    # class, and child 5 is the fixed 5 class (indices start at zero).
    a, m = 724, 840
    constraints = [{'residue': a, 'modulus': m}]
    for cert in partials:
        assert cert['order'][:2] == [2, 5]
        target = cert['progression']['residue']
        if a % 8 != target % 8:
            digit, factor = 4, 8
        else:
            assert a % 5 != target % 5
            digit, factor = 6, 5
        q = cert['fresh_prime']
        constraints.append({'residue': digit, 'modulus': q,
                            'missed_child_fixed_factor': factor})
        a += m * ((digit - a) * pow(m, -1, q) % q)
        m *= q
    assert a % 42 == 10
    assert all(a % v['modulus'] == v['residue'] for v in constraints)
    return {'residue': a, 'modulus': m, 'constraints': constraints,
            'scope': 'Uncovered by the six new partial-cover schemas, including their stated finite closures. Not certified uncovered by Owens retained classes.'}


def run(iterations, seed):
    ps = primes_through(89)
    fresh = [p for p in primes_through(200000) if p > 89]
    rng = random.Random(seed)
    results = []
    # Uniform splits of H. First isolate powers of 2, optionally also 3,5,7.
    for r in range(1, 7):
        for s in range(1, 3):
            for t in range(3):
                for u in range(1, 3):
                    fixed = {2: r, 3: s, 5: t, 7: u}
                    branches = 2**(r - 1) * 3**(s - 1) * 5**t * 7**(u - 1)
                    order, count = optimize_order(ps, fixed, rng, iterations)
                    bound = relaxed_capacity(ps, fixed)
                    need = fresh[branches - 1] - 1 if branches <= len(fresh) else None
                    results.append({
                        'fixed_exponents': {str(p): e for p, e in fixed.items()},
                        'equal_branches': branches, 'best_count_found': count,
                        'order': order, 'rational_order_independent_upper_bound': str(bound),
                        'largest_root_prime_needed': need + 1 if need is not None else None,
                        'root_inputs_needed': need,
                        'all_branches_fit_found_count': need is not None and count >= need,
                        'uniform_frozen_base_family_excluded': need is not None and bound < need,
                    })
    fixed = {2: 3, 3: 1, 5: 1, 7: 1}
    # Saved witness found by 12 seeds of 200,000 permutation mutations.
    # Restricting it to primes <q gives counts 109,110,...,114.
    witness_order = [2, 5, 11, 3, 7, 61, 67, 17, 59, 19, 37, 79, 73,
                     41, 71, 83, 43, 29, 23, 31, 47, 13, 103, 53, 97,
                     109, 89, 107, 101]
    qs = [97, 101, 103, 107, 109, 113]
    residues = [a for a in range(10, 840, 42)]
    partials = [make_partial([p for p in witness_order if p < q], fixed, a, q)
                for a, q in zip(residues, qs)]
    bound_127 = relaxed_capacity(primes_through(113), fixed)
    assert bound_127 < 126
    return {
        'status': 'NO_FULL_REPAIR_FOUND', 'minimum_modulus_target': 43,
        'seed': seed, 'iterations_per_configuration': iterations,
        'density_audit': density_audit(), 'uniform_split_experiments': results,
        'partial_replacements': partials,
        'partial_union_relative_density': str(F(len(partials), 20)),
        'partial_union_uncovered_progression': partial_union_miss(partials),
        'partial_fixed_shape_larger_prime_barrier': {
            'first_excluded_root_prime': 127,
            'rational_upper_bound_using_all_smaller_primes': str(bound_127),
            'required_inputs': 126,
            'scope': 'This one-stage-per-prime counting rule on a single class mod 840; all larger root primes also excluded by induction B(next)=B(q)*q/(q-1). Not an obstruction to other constructions.',
        },
        'limitations': [
            'Searches use a restricted one-stage-per-prime set-building rule.',
            'The uniform tests freeze base primes at 89 and do not reuse later fresh primes as bases.',
            'The partial union is not a full replacement for 42.',
            'Original later-stage overlaps and rearrangements are not encoded.',
            'Owens full construction has not been independently verified.',
        ],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--iterations', type=int, default=2500)
    parser.add_argument('--seed', type=int, default=43)
    parser.add_argument('--output', default='data/owens-42-branch-search.json')
    args = parser.parse_args()
    result = run(args.iterations, args.seed)
    (ROOT / args.output).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({
        'status': result['status'],
        'spare_7_smooth_mass': result['density_audit']['total'],
        'configurations': len(result['uniform_split_experiments']),
        'uniform_relaxation_exclusions': sum(v['uniform_frozen_base_family_excluded'] for v in result['uniform_split_experiments']),
        'partial_union_relative_density': result['partial_union_relative_density'],
        'partial_certificates': [v['verification'] for v in result['partial_replacements']],
    }, indent=2))
