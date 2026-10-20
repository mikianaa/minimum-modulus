"""Capacity screening and bounded sample-guided construction search.

Sample coverage is diagnostic only. Only the independent exact kernel
checker may certify a covering system. All successful-looking samples are
checked, and failures retain an explicit uncovered arithmetic progression.
"""
import argparse
from fractions import Fraction as F
import json
from math import prod
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tmp'/'construction-deps'))
from construction_kernel import resources, capacity, labels, sizes, verify
from owens_branch_search import primes_through, reserved_235, reserved_7


def depths_for(ps, minimum):
    ds = []
    for p in ps:
        e = 1
        while p**e < minimum:
            e += 1
        ds.append(e)
    return ds


def smooth_capacity(minimum, ps):
    smooth = {1}
    for p in ps:
        for n in list(smooth):
            v = n*p
            while v < minimum:
                smooth.add(v)
                v *= p
    return prod(F(p,p-1) for p in ps) - sum((F(1,n) for n in smooth), F())


def profile(minimum):
    early_moduli = []
    for n in range(42,minimum):
        rest, es = n, []
        for p in [2,3,5,7]:
            e = 0
            while rest % p == 0:
                rest //= p
                e += 1
            es.append(e)
        a,b,c,d = es
        if rest == 1 and ((d == 0 and reserved_235(a,b,c))
                          or (d > 0 and reserved_7(a,b,c))):
            early_moduli.append(n)
    support = []
    for p in primes_through(149):
        c = smooth_capacity(minimum, primes_through(p))
        support.append({'largest_allowed_prime': p, 'capacity': str(c),
                        'decimal_capacity': float(c), 'excluded': c < 1})
    finite = []
    for y in [59,73,97,149,311]:
        ps = primes_through(y)
        ds = depths_for(ps, minimum)
        for cap in [10000,30000,100000]:
            es = list(resources(ps, ds, minimum, cap))
            c = capacity(ps, ds, es)
            finite.append({'largest_prime': y, 'base_cap': cap,
                           'primes': ps, 'depths': ds,
                           'resources': len(es), 'capacity': str(c),
                           'decimal_capacity': float(c), 'excluded': c < 1})
    s, end = F(), minimum-1
    while s < 1:
        end += 1
        s += F(1,end)
    return {'minimum_target': minimum,
            'status': 'NECESSARY_CONDITIONS_ONLY_NOT_A_CONSTRUCTION',
            'owens_audited_early_moduli_below_target': early_moduli,
            'owens_audited_early_classes_to_change': len(early_moduli),
            'smooth_support_bounds': support, 'finite_kernel_profiles': finite,
            'cardinality_lower_bound': end-minimum+1,
            'maximum_modulus_lower_bound': end,
            'harmonic_sum_at_endpoint': str(s),
            'harmonic_sum_before_endpoint': str(s-F(1,end))}


def search(minimum, largest, cap, sample_count, sweeps, seconds, seed, uniform_weights=False):
    import numpy as np
    start = time.monotonic()
    ps = primes_through(largest)
    ds = depths_for(ps, minimum)
    exps = sorted(resources(ps, ds, minimum, cap),
                  key=lambda es: prod(p**e for p,e in zip(ps,es)))
    mass = capacity(ps, ds, exps)
    result = {'minimum': minimum, 'primes': ps, 'depths': ds,
              'maximum_resource_base': cap, 'resource_count': len(exps),
              'capacity': str(mass), 'seed': seed,
              'training_samples': sample_count, 'sweeps_budget': sweeps,
              'seconds_budget': seconds,
              'weight_policy': 'uniform' if uniform_weights else 'uncovered_times_1.3_each_sweep',
              'iterations': [], 'rows': []}
    if mass < 1:
        result['status'] = 'EXCLUDED_BY_EXACT_CAPACITY'
        return result
    if cap > 65535:
        raise ValueError('This compact cache requires cap <=65535')
    rng = np.random.default_rng(seed)
    dims = sizes(ps, ds)
    states = [rng.integers(0, n, size=sample_count, dtype=np.int64) for n in dims]
    projections = {(i,e): (states[i] if e == ds[i] else
                           (states[i]//(p-1)) % p**e)
                   for i,p in enumerate(ps) for e in range(1,ds[i]+1)}
    radices = [[labels(p,E,e) for p,E,e in zip(ps,ds,es)] for es in exps]
    def make_codes(j):
        es = exps[j]
        h = np.zeros(sample_count, dtype=np.int64)
        for i,e in enumerate(es):
            if e:
                h = h*radices[j][i] + projections[i,e]
        return h.astype(np.uint16)
    cached = len(exps)*sample_count*2 <= 600_000_000
    if cached:
        codes = np.empty((len(exps),sample_count), dtype=np.uint16)
        for j in range(len(exps)):
            codes[j] = make_codes(j)
        del projections
    else:
        codes = None
    result['cached_sample_codes'] = cached
    del states
    assignment = np.full(len(exps), -1, dtype=np.int32)
    cover_count = np.zeros(sample_count, dtype=np.int32)
    weights = np.ones(sample_count, dtype=np.float64)
    best_assignment = None
    best_uncovered = sample_count+1
    for sweep in range(sweeps):
        if time.monotonic()-start >= seconds:
            break
        order = range(len(exps)) if sweep == 0 else rng.permutation(len(exps))
        for j in order:
            row_codes = codes[j] if cached else make_codes(j)
            old = assignment[j]
            if old >= 0:
                cover_count[row_codes == old] -= 1
            gain = np.bincount(row_codes, weights=weights*(cover_count == 0),
                               minlength=prod(radices[j]))
            # Retain an existing choice in a tie, to avoid gratuitous churn.
            chosen = int(gain.argmax())
            if old >= 0 and gain[old] == gain[chosen]:
                chosen = int(old)
            assignment[j] = chosen
            cover_count[row_codes == chosen] += 1
        uncovered = int(np.count_nonzero(cover_count == 0))
        assert int(cover_count.min()) >= 0
        if uncovered < best_uncovered:
            best_uncovered = uncovered
            best_assignment = assignment.copy()
        info = {'sweep': sweep+1, 'training_uncovered': uncovered,
                'best_training_uncovered': best_uncovered,
                'elapsed_seconds': round(time.monotonic()-start,3)}
        result['iterations'].append(info)
        print(json.dumps(info), flush=True)
        if uncovered == 0:
            break
        if not uniform_weights:
            weights[cover_count == 0] *= 1.3
            weights = np.minimum(weights, 1000.)
    if best_assignment is None:
        result['status'] = 'UNKNOWN_SETUP_TIME_LIMIT'
        return result
    rows = []
    for es, radix, code in zip(exps, radices, best_assignment):
        rs = [0]*len(ps)
        code = int(code)
        for i in range(len(ps)-1,-1,-1):
            code, rs[i] = divmod(code, radix[i])
        assert code == 0
        rows.append({'exponents':list(es), 'residues':rs})
    result['rows'] = rows
    del codes
    # Fresh diagnostic samples never drive the choice of best candidate.
    validation_count = 32768
    states = [rng.integers(0,n,size=validation_count,dtype=np.int64) for n in dims]
    covered = np.zeros(validation_count, dtype=np.bool_)
    for row in rows:
        matches = np.ones(validation_count, dtype=np.bool_)
        for i,e in enumerate(row['exponents']):
            if e:
                values = states[i] if e == ds[i] else (states[i]//(ps[i]-1)) % ps[i]**e
                matches &= values == row['residues'][i]
        covered |= matches
    result['sample_diagnostics'] = {
        'best_training_uncovered': best_uncovered,
        'validation_samples': validation_count,
        'validation_uncovered': int(np.count_nonzero(~covered)),
        'interpretation': 'Samples are diagnostic; no claim of exact density or covering.'}
    result['verification'] = verify(result, node_limit=300000)
    result['status'] = ('VERIFIED' if result['verification']['status'] == 'VERIFIED'
                        else 'CANDIDATE_REJECTED' if result['verification']['status'] == 'UNCOVERED'
                        else 'UNKNOWN_VERIFICATION_LIMIT')
    result['elapsed_seconds'] = round(time.monotonic()-start,3)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', action='store_true')
    parser.add_argument('--minimum', type=int, default=1000)
    parser.add_argument('--largest-prime', type=int, default=149)
    parser.add_argument('--cap', type=int, default=10000)
    parser.add_argument('--samples', type=int, default=32768)
    parser.add_argument('--sweeps', type=int, default=12)
    parser.add_argument('--seconds', type=int, default=180)
    parser.add_argument('--seed', type=int, default=1000)
    parser.add_argument('--uniform-weights', action='store_true')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = (profile(args.minimum) if args.profile else
              search(args.minimum,args.largest_prime,args.cap,args.samples,
                     args.sweeps,args.seconds,args.seed,args.uniform_weights))
    (ROOT/args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in
                      ('rows','iterations','smooth_support_bounds','finite_kernel_profiles',
                       'harmonic_sum_at_endpoint','harmonic_sum_before_endpoint')},indent=2))
