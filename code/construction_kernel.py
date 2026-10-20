"""Exact finite kernels for distinct covering constructions.

Standard-library certificate checker. A coordinate (p,E) has
p**(E-1)*(p-1) states (prefix, first nonzero digit of the remaining tail).
Resource e<E fixes a residue modulo p**e; e=E fixes a tail state and
uses every exponent >=E. See src/construction-search.tex for the finite
completion proof. Different exponent resources must occur at most once.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import prod, isqrt
import argparse
import json


def prime(p):
    return type(p) is int and p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def sizes(primes, depths):
    return [p ** (E-1) * (p-1) for p, E in zip(primes, depths)]


def labels(p, E, e):
    return p**e if e < E else p**(E-1)*(p-1)


def project(x, p, E, e):
    return x if e == E else (x // (p-1)) % (p**e)


def resources(primes, depths, minimum, maximum=None):
    # Traverse only products under the cap. A Cartesian exponent product
    # becomes prohibitive once the search uses primes through 59 or beyond.
    suffix = [1] * (len(primes) + 1)
    for i in range(len(primes) - 1, -1, -1):
        suffix[i] = suffix[i + 1] * primes[i]**depths[i]
    cap = maximum if maximum is not None else suffix[0]
    def visit(i, n, es):
        if n > cap or n * suffix[i] < minimum:
            return
        if i == len(primes):
            yield tuple(es)
            return
        p = primes[i]
        for e in range(depths[i] + 1):
            if n > cap:
                break
            es.append(e)
            yield from visit(i + 1, n, es)
            es.pop()
            n *= p
    yield from visit(0, 1, [])


def capacity(primes, depths, exponents):
    return sum((Fraction(1, prod(labels(p, E, e)
                for p, E, e in zip(primes, depths, es)))
                for es in exponents), Fraction())


def matches(point, row, primes, depths):
    return all(project(x, p, E, e) == r
               for x, p, E, e, r in zip(point, primes, depths,
                                       row['exponents'], row['residues']))


def uncovered_progression(point, primes, depths):
    """CRT lift with the first nonzero tail digit at its earliest level.

    A point missing every row lifts to a progression missing every completed
    row too. Closing classes require x=a mod p**E, unlike this lift.
    """
    residue, modulus = 0, 1
    for state,p,E in zip(point,primes,depths):
        a,d0 = divmod(state,p-1)
        n = p**E
        r = a+(d0+1)*p**(E-1)
        residue += modulus * ((r-residue)*pow(modulus,-1,n) % n)
        modulus *= n
    return {'residue':residue, 'modulus':modulus}


def coordinate_successors(rows, active, i, p, E):
    """Exact sparse partition of one coordinate by a little-endian trie.

    Fixed conditions terminate at a prefix node; tails terminate after the
    final nonzero digit. Unmentioned child branches all have the same active
    rows, so one representative suffices even for enormous state spaces.
    """
    root = [[], {}]
    for j in active:
        e, r = rows[j]['exponents'][i], rows[j]['residues'][i]
        if e == E:
            a, d0 = divmod(r, p-1)
            path = []
            for _ in range(E-1):
                path.append(a % p)
                a //= p
            path.append(d0)
        else:
            path = []
            for _ in range(e):
                path.append(r % p)
                r //= p
        node = root
        for digit in path:
            node = node[1].setdefault(digit, [[], {}])
        node[0].append(j)
    successors = {}

    def walk(node, k, prefix, inherited):
        here = inherited + node[0]
        children = node[1]
        if not children:
            successors.setdefault(tuple(sorted(here)), prefix*(p-1))
            return
        alphabet = p if k < E-1 else p-1
        if len(children) < alphabet:
            missing = 0
            while missing in children:
                missing += 1
            state = ((prefix + missing*p**k)*(p-1)
                     if k < E-1 else prefix*(p-1)+missing)
            successors.setdefault(tuple(sorted(here)), state)
        for digit, child in children.items():
            if k < E-1:
                walk(child, k+1, prefix+digit*p**k, here)
            else:
                successors.setdefault(tuple(sorted(here+child[0])), prefix*(p-1)+digit)

    walk(root, 0, 0, [])
    return successors


def verify(certificate, node_limit=2_000_000):
    """Exact recursive union check; no floating point and no solver dependency.

    Returns a concrete uncovered state if one exists; node limit yields
    UNKNOWN, never a positive certificate. Identical successor sets are
    memoized, so enormous finite products need not be enumerated.
    """
    ps, ds = certificate['primes'], certificate['depths']
    rows, minimum = certificate['rows'], certificate['minimum']
    if type(minimum) is not int or minimum < 2:
        raise ValueError('Minimum must be an integer >=2')
    if not ps or len(ps) != len(ds) or len(set(ps)) != len(ps):
        raise ValueError('Invalid prime coordinates')
    if not all(prime(p) for p in ps) or not all(type(E) is int and E >= 1 for E in ds):
        raise ValueError('Invalid coordinate')
    seen = set()
    for row in rows:
        es, rs = row['exponents'], row['residues']
        if len(es) != len(ps) or len(rs) != len(ps):
            raise ValueError('Invalid row dimension')
        if not all(type(e) is int and 0 <= e <= E for e, E in zip(es, ds)):
            raise ValueError('Invalid exponent')
        if tuple(es) in seen:
            raise ValueError('Repeated exponent resource')
        seen.add(tuple(es))
        if prod(p**e for p, e in zip(ps, es)) < minimum:
            raise ValueError('Modulus below minimum')
        if not all(type(r) is int and 0 <= r < labels(p, E, e)
                   for r, p, E, e in zip(rs, ps, ds, es)):
            raise ValueError('Invalid residue label')
    dims = sizes(ps, ds)
    # Small coordinates first; grouping equal projections is exact.
    order = sorted(range(len(ps)), key=lambda i: dims[i])
    nodes = 0

    @lru_cache(None)
    def walk(k, active):
        nonlocal nodes
        nodes += 1
        if nodes > node_limit:
            raise RuntimeError('node limit')
        if not active:
            return (0,) * (len(ps)-k)
        if any(all(rows[j]['exponents'][i] == 0 for i in order[k:]) for j in active):
            return None
        if k == len(ps):
            return None
        i = order[k]
        successors = coordinate_successors(rows, active, i, ps[i], ds[i])
        for nxt, x in sorted(successors.items(), key=lambda t: len(t[0])):
            tail = walk(k+1, nxt)
            if tail is not None:
                return (x,) + tail
        return None

    try:
        bad = walk(0, tuple(range(len(rows))))
    except RuntimeError:
        return {'status': 'UNKNOWN', 'reason': 'verification node limit', 'nodes': nodes}
    if bad is None:
        return {'status': 'VERIFIED', 'minimum_at_least': minimum,
                'kernel_states': prod(dims), 'resources': len(rows), 'nodes': nodes}
    point = [0] * len(ps)
    for i, x in zip(order, bad):
        point[i] = x
    assert not any(matches(point, row, ps, ds) for row in rows)
    return {'status': 'UNCOVERED', 'point': point, 'nodes': nodes,
            'uncovered_progression':uncovered_progression(point,ps,ds)}


def finite_tail(p, E, state, q):
    """One finite completion, for testing small cases (may be enormous).

    A different auxiliary prime q is assigned to every (row,tail coordinate).
    All returned p-exponents are >=E, including the closing classes.
    """
    if not prime(p) or not prime(q) or p == q or E < 1:
        raise ValueError('Invalid completion parameters')
    a, d0 = divmod(state, p-1)
    if not 0 <= a < p**(E-1):
        raise ValueError('Invalid tail state')
    N = E + q - 1
    out = [(a+(d0+1)*p**(k-1), p**k) for k in range(E, N+1)]
    for j in range(q):
        power = p**(N-j)
        # x=a mod power and x=j mod q.
        residue = (a + power * ((j-a)*pow(power, -1, q) % q)) % (power*q)
        out.append((residue, power*q))
    return out


def completion_summary(certificate):
    """A reproducible finite expansion description, not a new coverage check."""
    ps, ds = certificate['primes'], certificate['depths']
    q = max(ps)
    total = 0
    allocation = []
    for j, row in enumerate(certificate['rows']):
        count = 1
        for i, (e,E) in enumerate(zip(row['exponents'],ds)):
            if e == E:
                q += 1
                while not prime(q):
                    q += 1
                allocation.append({'row':j, 'coordinate':i, 'auxiliary_prime':q,
                                   'last_exponent':E+q-1})
                count *= 2*q
        total += count
    return {'class_count':str(total), 'largest_auxiliary_prime':q,
            'allocation':allocation,
            'minimum_modulus':min(prod(p**e for p,e in zip(ps,row['exponents']))
                                  for row in certificate['rows'])}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('certificate')
    parser.add_argument('--node-limit', type=int, default=2_000_000)
    args = parser.parse_args()
    with open(args.certificate, encoding='utf-8') as f:
        result = verify(json.load(f), args.node_limit)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result['status'] == 'VERIFIED' else 1)
