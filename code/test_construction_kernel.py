"""Independent finite checks of the construction kernel and tail completion."""
from itertools import product
from math import lcm, prod, gcd
from pathlib import Path
import json
import random
from construction_kernel import (finite_tail, labels, matches, project,
                                 resources, sizes, verify, coordinate_successors)


def main():
    rng = random.Random(271828)
    partition_cases = 0
    for p in [2,3,5,7]:
        for E in [1,2,3,4]:
            for _ in range(10):
                rows = []
                for j in range(20):
                    e = rng.randrange(E+1)
                    rows.append({'exponents':[e],
                                 'residues':[rng.randrange(labels(p,E,e))]})
                active = tuple(range(len(rows)))
                expected = {tuple(j for j in active if project(x,p,E,rows[j]['exponents'][0])
                                  == rows[j]['residues'][0]) for x in range(p**(E-1)*(p-1))}
                got = coordinate_successors(rows,active,0,p,E)
                assert set(got) == expected
                for rows_at_state,x in got.items():
                    assert rows_at_state == tuple(j for j in active
                        if project(x,p,E,rows[j]['exponents'][0]) == rows[j]['residues'][0])
                partition_cases += 1
    # This coordinate has 2^49 states; sparse partition must not enumerate it.
    huge = verify({'minimum':2,'primes':[2],'depths':[50],
                   'rows':[{'exponents':[1],'residues':[0]}]})
    assert huge['status'] == 'UNCOVERED'
    # Compare capped recursion to independent full Cartesian enumeration.
    for ps, ds in [([2,3,5], [3,2,2]), ([7,3,2], [2,3,1])]:
        for minimum in [2,43,1000]:
            for maximum in [None,1,42,1000,10000]:
                brute_resources = [es for es in product(*(range(d+1) for d in ds))
                                   if prod(p**e for p,e in zip(ps,es)) >= minimum
                                   and (maximum is None or prod(p**e for p,e in zip(ps,es)) <= maximum)]
                assert list(resources(ps,ds,minimum,maximum)) == brute_resources
    cases = 0
    for ps, ds in [([2,3],[2,2]), ([2,3,5],[2,1,1]), ([3,5],[2,1])]:
        dims = sizes(ps, ds)
        exps = list(resources(ps, ds, 2))
        for _ in range(80):
            chosen = rng.sample(exps, rng.randrange(len(exps)+1))
            rows = [{'exponents': list(es), 'residues': [rng.randrange(labels(p,E,e))
                     for p,E,e in zip(ps,ds,es)]} for es in chosen]
            cert = {'primes':ps, 'depths':ds, 'minimum':2, 'rows':rows}
            brute = all(any(matches(x,row,ps,ds) for row in rows)
                        for x in product(*(range(d) for d in dims)))
            checked = verify(cert)
            assert (checked['status'] == 'VERIFIED') == brute
            if not brute:
                lift = checked['uncovered_progression']
                R,L = lift['residue'],lift['modulus']
                # Independent arithmetic: each row has a coordinate whose
                # every finite class is incompatible with the whole lift.
                for row in rows:
                    witnesses = []
                    for p,E,e,r in zip(ps,ds,row['exponents'],row['residues']):
                        classes = finite_tail(p,E,r,7) if e == E else [(r,p**e)]
                        witnesses.append(all((R-a) % gcd(L,n) != 0 for a,n in classes))
                    assert any(witnesses)
            cases += 1
    # Independent residue computation, including fixed prime powers.
    for p in [2,3,5]:
        for E in [1,2,3]:
            for x in range(1, p**(E+2)):
                prefix = x % p**(E-1)
                y = (x-prefix)//p**(E-1)
                if y == 0:
                    d = 1
                else:
                    while y % p == 0:
                        y //= p
                    d = y % p
                state = prefix*(p-1) + d-1
                for e in range(E):
                    assert project(state,p,E,e) == x % p**e
    tail_cases = 0
    for p,q in [(2,3),(3,2),(3,5),(5,2),(5,3)]:
        for E in [1,2,3]:
            states = range(p**(E-1)*(p-1))
            for state in rng.sample(list(states), min(4,len(states))):
                classes = finite_tail(p,E,state,q)
                mods = [n for a,n in classes]
                assert len(set(mods)) == len(mods) and min(mods) >= p**E
                period = lcm(*mods)
                a,d0 = divmod(state,p-1)
                N = E+q-1
                covered = bytearray(period)
                for r,n in classes:
                    assert r % p**(E-1) == a
                    covered[r::n] = b'\x01' * len(range(r,period,n))
                for x in range(a,period,p**(E-1)):
                    if (x-a) % p**N == 0:
                        assert covered[x]
                    else:
                        y = (x-a)//p**(E-1)
                        while y % p == 0:
                            y //= p
                        if y % p == d0+1:
                            assert covered[x]
                tail_cases += 1
    root = Path(__file__).resolve().parents[1]
    controls = {}
    for m in [3,6]:
        cert = json.loads((root/f'data/construction-kernel-control-{m}.json').read_text())
        controls[str(m)] = verify(cert)
        assert controls[str(m)]['status'] == 'VERIFIED'
        bad = dict(cert, rows=cert['rows']+[cert['rows'][0]])
        try:
            verify(bad)
        except ValueError:
            pass
        else:
            raise AssertionError('Duplicate resource accepted')
    summary = {'random_kernel_comparisons':cases, 'exhaustive_finite_tail_checks':tail_cases,
               'coordinate_partition_comparisons':partition_cases,
               'capped_resource_enumeration_comparisons':30,
               'huge_coordinate_states_checked_symbolically':2**49,
               'controls':controls, 'status':'PASS'}
    (root/'data/construction-kernel-tests.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
