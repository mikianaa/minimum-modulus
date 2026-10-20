"""Literal small-family checks independent of the count-only recurrence."""
from copy import deepcopy
from itertools import combinations, product
import json
from pathlib import Path

from owens_branch_search import density_audit, evaluate, verify_partial, partial_union_miss


def literal_step(families, p, e):
    # A row maps primes to (least exponent, is_tail, digit). Fixed factors
    # automatically match F; a tail selects one of its p-1 nonzero digits.
    out = list(families)
    for j in range(1, e + 1):
        out.extend([[dict(row, **{str(p): (j, False, None)}) for row in fam]
                    for fam in families])
    for start in range(0, len(families) - p + 2, p - 1):
        out.append([dict(row, **{str(p): (e + 1, True, digit)})
                    for digit, fam in enumerate(families[start:start+p-1], 1)
                    for row in fam])
    return out


def intervals_intersect(a, b):
    lo_a, tail_a, _ = a
    lo_b, tail_b, _ = b
    if tail_a and tail_b:
        return True
    if tail_a:
        return lo_b >= lo_a
    if tail_b:
        return lo_a >= lo_b
    return lo_a == lo_b


def check_literal(families, ps):
    # No repeated modulus anywhere, even at arbitrarily high tail powers.
    rows = [row for fam in families for row in fam]
    for a, b in combinations(rows, 2):
        assert not all(intervals_intersect(a.get(str(p), (0, False, None)),
                                          b.get(str(p), (0, False, None))) for p in ps)
    # Check every reduced digit state for every separately constructed set.
    for digits in product(*(range(1, p) for p in ps)):
        state = dict(zip(map(str, ps), digits))
        for fam in families:
            assert any(all(not tail or state[p] == digit
                           for p, (_, tail, digit) in row.items()) for row in fam)


def main():
    cases = 0
    for order in [(2, 3, 5, 7), (3, 2, 7, 5), (5, 2, 3, 7)]:
        for es in [(1, 1, 0, 0), (2, 1, 1, 0), (1, 2, 1, 1), (3, 1, 1, 1)]:
            fixed = dict(zip((2, 3, 5, 7), es))
            families = [[{}]]
            used = []
            for p in order:
                families = literal_step(families, p, fixed[p])
                used.append(p)
                check_literal(families, used)
                assert len(families) == evaluate(used, fixed)
            cases += 1
    data = json.loads((Path(__file__).resolve().parents[1] /
                       'data/owens-42-branch-search.json').read_text())
    for cert in data['partial_replacements']:
        verify_partial(cert)
    assert partial_union_miss(data['partial_replacements']) == data['partial_union_uncovered_progression']
    bad = deepcopy(data['partial_replacements'][0])
    bad['order'].append(bad['order'][0])
    try:
        verify_partial(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('duplicate prime accepted')
    audit = density_audit()
    assert audit == data['density_audit']
    print(json.dumps({'status': 'PASS', 'literal_small_constructions': cases,
                      'checked_actual_partial_schemas': len(data['partial_replacements']),
                      'inventory_tuples': audit['finite_inventory_cross_checks']}))


if __name__ == '__main__':
    main()
