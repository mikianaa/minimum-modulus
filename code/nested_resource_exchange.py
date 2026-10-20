"""Exact coupled relocation of two nested CRT-cylinder resources.

The selected digits of the coarse resource must be a subset of those of
the fine resource. Their placements are optimized jointly over the entire
residual, using two coarse maxima and two fine maxima. No sampling.
Independent kernel verification is the only positive coverage certificate.
"""
import argparse
from copy import deepcopy
from fractions import Fraction
import json
from math import prod
from pathlib import Path
import time

from construction_kernel import verify
from structural_cover_search import (Diagram, DiagramLimit, ONE, coordinates,
    condition_from_row, row_from_condition)
from verify_construction_searches import verify_uncovered_progression


def best_nested_pair(diagram, root, coarse, fine):
    """Return globally maximum union gain and its two assignments.

    A fine cylinder is either inside a coarse cylinder or disjoint from it.
    For disjoint placements the gain is H(a)+G(b). For any fixed b, an
    optimal a is among the top two H values, since only one a is forbidden.
    Nested placements cannot beat the best coarse gain, already included.
    """
    coarse, fine = frozenset(coarse), frozenset(fine)
    if not coarse <= fine:
        raise ValueError('The resource supports must be nested')
    h, a, _ = diagram.maximum_cofactor(root, coarse)
    best = (h, a, {v: a.get(v, 0) for v in fine})
    if not coarse or prod(diagram.domains[v] for v in coarse) == 1:
        return best
    after_first = diagram.erase(root, a)
    h2, a2, _ = diagram.maximum_cofactor(after_first, coarse)
    for score, assignment in [(h, a), (h2, a2)]:
        rest = diagram.erase(root, assignment)
        g, b, _ = diagram.maximum_cofactor(rest, fine)
        # If g=0 the fine row can be nested without affecting the union.
        if score + g > best[0]:
            best = (score + g, assignment, b)
    return best


def rebuild(rows, ps, ds, domains, blocks, layout, limit, deadline, omit=()):
    diagram, root = Diagram(domains, limit, deadline), ONE
    omit = set(omit)
    for j, row in enumerate(rows):
        if j in omit:
            continue
        root = diagram.erase(root, condition_from_row(row, ps, ds, blocks))
        if len(diagram.nodes) > limit // 3:
            old = diagram
            diagram, root = old.compact(root)
            old.count.cache_clear()
    return diagram, root


def exchange(data, seconds=60, limit=200000, sweeps=2, pair_limit=100):
    result = deepcopy(data)
    ps, ds, rows = result['primes'], result['depths'], result['rows']
    layout = result.get('layout', 'large-prime-first')
    domains, blocks = coordinates(ps, ds, layout)
    start = time.monotonic()
    deadline = start + seconds
    diagram, root = rebuild(rows, ps, ds, domains, blocks, layout, limit, deadline)
    before = diagram.count(root)
    total = prod(domains)
    if 'exact_uncovered_states' in data:
        assert before == data['exact_uncovered_states']
    initial = before
    diagram.count.cache_clear()
    supports = [frozenset(v for b,e in zip(blocks,r['exponents']) for v in b[:e])
                for r in rows]
    pairs = []
    for i in range(len(rows)):
        for j in range(i+1,len(rows)):
            a,b = supports[i],supports[j]
            if a <= b:
                pairs.append((i,j))
            elif b <= a:
                pairs.append((j,i))
    # Prioritize large cylinders and small changes of exponent resource.
    pairs.sort(key=lambda ij:(prod(domains[v] for v in supports[ij[0]]),
                             len(supports[ij[1]]-supports[ij[0]]),ij))
    pairs = pairs[:pair_limit]
    trace = []
    attempts = 0
    stop = 'SWEEPS_EXHAUSTED'
    try:
        for sweep in range(sweeps):
            improved = 0
            for i,j in pairs:
                if time.monotonic() >= deadline:
                    raise DiagramLimit('time limit')
                diagram, root = rebuild(rows, ps, ds, domains, blocks, layout,
                                        limit, deadline, (i,j))
                residual = diagram.count(root)
                gain, a, b = best_nested_pair(diagram, root, supports[i], supports[j])
                attempts += 1
                after = residual - gain
                assert after <= before, 'Joint maximum cannot worsen incumbent'
                if after < before:
                    proposed_i = row_from_condition(rows[i]['exponents'],a,ps,ds,blocks)
                    proposed_j = row_from_condition(rows[j]['exponents'],b,ps,ds,blocks)
                    check_root = diagram.erase(diagram.erase(root,a),b)
                    assert diagram.count(check_root) == after
                    rows[i], rows[j] = proposed_i, proposed_j
                    info = dict(sweep=sweep+1, rows=[i,j],
                                before=before, after=after,
                                uncovered_fraction=str(Fraction(after,total)))
                    trace.append(info)
                    print(json.dumps(info),flush=True)
                    before = after
                    improved += 1
                diagram.count.cache_clear()
                if before == 0:
                    stop = 'CANDIDATE_COVER'
                    break
            if before == 0:
                break
            if improved == 0:
                stop = 'LOCAL_OPTIMUM_FOR_TESTED_PAIRS'
                break
    except DiagramLimit as exc:
        stop = str(exc).upper().replace(' ','_')
    finally:
        diagram.count.cache_clear()
    result.update(method='exact_nested_resource_pair_exchange',
                  parent_method=data.get('method'), initial_uncovered_states=initial,
                  exact_uncovered_states=before, kernel_states=total,
                  exact_uncovered_fraction=str(Fraction(before,total)),
                  exchange_trace=trace, tested_pairs=attempts,
                  pair_pool_size=len(pairs), exchange_outcome=stop,
                  exchange_elapsed_seconds=time.monotonic()-start)
    # Recompute the union from scratch, separate from the gain calculation.
    diagram, root = rebuild(rows,ps,ds,domains,blocks,layout,limit,float('inf'))
    assert diagram.count(root) == before
    diagram.count.cache_clear()
    result['verification'] = verify(result)
    result['status'] = result['verification']['status']
    if result['status'] == 'UNCOVERED':
        result['independent_witness_check'] = verify_uncovered_progression(result)
    if before == 0:
        assert result['status'] == 'VERIFIED'
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidate')
    parser.add_argument('--seconds',type=float,default=60)
    parser.add_argument('--nodes',type=int,default=200000)
    parser.add_argument('--sweeps',type=int,default=2)
    parser.add_argument('--pairs',type=int,default=100)
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    data = json.loads(Path(args.candidate).read_text())
    result = exchange(data,args.seconds,args.nodes,args.sweeps,args.pairs)
    result['seed_file'] = args.candidate
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','exchange_outcome','tested_pairs',
          'initial_uncovered_states','exact_uncovered_states','exact_uncovered_fraction']},indent=2))
