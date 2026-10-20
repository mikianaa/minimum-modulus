"""Exact squarefree construction experiment: ordinary p-state coordinates.

Each prime occurs at most once in every modulus. The existing independent
kernel verifier is used with E=2 and exponents only 0 or 1, so NO tail row
or finite completion is present. Successful certificates are literal finite
congruence lists. Failed ones receive a directly checked CRT progression.
"""
import argparse
from fractions import Fraction as F
import json
from math import prod
from pathlib import Path
import time

from construction_kernel import resources, verify
from owens_branch_search import primes_through
from structural_cover_search import Diagram, DiagramLimit, ONE, ZERO
from structural_prime_patches import crt_prefix

ROOT = Path(__file__).resolve().parents[1]


def profile(minimum, largest):
    out = []
    for p in primes_through(largest):
        ps = primes_through(p)
        small = [1]
        for q in ps:
            small += [q*n for n in small if q*n < minimum]
        cap = prod(F(p+1,p) for p in ps)-sum((F(1,n) for n in small),F())
        out.append({'largest_prime':p,'full_squarefree_capacity':str(cap),
                    'excluded_by_density':cap < 1})
    return out


def search(minimum,largest,seconds,node_limit):
    start = time.monotonic()
    ps = primes_through(largest)
    pool = list(resources(ps,[1]*len(ps),minimum))
    pool.sort(key=lambda es:prod(p**e for p,e in zip(ps,es)))
    cap = sum((F(1,prod(p**e for p,e in zip(ps,es))) for es in pool),F())
    result = {'minimum':minimum,'primes':ps,'depths':[2]*len(ps),
              'model':'ordinary_residues_squarefree_no_tails',
              'resource_count':len(pool),'full_squarefree_capacity':str(cap),
              'period':prod(ps),'node_limit':node_limit,'seconds_budget':seconds,
              'rows':[],'congruences':[],'capacity_profile':profile(minimum,largest)}
    if cap < 1:
        result['status'] = 'EXCLUDED_BY_EXACT_CAPACITY'
        return result
    diagram = Diagram(ps,node_limit,start+seconds)
    root = ONE
    before = prod(ps)
    outcome = 'RESOURCES_EXHAUSTED'
    processed = 0
    try:
        for es in pool:
            if time.monotonic()-start >= seconds:
                raise DiagramLimit('time limit')
            selected = [i for i,e in enumerate(es) if e]
            gain,assignment,_ = diagram.maximum_cofactor(root,selected)
            if gain:
                new_root = diagram.erase(root,assignment)
                after = diagram.count(new_root)
                assert before-after == gain
                rs = [assignment.get(i,0) for i in range(len(ps))]
                a,n = crt_prefix(es,rs,ps)
                assert n >= minimum and all(e in (0,1) for e in es)
                result['rows'].append({'exponents':list(es),'residues':rs})
                result['congruences'].append({'residue':a,'modulus':n})
                root,before = new_root,after
            processed += 1
            if root == ZERO:
                outcome = 'EXACT_COVER_FOUND'
                break
            if len(diagram.nodes) > node_limit//4:
                old = diagram
                diagram,root = old.compact(root)
                old.count.cache_clear()
    except DiagramLimit as exc:
        outcome = str(exc).upper().replace(' ','_')
    result['search_outcome'] = outcome
    result['resources_processed'] = processed
    result['uncovered_residues_per_period'] = before
    result['exact_uncovered_density'] = str(F(before,prod(ps)))
    result['verification'] = verify(result,node_limit=300000)
    result['status'] = result['verification']['status']
    assert len({c['modulus'] for c in result['congruences']}) == len(result['congruences'])
    assert all(prod(ps) % c['modulus'] == 0 for c in result['congruences'])
    if result['status'] == 'UNCOVERED':
        L = prod(ps)
        R = result['verification']['uncovered_progression']['residue'] % L
        for c in result['congruences']:
            assert R % c['modulus'] != c['residue']
        result['direct_modular_check'] = {'status':'UNCOVERED_PROGRESSION_VERIFIED',
            'rows_checked':len(result['rows']),'residue':R,'modulus':L,
            'all_moduli_distinct_squarefree':True}
    if root == ZERO:
        assert result['status'] == 'VERIFIED'
    result['elapsed_seconds'] = round(time.monotonic()-start,3)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum',type=int,default=30)
    parser.add_argument('--largest-prime',type=int,default=37)
    parser.add_argument('--seconds',type=int,default=60)
    parser.add_argument('--nodes',type=int,default=100000)
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    out = search(args.minimum,args.largest_prime,args.seconds,args.nodes)
    (ROOT/args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('rows','congruences','capacity_profile')},indent=2))
