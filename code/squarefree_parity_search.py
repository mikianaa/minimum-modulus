"""Joint parity construction after a lossless prime/twice-prime normalization.

For minimum 3, add every missing p and 2p resource (p odd), normalize p to
residue 0, move a redundant 2p class off 0, and relabel its nonzero digit 1.
Only the parity of 2p remains free. All non-unit coordinates are covered by
the prime classes. Optimize every remaining unit CRT point simultaneously.
Solver status is separate from independent coverage checks; no global
nonexistence claim or independently checked solver optimality is made.
"""
import argparse
from fractions import Fraction
from itertools import combinations, product
import json
from math import prod
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tmp'/'construction-deps'))
from ortools.sat.python import cp_model
from construction_kernel import verify
from owens_branch_search import primes_through
from structural_prime_patches import crt_prefix


def search(largest=11,minimum=3,seconds=60):
    assert minimum in (2,3)
    start = time.monotonic()
    ps = [p for p in primes_through(largest) if p > 2]
    assert ps
    model = cp_model.CpModel()
    colors = [model.new_bool_var(f'parity_{p}') for p in ps]
    # Swapping the parity labels is a symmetry when modulus 2 is absent.
    if minimum == 3:
        model.add(colors[0] == 0)
    supports = [S for k in range(2,len(ps)+1) for S in combinations(range(len(ps)),k)]
    choices = {}
    for S in supports:
        assignments = list(product(*(range(1,ps[i]) for i in S)))
        for layer in [-1,0,1]:
            # The odd modulus has layer=-1; the even modulus has both layers
            # in a single exactly-one constraint below.
            for a in assignments:
                choices[S,layer,a] = model.new_bool_var('')
        model.add_exactly_one([choices[S,-1,a] for a in assignments])
        model.add_exactly_one([choices[S,b,a] for b in [0,1] for a in assignments])
    holes = []
    for x in product(*(range(1,p) for p in ps)):
        for b in ([0,1] if minimum == 3 else [1]):
            cover = [colors[i] if b else colors[i].Not() for i,t in enumerate(x) if t == 1]
            for S in supports:
                a = tuple(x[i] for i in S)
                cover.extend([choices[S,-1,a],choices[S,b,a]])
            h = model.new_bool_var('')
            model.add_bool_or(cover+[h])
            model.add_bool_and([v.Not() for v in cover]).only_enforce_if(h)
            holes.append(h)
    model.minimize(sum(holes))
    result = {'minimum':minimum,'primes':[2]+ps,'depths':[2]*(len(ps)+1),
              'method':'joint_parity_normalized_exact_period_model',
              'period':2*prod(ps),'unit_points_constrained':len(holes),
              'selection_variables':len(choices),'seconds_solver_budget':seconds,
              'normalization':'p:0; 2p:1 modulo p and free parity; no prime powers',
              'rows':[],'congruences':[]}
    result['model_build_seconds'] = round(time.monotonic()-start,3)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 3
    status = solver.solve(model)
    result['solver_status'] = solver.status_name(status)
    result['solver_best_bound'] = solver.best_objective_bound
    if status not in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        result['status'] = 'NO_CANDIDATE_RETURNED'
        result['elapsed_seconds'] = round(time.monotonic()-start,3)
        return result
    result['solver_objective'] = solver.objective_value
    result['anchor_parities'] = [solver.value(c) for c in colors]
    def add_class(selected):
        es = [int(i in selected) for i in range(len(ps)+1)]
        rs = [selected.get(i,0) for i in range(len(ps)+1)]
        a,n = crt_prefix(es,rs,[2]+ps)
        result['rows'].append({'exponents':es,'residues':rs})
        result['congruences'].append({'residue':a,'modulus':n})
    if minimum == 2:
        add_class({0:0})
    for i,p in enumerate(ps):
        add_class({i+1:0})
        add_class({0:solver.value(colors[i]),i+1:1})
    for (S,layer,a),v in choices.items():
        if solver.value(v):
            selected = {i+1:t for i,t in zip(S,a)}
            if layer >= 0:
                selected[0] = layer
            add_class(selected)
    ns = [c['modulus'] for c in result['congruences']]
    L = result['period']
    assert len(ns) == len(set(ns)) and min(ns) == minimum
    assert all(L % n == 0 for n in ns)
    # Independent direct integer-period check, not using the unit-box model.
    missed = [x for x in range(L) if not any(x % c['modulus'] == c['residue']
                                            for c in result['congruences'])]
    result['direct_period_check'] = {
        'status':'VERIFIED' if not missed else 'UNCOVERED',
        'period_checked':L,'uncovered_count':len(missed),
        'uncovered_density':str(Fraction(len(missed),L)),
        'first_uncovered_residues':missed[:100],
        'all_moduli_distinct_squarefree':True}
    assert len(missed) == round(solver.objective_value)
    result['verification'] = verify(result,node_limit=1000000)
    result['status'] = 'VERIFIED' if not missed else 'UNCOVERED'
    assert result['verification']['status'] == result['status']
    if missed:
        R = missed[0]
        assert all(R % c['modulus'] != c['residue'] for c in result['congruences'])
        result['uncovered_progression'] = {'residue':R,'modulus':L}
    result['elapsed_seconds'] = round(time.monotonic()-start,3)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--largest-prime',type=int,default=11)
    parser.add_argument('--minimum',type=int,choices=[2,3],default=3)
    parser.add_argument('--seconds',type=int,default=60)
    parser.add_argument('--output',required=True)
    a = parser.parse_args()
    d = search(a.largest_prime,a.minimum,a.seconds)
    (ROOT/a.output).write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({k:v for k,v in d.items() if k not in ('rows','congruences')},indent=2))
