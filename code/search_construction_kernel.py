"""Bounded CP-SAT search; only construction_kernel.verify certifies success.

Install ortools normally or in tmp/construction-deps. Search failures apply
only to the explicitly recorded finite family, never to all coverings.
"""
import argparse
from itertools import product
import json
from math import prod
from pathlib import Path
import random
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tmp'/'construction-deps'))
from ortools.sat.python import cp_model
from construction_kernel import resources, capacity, sizes, project, labels, verify


def search(ps, ds, minimum, seconds, cap, seed, initial, batch):
    start = time.monotonic()
    rng = random.Random(seed)
    exps = list(resources(ps, ds, minimum, cap))
    mass = capacity(ps, ds, exps)
    dims = sizes(ps, ds)
    result = {'minimum': minimum, 'primes': ps, 'depths': ds,
              'maximum_resource_base': cap, 'seed': seed, 'seconds_budget': seconds,
              'kernel_states': prod(dims), 'resource_count': len(exps),
              'capacity': str(mass), 'rows': [], 'iterations': []}
    if mass < 1:
        result['status'] = 'EXCLUDED_BY_EXACT_CAPACITY'
        return result
    model = cp_model.CpModel()
    variables = {}
    equalities = {}
    points = set()
    for j, es in enumerate(exps):
        for i, e in enumerate(es):
            if e:
                variables[j, i] = model.new_int_var(0, labels(ps[i], ds[i], e)-1, f'r{j}_{i}')

    def add_point(point):
        point = tuple(point)
        if point in points:
            return
        points.add(point)
        covers = []
        for j, es in enumerate(exps):
            parts = []
            for i, e in enumerate(es):
                if not e:
                    continue
                r = project(point[i], ps[i], ds[i], e)
                key = (j, i, r)
                if key not in equalities:
                    eq = model.new_bool_var(f'e{j}_{i}_{r}')
                    model.add(variables[j, i] == r).only_enforce_if(eq)
                    model.add(variables[j, i] != r).only_enforce_if(eq.Not())
                    equalities[key] = eq
                parts.append(equalities[key])
            if len(parts) == 1:
                covers.append(parts[0])
            else:
                c = model.new_bool_var(f'c{len(points)}_{j}')
                model.add_bool_and(parts).only_enforce_if(c)
                covers.append(c)
        model.add_bool_or(covers)

    if prod(dims) <= initial:
        for point in product(*(range(n) for n in dims)):
            add_point(point)
    else:
        add_point([0]*len(ps))
        for _ in range(initial):
            add_point([rng.randrange(n) for n in dims])
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 4
    solver.parameters.random_seed = seed
    # Keep each solve short enough for progress reports and cancellation.
    while time.monotonic()-start < seconds:
        solver.parameters.random_seed = seed + len(result['iterations'])
        solver.parameters.max_time_in_seconds = min(30, seconds-(time.monotonic()-start))
        status = solver.solve(model)
        info = {'samples': len(points), 'solver_status': solver.status_name(status),
                'elapsed_seconds': round(time.monotonic()-start, 3)}
        result['iterations'].append(info)
        print(json.dumps(info), flush=True)
        if status == cp_model.INFEASIBLE:
            result['status'] = 'SOLVER_INFEASIBLE_RESTRICTED_FAMILY'
            break
        if status == cp_model.MODEL_INVALID:
            result['status'] = 'MODEL_INVALID'
            result['solver_diagnostic'] = solver.solution_info()
            break
        if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
            result['status'] = 'UNKNOWN_TIME_LIMIT'
            continue
        result['rows'] = [{'exponents': list(es), 'residues': [
            solver.value(variables[j, i]) if e else 0 for i, e in enumerate(es)]}
            for j, es in enumerate(exps)]
        check = verify(result, node_limit=300000)
        result['verification'] = check
        if check['status'] == 'VERIFIED':
            result['status'] = 'VERIFIED'
            break
        if check['status'] == 'UNKNOWN':
            result['status'] = 'UNKNOWN_VERIFICATION_LIMIT'
            break
        add_point(check['point'])
        # Add uncovered samples, not arbitrary previously covered states.
        added = 0
        for _ in range(batch*50):
            point = tuple(rng.randrange(n) for n in dims)
            if point in points:
                continue
            covered = any(all(project(point[i], ps[i], ds[i], e) == row['residues'][i]
                          for i, e in enumerate(row['exponents'])) for row in result['rows'])
            if not covered:
                add_point(point)
                added += 1
                if added >= batch:
                    break
        model.clear_hints()
        for var in variables.values():
            model.add_hint(var, solver.value(var))
        result['status'] = 'INCOMPLETE'
    result['elapsed_seconds'] = round(time.monotonic()-start, 3)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--primes', default='2,3,5,7,11,13')
    parser.add_argument('--depths', default='6,4,3,2,2,2')
    parser.add_argument('--minimum', type=int, default=43)
    parser.add_argument('--seconds', type=float, default=180)
    parser.add_argument('--cap', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--initial', type=int, default=64)
    parser.add_argument('--batch', type=int, default=32)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    ps = list(map(int, args.primes.split(',')))
    ds = list(map(int, args.depths.split(',')))
    if len(ps) != len(ds) or any(E < 1 for E in ds):
        parser.error('Primes and positive depths must have equal length')
    result = search(ps, ds, args.minimum, args.seconds, args.cap,
                    args.seed, args.initial, args.batch)
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'iterations')}, indent=2))
