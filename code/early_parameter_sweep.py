"""Experimental coordinated perturbations of the fixed early prefix.

This finite sweep is not a global optimality claim or a proof certificate.
The later parameters are fixed, and the complete terminal loss is evaluated.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from frozen_prefix_search import FrozenPrefixProblem


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('witness')
    ap.add_argument('--K', type=int, required=True)
    ap.add_argument('--terminal', type=int, default=50000)
    ap.add_argument('--freeze', type=int, default=100)
    ap.add_argument('--shifts', nargs='+', type=float, default=[-.5, -.25, .25, .5])
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    base = json.loads(Path(args.witness).read_text())
    rows = []
    for shift in args.shifts:
        witness = dict(base)
        witness['deltas'] = list(base['deltas'])
        for j in range(args.freeze - witness['start']):
            d = witness['deltas'][j] / 1e12
            if 0 < d < .5:
                z = math.log(2*d/(1-2*d)) + shift
                witness['deltas'][j] = round(1e12*.5/(1+math.exp(-z)))
        problem = FrozenPrefixProblem(witness, K=args.K, terminal=args.terminal,
                                      freeze=args.freeze, knots=2)
        objective, result = problem.evaluate(np.zeros(2))
        row = dict(shift=shift, mu=result['mu'], objective=objective,
                   tail=result['terminal_cost'])
        rows.append(row)
        print(row, flush=True)
        Path(args.output).write_text(json.dumps(dict(status='EXPERIMENTAL',
            K=args.K, terminal=args.terminal, freeze=args.freeze,
            source=args.witness, rows=rows), indent=2) + '\n')
        del problem
