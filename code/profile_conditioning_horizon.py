"""Experimental terminal profiling; output is not a proof certificate.

Keep the supplied early parameters, extend missing parameters by 0.4, and
evaluate the existing finite-conditioning model at the requested horizon.
Use extract_sparse_witness.py and the integer verifier to certify a result.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from frozen_prefix_search import FrozenPrefixProblem


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('witness')
    ap.add_argument('--K', type=int)
    ap.add_argument('--terminal', type=int, default=150000)
    ap.add_argument('--freeze', type=int, default=100)
    ap.add_argument('--resolution', type=int)
    ap.add_argument('--state-cap', type=int)
    ap.add_argument('--tail-delta-numerator', type=int, default=493000000000)
    ap.add_argument('--output', required=True)
    args = vars(ap.parse_args())
    witness = json.loads(Path(args.pop('witness')).read_text())
    output = args.pop('output')
    problem = FrozenPrefixProblem(witness, knots=2, **args)
    _, result = problem.evaluate(np.zeros(2), save=True)
    result['status'] = 'EXPERIMENTAL'
    Path(output).write_text(json.dumps(result) + '\n')
    print(json.dumps({key: result[key] for key in
                     ('K', 'mu', 'minimum_mu', 'objective', 'terminal_checkpoints')}, indent=2))
