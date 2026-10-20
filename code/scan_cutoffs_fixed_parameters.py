"""Experimental cutoff transfer at fixed parameters, never a certificate.

After every tested cutoff, the finite losses and terminal tail are independent
of K. Thus a short common-horizon survivor difference transfers to the saved
long-horizon margin. Each row still requires independent integer verification.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from frozen_prefix_search import FrozenPrefixProblem


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source')
    ap.add_argument('--cutoffs', type=int, nargs='+', required=True)
    ap.add_argument('--terminal', type=int, default=8000)
    ap.add_argument('--freeze', type=int, default=150)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    raw = Path(a.source).read_bytes(); w = json.loads(raw)
    assert w.get('search_tail_weight', 1) == 1
    assert w['terminal_index'] > a.terminal
    rows = []; base = None
    short_witness = dict(w, search_stratified_tail=False)
    # Only the finite survivor is used. High tail orders valid at the long
    # endpoint may be inadmissible here; no short tail is being certified.
    for K in dict.fromkeys([w['K'], *a.cutoffs]):
        p = FrozenPrefixProblem(short_witness, K=K, terminal=a.terminal, freeze=a.freeze,
                                knots=2, tail_delta_numerator=w['tail_delta_numerator'])
        assert p.primes[-1] >= max([w['K'], *a.cutoffs])
        _, r = p.evaluate(np.zeros(2))
        if base is None: base = r['mu']
        change = r['mu'] - base
        row = dict(K=K, short_mu=r['mu'], survivor_change=change,
                   predicted_long_margin=-w['objective']+change)
        rows.append(row); print(json.dumps(row), flush=True)
        result = dict(status='EXPERIMENTAL', source=a.source,
                      source_sha256=hashlib.sha256(raw).hexdigest(),
                      short_terminal=a.terminal, long_terminal=w['terminal_index'],
                      note='Fixed parameters; cutoff transfer only. No impossibility or optimality claim.',
                      rows=rows)
        Path(a.output).write_text(json.dumps(result, indent=2)+'\n')
