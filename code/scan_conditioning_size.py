"""Compare larger conditioned prime sets at fixed parameters; discovery only."""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from frozen_prefix_search import FrozenPrefixProblem


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('witness');ap.add_argument('--K',type=int,nargs='+',required=True)
    ap.add_argument('--small',type=int,nargs='+',default=[8,9,10])
    ap.add_argument('--terminal',type=int,default=8000)
    ap.add_argument('--output',required=True)
    a=ap.parse_args();raw=Path(a.witness).read_bytes();base=json.loads(raw);rows=[]
    result=dict(status='EXPERIMENTAL',source=a.witness,source_sha256=hashlib.sha256(raw).hexdigest(),
                terminal_index=a.terminal,rows=rows,
                note='Fixed rational parameters; finite survivor comparisons only. No modulus certificate or global optimization claim.')
    for K in a.K:
        for small in a.small:
            w=dict(base,small=small,search_stratified_tail=False);t=time.time()
            problem=FrozenPrefixProblem(w,K=K,terminal=a.terminal,freeze=150,knots=2,
                                        tail_delta_numerator=500000000000)
            _,record=problem.evaluate(np.zeros(2))
            row=dict(K=K,small=small,conditioning_states=len(problem.sm.ns),
                     state_cap=w['state_cap'],mu=record['mu'],seconds=time.time()-t)
            rows.append(row);Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps(row),flush=True)
