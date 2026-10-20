"""Finite policy screening for the sparse conditioning bound, not proof."""
import json,time
from pathlib import Path
from sparse_conditioning_experiment import run

results=[]
for policy in (.75,.9,1.,1.1,1.25):
    t=time.time()
    r=run(K=20000,start=16,small=16,state_cap=3000000,policy=policy,max_steps=15000)
    r.update(policy=policy,seconds=time.time()-t);results.append(r)
    Path('data/sparse-policy-sweep.json').write_text(json.dumps(results,indent=2))
    print(policy,r['success'],r.get('p'),r.get('terminal',{}).get('ratio'),r['best_ratio'],r['seconds'],flush=True)
