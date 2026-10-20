"""Reproducible parameter-policy screening; floating point only."""
import json,time
from pathlib import Path
from strong_conditioning_experiment import run

results=[]
for K in (10000,20000,30000):
    for start in (6,12,20):
        for policy in (.3,1.,3.,10.):
            t=time.time()
            r=run(K=K,start=start,policy=policy,max_steps=15000,
                  small_exponents=(26,17,11,9),resolution=128)
            r.update(seconds=time.time()-t,policy=policy)
            results.append(r)
            Path('data/strong-conditioning-sweep.json').write_text(json.dumps(results,indent=2))
            print(K,start,policy,r['success'],r.get('terminal',{}).get('ratio'),
                  r.get('p'),r['best_ratio'],round(r['seconds'],2),flush=True)
