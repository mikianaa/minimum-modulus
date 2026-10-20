"""Screen survivor-cost policies and stronger conditioning; experimental only."""
import json,time
from pathlib import Path
from strong_conditioning_experiment import run

results=[]
for es,K in [((26,17,11,9),25000),((22,14,9,7,5,5),22000)]:
    for policy in (.65,.8,.9,1.,1.1,1.3):
        t=time.time()
        r=run(K=K,start=20,policy=policy,max_steps=30000,cap=1000000,
              small_exponents=es,resolution=256)
        r.update(seconds=time.time()-t,policy=policy)
        results.append(r)
        Path('data/strong-policy-sweep.json').write_text(json.dumps(results,indent=2))
        print(K,len(es),policy,r['success'],r.get('terminal',{}).get('ratio'),
              r.get('p'),r['best_ratio'],round(r['seconds'],2),flush=True)
