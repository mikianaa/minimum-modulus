"""Reproducible low-dimensional numerical search; results are experimental."""
import json
from pathlib import Path
from minimum_modulus_experiment import run


def optimize(K=615000,passes=2):
    blocks=[[400,1.],[800,1.],[2000,1.],[10000,1.],[100000,1.],[630000,1.]]
    trials=[]
    def score(j,x):
        blocks[j][1]=x
        r=run(K=K,blocks=blocks,limit=630000)
        cs=r.get('checkpoints',[])
        value=cs[-1]['f'] if cs and cs[-1]['i']==51000 else 1e100
        trials.append(dict(blocks=[b[:] for b in blocks],f51000=value))
        return value
    for sweep in range(passes):
        for j in range(len(blocks)):
            center=blocks[j][1];lo=max(.8,center-.06);hi=min(1.2,center+.06)
            g=(5**.5-1)/2
            x=hi-g*(hi-lo);y=lo+g*(hi-lo)
            fx=score(j,x);fy=score(j,y)
            for _ in range(11):
                if fx<fy:
                    hi,y,fy=y,x,fx;x=hi-g*(hi-lo);fx=score(j,x)
                else:
                    lo,x,fx=x,y,fy;y=lo+g*(hi-lo);fy=score(j,y)
            blocks[j][1]=x if fx<fy else y
            print(sweep,j,blocks[j][1],min(fx,fy),flush=True)
    result=dict(K=K,blocks=blocks,trials=trials)
    Path('data/minimum-modulus-optimization.json').write_text(json.dumps(result,indent=2))
    return result


if __name__=='__main__':optimize()
