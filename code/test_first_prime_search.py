"""Starting at prime 2 agrees at zero tilt and has the correct derivatives."""
import json,tempfile
from pathlib import Path
import numpy as np
from joint_distortion_search import JointProblem
from verify_sparse_conditioning import verify
from verify_prefix_refinement import scaled


if __name__=='__main__':
    w=json.loads(Path('data/joint-eight-optimized-seed.json').read_text())
    w.update(exact_early_prime=53,conditioned_complement_prime=859,prefix_difference='central')
    args=dict(terminal=8000,freeze=150,state_cap=1000000,resolution=32,early_knots=2,knots=2)
    a=JointProblem(w,start=2,**args);b=JointProblem(w,start=0,**args)
    a.set_prefix(np.zeros(a.early_count));b.set_prefix(np.zeros(b.early_count))
    assert abs(a.mu0-b.mu0)<1e-12
    assert abs(a.C0-b.C0)<1e-12 and abs(a.U0-b.U0)<1e-8
    assert all(np.max(np.abs(x-y))<1e-12 for x,y in zip(a.initial,b.initial))
    x=np.zeros(b.dimension);x[0]=.025;x[1]=.015
    _,g,_=b.evaluate(x,True);h=5e-5;errors=[]
    for j in (0,1):
        e=np.eye(b.dimension)[j]
        fd=(b.evaluate(x+h*e)[0]-b.evaluate(x-h*e)[0])/(2*h)
        assert abs(fd-g[j])<2e-5,(j,fd,g[j]);errors.append(abs(fd-g[j]))
    scratch=(Path(__file__).resolve().parent.parent/'data').resolve()
    with tempfile.TemporaryDirectory(prefix='first-prime-check-',dir=scratch) as tmp:
        root=Path(tmp);assert root.resolve().parent==scratch
        base=dict(K=20000,start=2,small=2,dynamic=True,state_cap=1000,cap=100,resolution=8,deltas=[0]*3)
        results=[]
        for name,v in [('old',base),('zero',dict(base,start=0,deltas=[0]*5)),
                       ('tilted',dict(base,start=0,deltas=[100000000,70000000,0,0,0],exact_early_prime=11,recursive_atom_cap=128,recursive_atom_prime=11,recursive_density_degree=2))]:
            p=root/(name+'.json');p.write_text(json.dumps(v));results.append(verify(p,root/(name+'-report.json'),True))
        assert abs(scaled(results[0]['mu_lower'])-scaled(results[1]['mu_lower']))<10**10
        assert all(r['status']=='VERIFIED_FINITE_PREFIX' for r in results)
    print('First-prime search: zero-tilt identity, 2 derivatives, 3 integer-prefix checks; max error',max(errors))
