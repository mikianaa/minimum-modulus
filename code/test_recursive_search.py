"""Check the floating atom backend and its joint-search pullback."""
import json
from pathlib import Path
import numpy as np
from recursive_atom_bound import RecursiveAtoms
from recursive_atom_experiment import FloatRecursiveAtoms
from verify_minimum_modulus import S, primes_to
from joint_distortion_search import JointProblem


if __name__=='__main__':
    checks=0
    for K,degree in ((30,1),(100,1),(1000,1),(30,2),(100,2),(1000,2)):
        a=RecursiveAtoms(K,128,degree);b=FloatRecursiveAtoms(K,128,degree)
        for p in primes_to(43):
            x=a.moment_gains(p);y=b.gains(p)
            assert max(abs(x[r]/S-y[r-1]) for r in (1,2))<1e-12
            D=0 if p<5 else 230000000000
            a.update(p,D);b.update(p,D/1e12)
            assert max(abs(v/S-b.F[n]) for n,v in a.F.items())<1e-12
            checks+=1
    w=json.loads(Path('data/joint-eight-optimized-seed.json').read_text())
    w.update(exact_early_prime=53,conditioned_complement_prime=859,
             recursive_atom_cap=128,recursive_atom_prime=859,recursive_density_degree=2,
             prefix_difference='forward',individual_prefix_count=16)
    p=JointProblem(w,terminal=8000,freeze=150,state_cap=1000000,
                   resolution=32,early_knots=2,knots=2)
    x=np.full(p.dimension,.02);_,g,_=p.evaluate(x,True);h=5e-5;errors=[]
    for j in (3,7,16):
        e=np.eye(p.dimension)[j]
        fd=(p.evaluate(x+h*e)[0]-p.evaluate(x-h*e)[0])/(2*h)
        assert abs(fd-g[j])<2e-5,(j,fd,g[j]);errors.append(abs(fd-g[j]))
    print('Recursive search:',checks,'floating/integer comparisons, 3 joint derivatives; max error',max(errors))
