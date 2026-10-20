"""Numerical checks of search derivatives; these do not certify a theorem."""
import json
from pathlib import Path
import numpy as np
from global_distortion_search import GlobalProblem


def check():
    w=json.loads(Path('data/sparse-certificate-20000-50000.json').read_text())
    cases=0
    for mode in ('uniform','reachable'):
        for terminal,control in ((80,300),(130,1000),(400,10000)):
            problem=GlobalProblem(w,K=18000,terminal=terminal,control_prime=control,
                knots=4,resolution=32,state_cap=10000,support_mode=mode)
            x=np.asarray([.13,-.24,.06,.15]);f,g,_=problem.evaluate(x,True)
            h=1e-5;eye=np.eye(4)
            fd=np.asarray([(problem.evaluate(x+h*e)[0]-problem.evaluate(x-h*e)[0])/(2*h) for e in eye])
            assert np.allclose(g,fd,rtol=3e-4,atol=2e-6),(mode,g,fd)
            print(mode,terminal,'max derivative error',np.max(np.abs(g-fd)))
            cases+=4
    print('PASS:',cases,'adjoint derivative comparisons.')
    problem=GlobalProblem(w,K=18000,terminal=8000,control_prime=1000,
        knots=2,resolution=32,state_cap=10000,support_mode='reachable',
        stopping_method='third_moment')
    x=np.asarray([.13,-.24]);_,g,_=problem.evaluate(x,True);h=1e-5
    fd=np.asarray([(problem.evaluate(x+h*e)[0]-problem.evaluate(x-h*e)[0])/(2*h) for e in np.eye(2)])
    assert np.allclose(g,fd,rtol=3e-4,atol=2e-6),(g,fd)
    print('PASS: 2 third-moment terminal-objective derivative comparisons; error',np.max(np.abs(g-fd)))


if __name__=='__main__':check()
