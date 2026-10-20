"""Finite-difference check of the post-prefix adjoint optimizer."""
import json
from pathlib import Path
import numpy as np
from frozen_prefix_search import FrozenPrefixProblem


if __name__=='__main__':
    w=json.loads(Path('data/sparse-certificate-15500.json').read_text())
    assert w['start']<w['small'] and any(w['deltas'][:w['small']-w['start']])
    p=FrozenPrefixProblem(w,terminal=8000,state_cap=10000,resolution=32,knots=2)
    x=np.asarray([.1,-.2]);_,g,_=p.evaluate(x,True);h=1e-5
    fd=np.asarray([(p.evaluate(x+h*e)[0]-p.evaluate(x-h*e)[0])/(2*h) for e in np.eye(2)])
    assert np.allclose(g,fd,rtol=5e-4,atol=3e-6),(g,fd)
    _,record=p.evaluate(x,save=True)
    assert record['start']+len(record['deltas'])==8000
    assert record['dynamic'] and record['joint_presence'] and record['interpolate_tail']
    print('PASS: 2 post-prefix adjoint derivatives; parameter serialization; maximum error',max(abs(g-fd)))
