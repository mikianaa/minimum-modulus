"""Independent finite differences of the joint early/late pullback."""
import json
from pathlib import Path
import numpy as np
from joint_distortion_search import JointProblem


if __name__=='__main__':
    w=json.loads(Path('data/sparse-certificate-14900.json').read_text())
    w.update(shifted_tail=True,statewise_presence=True)
    p=JointProblem(w,terminal=8000,state_cap=10000,resolution=32,early_knots=2,knots=2)
    x=np.full(p.dimension,.02);_,g,_=p.evaluate(x,True);h=5e-5;errors=[]
    for j in (0,3,7,13,14,15,16,17):
        e=np.eye(p.dimension)[j]
        fd=(p.evaluate(x+h*e)[0]-p.evaluate(x-h*e)[0])/(2*h)
        assert abs(fd-g[j])<2e-5,(j,fd,g[j])
        errors.append(abs(fd-g[j]))
    _,r=p.evaluate(x,save=True)
    assert len(r['deltas'])+r['start']==8000 and r['statewise_presence'] and r['shifted_tail']
    print('PASS: 8 joint early/late derivative checks; maximum error',max(errors))
    w.update(exact_early_prime=53,stopping_method='fractional_moment',moment_numerator=11)
    p=JointProblem(w,terminal=8000,state_cap=10000,resolution=32,early_knots=2,knots=2)
    x=np.full(p.dimension,.02);_,g,_=p.evaluate(x,True);errors=[]
    for j in (0,7,14,16):
        e=np.eye(p.dimension)[j]
        fd=(p.evaluate(x+h*e)[0]-p.evaluate(x-h*e)[0])/(2*h)
        assert abs(fd-g[j])<2e-5,(j,fd,g[j]);errors.append(abs(fd-g[j]))
    print('PASS: 4 exact-early/fractional/packed derivative checks; maximum error',max(errors))
    p.joint_product=True
    _,g,_=p.evaluate(x,True);errors=[]
    for j in (7,14):
        e=np.eye(p.dimension)[j]
        fd=(p.evaluate(x+h*e)[0]-p.evaluate(x-h*e)[0])/(2*h)
        assert abs(fd-g[j])<2e-5,(j,fd,g[j]);errors.append(abs(fd-g[j]))
    print('PASS: 2 full-product conditioning derivative checks; maximum error',max(errors))
