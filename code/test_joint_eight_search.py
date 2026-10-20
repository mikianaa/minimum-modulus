"""Check the early/late pullback with eight-prime conditioning."""
import json
import argparse
from pathlib import Path
import numpy as np
from joint_distortion_search import JointProblem


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--forward',action='store_true');ap.add_argument('--complement',action='store_true');ap.add_argument('--stratified',action='store_true');args=ap.parse_args()
    w=json.loads(Path('data/joint-eight-optimized-seed.json').read_text())
    w['exact_early_prime']=53
    if args.complement:w['conditioned_complement_prime']=859
    if args.stratified:w.update(search_stratified_tail=True,moment_numerator=5,moment_denominator=4,stratified_tail_orders=[[5,4],[10,4],[12,4],[16,4]])
    if args.forward:w['prefix_difference']='forward'
    p=JointProblem(w,terminal=8000,freeze=150,state_cap=1000000 if args.complement else 10000,resolution=32,early_knots=2,knots=2)
    x=np.full(p.dimension,.02);_,g,_=p.evaluate(x,True);h=5e-5;errors=[]
    for j in (3,7,16):
        e=np.eye(p.dimension)[j]
        fd=(p.evaluate(x+h*e)[0]-p.evaluate(x-h*e)[0])/(2*h)
        assert abs(fd-g[j])<2e-5,(j,fd,g[j]);errors.append(abs(fd-g[j]))
    print('PASS: 3 eight-prime joint-product derivative checks;',p.prefix_difference,'maximum error',max(errors))
