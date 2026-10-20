"""Experimental distribution-stratified tail, with independently bounded overflow."""
import argparse,json
from pathlib import Path
import numpy as np
from frozen_prefix_search import FrozenPrefixProblem
from fractional_moment_tail import fractional_moment,fractional_tail_bound
from verify_minimum_modulus import primes_to,ONE,S


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness');ap.add_argument('--output',required=True)
    ap.add_argument('--terminal',type=int,default=125000);ap.add_argument('--K',type=int)
    ap.add_argument('--tail-deltas',type=int,nargs='+',default=[500000000000])
    ap.add_argument('--witness-output')
    a=ap.parse_args();w=json.loads(Path(a.witness).read_text())
    p=FrozenPrefixProblem(w,K=a.K,terminal=a.terminal,freeze=150,knots=2,tail_delta_numerator=500000000000)
    _,record=p.evaluate(np.zeros(2),save=True);P,M=p.terminal_distribution
    ps=primes_to(max(2000000,20*a.terminal));prepared=[]
    for num,den in ((5,4),(6,4),(7,4),(8,4),(9,4),(5,2),(21,8),(43,16),(11,4),(23,8),(3,1),(13,4),(14,4),(15,4),(4,1)):
        while den<4:num*=2;den*=2
        s=num/den
        U=fractional_moment(ps,record['start'],record['deltas'],num,den).hi/S
        means=np.divide(M,P,out=np.zeros_like(M),where=P>0)
        retained=float(np.maximum(p.values,means)**(s-1)@M)
        outside=max(0.,U-retained)
        width=p.upper-p.values
        aa=np.divide(p.upper*P-M,width,out=np.zeros_like(P),where=width>0)
        bb=np.divide(M-p.values*P,width,out=np.zeros_like(P),where=width>0)
        upper=np.where(width>0,p.values**s*np.maximum(0,aa)+p.upper**s*np.maximum(0,bb),p.values**s*P)
        prepared.append((num,den,U,retained,outside,upper))
    scan=[];best=None
    for D in a.tail_deltas:
        rows=[];candidates=[];overflow=[]
        for num,den,U,retained,outside,upper in prepared:
            coefficient=fractional_tail_bound(ps[a.terminal-1],a.terminal,ONE,D,num,den)[0]/S
            candidates.append(coefficient*upper);overflow.append(coefficient*outside)
            rows.append(dict(order=f'{num}/{den}',moment=U,retained_lower=retained,overflow_tail=overflow[-1],unstratified_tail=coefficient*U))
        tail=float(np.min(candidates,axis=0).sum())+min(overflow)
        scan.append(dict(tail_delta_numerator=D,tail=tail,margin=record['mu']-tail))
        if best is None or tail<best[0]:best=(tail,D,rows)
    tail,D,rows=best
    result=dict(status='EXPERIMENTAL',K=p.K,source=a.witness,terminal=a.terminal,mu=record['mu'],stratified_tail=tail,margin=record['mu']-tail,tail_delta_numerator=D,tail_parameter_scan=scan,rows=rows)
    if a.witness_output:
        record.update(stratified_tail_orders=[[v[0],v[1]] for v in prepared],tail_delta_numerator=D,stratified_tail=tail,stratified_margin=record['mu']-tail)
        Path(a.witness_output).write_text(json.dumps(record)+'\n')
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
