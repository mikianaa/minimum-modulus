"""Reoptimize one common future distortion from a completed base report.

The finite distribution is not reconstructed here. As with prefix-report
composition, the caller must first complete the base witness verification.
Every tested future coefficient and every sum below uses integer intervals.
"""
import argparse, hashlib, json
from pathlib import Path
from fractional_moment_tail import fractional_tail_bound,dyadic_power
from verify_minimum_modulus import I,ONE,S,ceildiv,fixed_text
from verify_prefix_refinement import scaled


def prepare(report):
    """Rebuild per-bin moment enclosures once from the certified integer state."""
    state=report.get('stratified_state')
    if state is None:return None
    mass=list(map(int,state['mass']));moment=list(map(int,state['moment']))
    ep=int(state['mass_error']);em=int(state['moment_error']);table=[]
    for row in report['stratified_orders']:
        num,den=map(int,row['order'].split('/'));bounds=[]
        for b,u,P,M in zip(state['bins'],state['upper'],mass,moment):
            up=dyadic_power(I.rat(u),num,den).hi
            value=ceildiv(up*P,S)
            if u>b:
                bp=dyadic_power(I.rat(b),num,den).hi
                aa=max(0,u*P-M);bb=max(0,M+em-b*max(0,P-ep))
                value=min(value,ceildiv(bp*aa+up*bb,(u-b)*S))
            bounds.append(value)
        table.append(bounds)
    return table


def evaluate(report,D,table=None):
    P=report['terminal_prime'];k=report['terminal_index'];total=0
    if table is None:
        for group in report['stratified_selected_groups']:
            num,den=group['order']
            b=fractional_tail_bound(P,k,ONE,D,num,den)[0]
            total+=ceildiv(b*scaled(group['moment_upper']),S)
    else:
        coefficients=[]
        for row in report['stratified_orders']:
            num,den=map(int,row['order'].split('/'))
            coefficients.append(fractional_tail_bound(P,k,ONE,D,num,den)[0])
        total=sum(min(ceildiv(b*u,S) for b,u in zip(coefficients,cell)) for cell in zip(*table))
    available=[]
    for row in report['stratified_orders']:
        if row['moment_upper'] is None:continue
        num,den=map(int,row['order'].split('/'))
        b=fractional_tail_bound(P,k,ONE,D,num,den)[0]
        outside=max(0,scaled(row['moment_upper'])-scaled(row['retained_moment_lower']))
        available.append(ceildiv(b*outside,S))
    assert available
    return total+min(available)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('witness');ap.add_argument('report');ap.add_argument('--output',required=True)
    ap.add_argument('--lower',type=int,default=250);ap.add_argument('--upper',type=int,default=500)
    ap.add_argument('--step',type=int,default=5)
    a=ap.parse_args();wr=Path(a.witness).read_bytes();rr=Path(a.report).read_bytes()
    w=json.loads(wr);r=json.loads(rr)
    assert r['status']=='VERIFIED' and r['witness_sha256']==hashlib.sha256(wr).hexdigest()
    assert r['cutoff']==w['K'] and r['terminal_index']==w['start']+len(w['deltas'])
    assert 0<a.lower<=a.upper<=500 and a.step>0
    assert sum(scaled(g['tail_upper']) for g in r['stratified_selected_groups'])==scaled(r['stratified_retained_tail_upper'])
    table=prepare(r)
    if table is not None:
        assert evaluate(r,w['tail_delta_numerator'],table)==scaled(r['stratified_retained_tail_upper'])+scaled(r['stratified_overflow_tail_upper'])
    rows=[]
    for d in range(a.lower,a.upper+1,a.step):
        tail=evaluate(r,d*10**9,table)
        rows.append(dict(tail_delta_numerator=d*10**9,tail_upper=fixed_text(tail)))
    best=min(rows,key=lambda x:scaled(x['tail_upper']))
    original=scaled(r['tail_loss_upper']);tail=min(original,scaled(best['tail_upper']))
    margin=scaled(r['mu_lower'])-tail;assert margin>0
    result=dict(status='VERIFIED_TAIL_REFINEMENT_FROM_COMPLETED_REPORT',cutoff=w['K'],
                minimum_modulus_at_most=w['K']-1,decimal_places=40,
                base_witness=a.witness,base_report=a.report,
                base_witness_sha256=hashlib.sha256(wr).hexdigest(),base_report_sha256=hashlib.sha256(rr).hexdigest(),
                original_tail_upper=fixed_text(original),best_tested=best,
                tail_loss_upper=fixed_text(tail),infinite_survivor_lower=fixed_text(margin),
                finite_computation_replayed=False,bin_moment_orders_reoptimized=table is not None,rows=rows,
                note='Conditional on the separately completed base verification. One common future distortion; no adaptive bin-specific parameters and no global optimality claim.')
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
