"""Extract only rational proof inputs, optionally shortening a search run."""
import argparse,json
from pathlib import Path

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('output')
    ap.add_argument('--terminal',type=int,required=True)
    ap.add_argument('--tail-delta-numerator',type=int);a=ap.parse_args()
    w=json.loads(Path(a.input).read_text())
    assert w['start']<a.terminal<=w['start']+len(w['deltas'])
    result={k:w[k] for k in ('K','start','small','state_cap','cap','resolution')}
    for key in ('support_mode','stopping_method','tail_delta_numerator','normalize_conditioning','interpolate_tail','joint_presence','statewise_presence','joint_product','dynamic','shifted_tail','exact_early_prime','conditioned_complement_prime','moment_numerator','moment_denominator','recursive_atom_cap','recursive_atom_prime','recursive_density_degree'):
        if key in w:result[key]=w[key]
    result['deltas']=w['deltas'][:a.terminal-w['start']]
    if 'stratified_tail_orders' in w:result['stratified_tail_orders']=w['stratified_tail_orders']
    if 'stratified_overflow_orders' in w:result['stratified_overflow_orders']=w['stratified_overflow_orders']
    for key in ('prefix_resolution','prefix_index','product_tail_power'):
        if key in w:result[key]=w[key]
    if a.tail_delta_numerator is not None:
        assert 0<a.tail_delta_numerator<=500000000000
        assert result.get('stopping_method') in ('third_moment','fractional_moment')
        result['tail_delta_numerator']=a.tail_delta_numerator
    Path(a.output).write_text(json.dumps(result)+'\n')
