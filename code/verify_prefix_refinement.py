"""Replay and combine a refined finite prefix with an unchanged certified tail.

Default mode independently recomputes all three component certificates.
--combine-completed-reports only checks compatibility, hashes and arithmetic
of existing component reports; use it only after running their verifiers.
"""
import argparse
import hashlib
import json
from decimal import Decimal,localcontext
from pathlib import Path
from verify_minimum_modulus import S,fixed_text,primes_to
from verify_sparse_conditioning import verify


def scaled(value):
    with localcontext() as ctx:
        ctx.prec=100
        v=Decimal(value)*S
        assert v==int(v)
        return int(v)


def combine(spec):
    inputs={};reports={};hashes={}
    for role in ('base','old_prefix','new_prefix'):
        raw=Path(spec[role]['witness']).read_bytes();w=json.loads(raw)
        r=json.loads(Path(spec[role]['report']).read_bytes())
        assert r['witness_sha256']==hashlib.sha256(raw).hexdigest()
        assert r['cutoff']==w['K'] and r['terminal_index']==w['start']+len(w['deltas'])
        n=r['terminal_index'];limit=max(100,20*n);ps=primes_to(limit)
        while len(ps)<n:limit*=2;ps=primes_to(limit)
        assert r['terminal_prime']==ps[n-1]
        assert r['decimal_places']==40 and scaled(r['mu_lower'])>0
        assert r['status']==('VERIFIED' if role=='base' else 'VERIFIED_FINITE_PREFIX')
        inputs[role]=w;reports[role]=r
        hashes[role]=hashlib.sha256(Path(spec[role]['report']).read_bytes()).hexdigest()
    base=inputs['base'];old=inputs['old_prefix'];new=inputs['new_prefix']
    expected=dict(base);expected['deltas']=base['deltas'][:len(old['deltas'])]
    assert old==expected, 'The subtracted prefix must reproduce the base inputs exactly.'
    assert len(old['deltas'])<len(base['deltas']) and new['K']<=base['K']
    expected=dict(old);expected.update(K=new['K'],resolution=new['resolution'])
    assert new==expected, 'Only cutoff and bin resolution may change in the refined prefix.'
    assert reports['old_prefix']['terminal_prime']==reports['new_prefix']['terminal_prime']
    assert reports['old_prefix']['terminal_prime']>=base['K']
    b=reports['base'];o=reports['old_prefix'];f=reports['new_prefix']
    tail=scaled(b['tail_loss_upper'])
    assert scaled(b['infinite_survivor_lower'])==scaled(b['mu_lower'])-tail
    block_details={}
    if 'block_tail' in spec:
        settings=spec['block_tail'];raw=Path(settings['report']).read_bytes();t=json.loads(raw)
        assert t['status']=='VERIFIED_BLOCK_TAIL_FROM_COMPLETED_REPORT'
        assert t['base_witness_sha256']==b['witness_sha256']
        assert t['base_report_sha256']==hashes['base']
        assert t['cutoff']==base['K'] and t['base_index']==b['terminal_index']
        assert t['mu_lower']==b['mu_lower'] and t['decimal_places']==40
        for key in ('terminal_index','block_size','future_delta_numerator','exponent_truncation'):
            assert t[key]==settings[key]
        block_tail=scaled(t['tail_loss_upper'])
        assert scaled(t['infinite_survivor_lower'])==scaled(b['mu_lower'])-block_tail>0
        tail=min(tail,block_tail)
        block_details=dict(block_tail_report=settings['report'],block_tail_report_sha256=hashlib.sha256(raw).hexdigest(),
                           block_tail_terminal_index=t['terminal_index'],block_tail_terminal_prime=t['terminal_prime'])
    adjustment=scaled(f['mu_lower'])-scaled(o['mu_lower'])
    survivor=scaled(b['mu_lower'])+adjustment;gap=survivor-tail
    assert gap>0, ('refined prefix does not close the infinite tail',fixed_text(gap))
    return dict(status='VERIFIED_COMPOSITION',cutoff=new['K'],minimum_modulus_at_most=new['K']-1,
        decimal_places=40,method='refined finite prefix with independently computed base suffix and tail',
        terminal_index=b['terminal_index'],terminal_prime=b['terminal_prime'],
        prefix_index=o['terminal_index'],prefix_prime=o['terminal_prime'],
        prefix_survivor_adjustment=fixed_text(adjustment),mu_lower=fixed_text(survivor),
        tail_loss_upper=fixed_text(tail),infinite_survivor_lower=fixed_text(gap),
        component_report_sha256=hashes,
        note='The mathematical conclusion uses all independently verified components, including the optional block tail. Their fresh replay is the default command.',**block_details)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('witness')
    ap.add_argument('--output',required=True)
    ap.add_argument('--combine-completed-reports',action='store_true')
    a=ap.parse_args();raw=Path(a.witness).read_bytes();spec=json.loads(raw)
    if not a.combine_completed_reports:
        for role in ('base','old_prefix','new_prefix'):
            verify(spec[role]['witness'],spec[role]['report'],finite_only=(role!='base'))
        if 'block_tail' in spec:
            from block_tail_bound import verify as verify_block_tail
            t=spec['block_tail']
            r=verify_block_tail(spec['base']['witness'],spec['base']['report'],t['terminal_index'],
                                t['block_size'],t['future_delta_numerator'],t['exponent_truncation'])
            Path(t['report']).write_text(json.dumps(r,indent=2)+'\n')
    result=combine(spec)
    result['witness_sha256']=hashlib.sha256(raw).hexdigest()
    result['components_replayed_in_this_run']=not a.combine_completed_reports
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
