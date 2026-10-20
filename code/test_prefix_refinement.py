"""Check prefix-composition arithmetic, compatibility and finite-only status."""
import copy
import hashlib
import json
import tempfile
from pathlib import Path
from verify_prefix_refinement import combine
from verify_sparse_conditioning import verify


if __name__=='__main__':
    scratch=(Path(__file__).resolve().parent.parent/'data').resolve()
    with tempfile.TemporaryDirectory(prefix='prefix-check-',dir=scratch) as tmp:
        root=Path(tmp)
        assert root.resolve().parent==scratch
        base=dict(K=20,start=2,deltas=[0]*12,resolution=128)
        old=dict(base,deltas=[0]*8);new=dict(old,K=19,resolution=512)
        witnesses=dict(base=base,old_prefix=old,new_prefix=new)
        reports=dict(base=dict(status='VERIFIED',cutoff=20,terminal_index=14,terminal_prime=43,
                              mu_lower='0.2',tail_loss_upper='0.1',infinite_survivor_lower='0.1'),
                     old_prefix=dict(status='VERIFIED_FINITE_PREFIX',cutoff=20,terminal_index=10,terminal_prime=29,mu_lower='0.4'),
                     new_prefix=dict(status='VERIFIED_FINITE_PREFIX',cutoff=19,terminal_index=10,terminal_prime=29,mu_lower='0.35'))
        spec={}
        def save(ws,rs):
            for role in ws:
                p=root/(role+'.json');p.write_text(json.dumps(ws[role]))
                r=dict(rs[role],witness_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),decimal_places=40)
                q=root/(role+'-report.json');q.write_text(json.dumps(r))
                spec[role]=dict(witness=str(p),report=str(q))
        save(witnesses,reports);r=combine(spec)
        assert r['mu_lower']=='0.1500000000000000000000000000000000000000'
        assert r['infinite_survivor_lower']=='0.0500000000000000000000000000000000000000'
        for case in range(5):
            ws=copy.deepcopy(witnesses);rs=copy.deepcopy(reports)
            if case==0:ws['old_prefix']['deltas'][0]=1
            if case==1:ws['new_prefix']['deltas'][0]=1
            if case==2:rs['old_prefix']['terminal_prime']=19
            if case==3:rs['new_prefix']['mu_lower']='0.2'
            if case==4:rs['base']['infinite_survivor_lower']='0.2'
            save(ws,rs)
            try:combine(spec)
            except AssertionError:pass
            else:raise AssertionError(('invalid composition accepted',case))
        save(witnesses,reports)
        settings=dict(report=str(root/'block.json'),terminal_index=100,block_size=10,
                      future_delta_numerator=400000000000,exponent_truncation=8)
        block=dict(status='VERIFIED_BLOCK_TAIL_FROM_COMPLETED_REPORT',cutoff=20,base_index=14,
                   base_witness_sha256=hashlib.sha256(Path(spec['base']['witness']).read_bytes()).hexdigest(),
                   base_report_sha256=hashlib.sha256(Path(spec['base']['report']).read_bytes()).hexdigest(),
                   mu_lower='0.2',tail_loss_upper='0.08',infinite_survivor_lower='0.12',decimal_places=40,
                   terminal_prime=541,**{k:v for k,v in settings.items() if k!='report'})
        Path(settings['report']).write_text(json.dumps(block));spec['block_tail']=settings
        assert combine(spec)['infinite_survivor_lower']=='0.0700000000000000000000000000000000000000'
        for key,value in [('base_report_sha256','wrong'),('future_delta_numerator',1),('infinite_survivor_lower','0.2')]:
            bad=dict(block);bad[key]=value;Path(settings['report']).write_text(json.dumps(bad))
            try:combine(spec)
            except AssertionError:pass
            else:raise AssertionError(('invalid block composition accepted',key))
        w=dict(K=20000,start=2,small=2,state_cap=1000,cap=100,resolution=8,deltas=[0]*3)
        p=root/'finite.json';p.write_text(json.dumps(w))
        result=verify(p,root/'finite-report.json',finite_only=True)
        assert result['status']=='VERIFIED_FINITE_PREFIX' and 'minimum_modulus_at_most' not in result
    print('PASS: exact prefix and block-tail compositions, eight invalid-composition rejections, and real finite-only verification.')
