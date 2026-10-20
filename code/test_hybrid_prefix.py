"""One-pass hybrid must equal the exact sum of independent finite runs."""
import json
import tempfile
from pathlib import Path
from verify_sparse_conditioning import verify
from verify_prefix_refinement import scaled


if __name__=='__main__':
    scratch=(Path(__file__).resolve().parent.parent/'data').resolve()
    source=json.loads((scratch/'sparse-certificate-14498.json').read_text())
    source.update(K=20000,resolution=16,cap=100000,state_cap=10000,
                  exact_early_prime=0,conditioned_complement_prime=0,
                  recursive_atom_cap=128,recursive_atom_prime=97,recursive_density_degree=2)
    source['deltas']=source['deltas'][:98]
    for key in ('stratified_tail_orders','stratified_overflow_orders'):
        source.pop(key,None)
    with tempfile.TemporaryDirectory(prefix='hybrid-check-',dir=scratch) as tmp:
        root=Path(tmp);assert root.resolve().parent==scratch
        def check(w,name):
            p=root/(name+'.json');p.write_text(json.dumps(w))
            return verify(p,root/(name+'-result.json'),finite_only=True)
        coarse=check(source,'coarse')
        old=dict(source,deltas=source['deltas'][:48]);oldr=check(old,'old')
        fine=dict(old,resolution=64);finer=check(fine,'fine')
        hybrid=dict(source,prefix_index=50,prefix_resolution=64)
        result=check(hybrid,'hybrid')
        assert scaled(result['mu_lower'])==scaled(coarse['mu_lower'])-scaled(oldr['mu_lower'])+scaled(finer['mu_lower'])
        assert result['status']=='VERIFIED_FINITE_PREFIX'
    print('Hybrid prefix: exact equality with four independent finite runs, including recursive caps.')
