"""Run the independent checks; this does not certify a modulus cutoff."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--with-search',action='store_true',help='Also check NumPy search derivatives.')
    ap.add_argument('--output',default='data/conditioning-method-checks.json')
    args=ap.parse_args();root=Path(__file__).resolve().parent.parent
    names=['test_dynamic_joint.py','test_statewise_presence.py','test_early_exact.py',
           'test_fractional_tail.py','test_third_moment_tail.py','test_joint_product.py','test_conditioned_complement.py','test_stratified_tail.py','test_prefix_refinement.py','test_recursive_atoms.py','test_hybrid_prefix.py','test_product_truncation.py','test_block_tail.py']
    if args.with_search:names.extend(['test_recursive_search.py','test_first_prime_search.py','test_joint_search.py','test_joint_eight_search.py','test_product_search.py',('test_joint_eight_search.py','--forward','--complement','--stratified')])
    rows=[]
    for name in names:
        argv=list(name) if isinstance(name,tuple) else [name]
        run=subprocess.run([sys.executable,str(root/'code'/argv[0]),*argv[1:]],cwd=root,text=True,capture_output=True)
        rows.append(dict(file='code/'+argv[0],arguments=argv[1:],returncode=run.returncode,stdout=run.stdout,stderr=run.stderr))
        print(run.stdout,end='',flush=True)
        if run.returncode:print(run.stderr,file=sys.stderr);break
    result=dict(status='PASSED' if len(rows)==len(names) and all(r['returncode']==0 for r in rows) else 'FAILED',
                note='Independent implementation checks, not a modulus certificate or formal proof.',
                checks=rows,source_sha256={str(p.relative_to(root)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in sorted((root/'code').glob('*.py'))})
    (root/args.output).write_text(json.dumps(result,indent=2)+'\n')
    if result['status']!='PASSED':sys.exit(1)
