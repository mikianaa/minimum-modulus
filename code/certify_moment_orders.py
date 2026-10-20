"""Verify several tail bounds for fixed rational distortion parameters.

This checks the infinite tail only. It is not a stand-alone cutoff
certificate: finite losses still require verify_sparse_conditioning.py.
"""
import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction
from fractional_moment_tail import fractional_moment,fractional_tail_bound
from verify_minimum_modulus import primes_to
from verify_stoploss import L


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled.')
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('witness');ap.add_argument('--output',required=True)
    ap.add_argument('--orders',nargs='+',default=['85/32','43/16','87/32','11/4']);args=ap.parse_args()
    raw=Path(args.witness).read_bytes();w=json.loads(raw);start=w['start'];ds=w['deltas'];n=start+len(ds)
    assert start>=2 and all(isinstance(d,int) and 0<=d<=L//2 for d in ds)
    ps=primes_to(max(2000000,20*n));assert len(ps)>=n
    rows=[]
    for order in args.orders:
        s=Fraction(order);num,den=s.numerator,s.denominator
        # The implementation accepts denominator 4 for integral/half orders.
        while den<4:num*=2;den*=2
        tail,row=fractional_tail_bound(ps[n-1],n,fractional_moment(ps,start,ds,num,den),w['tail_delta_numerator'],num,den)
        rows.append(row);print(row['moment_order'],row['tail_loss_upper'],flush=True)
    result=dict(status='VERIFIED_TAIL_BOUNDS_ONLY',note='Finite losses are not checked by this script.',witness=args.witness,witness_sha256=hashlib.sha256(raw).hexdigest(),terminal_index=n,terminal_prime=ps[n-1],rows=rows)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
