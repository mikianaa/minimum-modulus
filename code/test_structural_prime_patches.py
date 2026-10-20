"""Exhaustive checks of the fresh-prime projection identity on small kernels."""
from itertools import product
from math import lcm
import random

from construction_kernel import matches, verify
from structural_prime_patches import divisors, factor_over, crt_prefix, capacity_audit


def main():
    rng = random.Random(100)
    ps,ds = [2,3],[4,2]
    checks = 0
    for _ in range(60):
        blocks = []
        rows = []
        for j,q in enumerate([5,7]):
            ns = sorted(rng.sample(divisors([3,1],ps),q-1))
            m = lcm(*ns)
            a = rng.randrange(m)
            ex = factor_over(m,ps)
            assert crt_prefix(ex,[a % p**e for p,e in zip(ps,ex)],ps) == (a,m)
            blocks.append((a,m))
            for digit,n in enumerate(ns,1):
                es = factor_over(n,ps)
                rows.append({'exponents':es+[int(j == 0),int(j == 1)],
                             'residues':[a % p**e for p,e in zip(ps,es)]+
                                 [digit-1 if j == 0 else 0,digit-1 if j == 1 else 0]})
        # The first nonzero tail digit is immaterial to all fixed divisors.
        for s2,s3 in product(range(8),range(6)):
            a2,a3 = s2,s3//2
            x,_ = crt_prefix([3,1],[a2,a3],ps)
            predicted = all(x % m != a for a,m in blocks)
            actual = any(not any(matches([s2,s3,u,v],r,ps+[5,7],ds+[1,1]) for r in rows)
                         for u,v in product(range(4),range(6)))
            assert actual == predicted
            checks += 1
    # An explicit positive control: parity plus a root-3 divisor block.
    cert = {'minimum':2,'primes':[2,3],'depths':[2,1], 'rows':[
        {'exponents':[1,0],'residues':[1,0]},
        {'exponents':[0,1],'residues':[0,0]},
        {'exponents':[1,1],'residues':[0,1]}]}
    assert verify(cert)['status'] == 'VERIFIED'
    for minimum in range(2,16):
        b = {'minimum':minimum,'primes':[2,3],'depths':[4,3],
             'exact_uncovered_states':1,'kernel_states':1}
        audit = capacity_audit(b)
        for entry in audit['root_capacities']:
            q = entry['root_prime']
            brute = min(2**a*3**b for a,b in product(range(4),range(3))
                        if sum(n*q >= minimum for n in divisors([a,b],[2,3])) >= q-1)
            assert brute == entry['smallest_cylinder_modulus']
    print({'status':'PASS','exhaustive_projected_state_checks':checks,
           'finite_kernel_positive_control':2,'capacity_audits_checked':14})


if __name__ == '__main__':
    main()
