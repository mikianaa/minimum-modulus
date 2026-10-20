"""Fresh-prime divisor blocks, with exact control of their overlap.

A block at q uses q-1 distinct divisors n_d of m. Its row for tail digit d
is x=a (mod n_d), x=d*q^(k-1) (mod q^k), k>=1. All its q states are covered
exactly when x=a (mod lcm(n_d)). Distinct fresh q's are independent, so a
base state has an uncovered extension iff it misses every such cylinder.
The emitted rows are checked by the independent finite-kernel verifier.
"""
import argparse
from fractions import Fraction
from itertools import product
import json
from math import lcm, prod
from pathlib import Path
import time

from construction_kernel import verify, uncovered_progression, prime
from owens_branch_search import primes_through
from structural_cover_search import (Diagram, DiagramLimit, ONE, ZERO,
    coordinates, condition_from_row, row_from_condition)
from verify_construction_searches import verify_uncovered_progression

ROOT = Path(__file__).resolve().parents[1]


def divisors(exponents, ps):
    return sorted(prod(p**e for p,e in zip(ps,es))
                  for es in product(*(range(e+1) for e in exponents)))


def factor_over(n, ps):
    out = []
    for p in ps:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        out.append(e)
    assert n == 1
    return out


def crt_prefix(es, rs, ps):
    a,m = 0,1
    for p,e,r in zip(ps,es,rs):
        n = p**e
        if n > 1:
            a += m*((r-a)*pow(m,-1,n) % n)
            m *= n
    return a,m


def shapes(q, minimum, ps, ds):
    """All 2,3,5,7 prefix shapes in the original kernel, no random sampling."""
    found = {}
    ranges = [range(E) if p <= 7 else range(1) for p,E in zip(ps,ds)]
    for es in product(*ranges):
        good = [n for n in divisors(es,ps) if n*q >= minimum]
        if len(good) < q-1:
            continue
        chosen = good[:q-1]
        m = lcm(*chosen)
        found.setdefault(m,(factor_over(m,ps),chosen))
    return [(m,*found[m]) for m in sorted(found)]


def capacity_audit(base):
    """Exact upper bound for one fixed-divisor block per fresh prime.

    This screen allows ALL original base primes, not only 2,3,5,7.
    For q>=minimum every divisor is admissible, so a suffix minimum indexed
    by divisor count finds the smallest possible cylinder modulus.
    """
    ps,ds,minimum = base['primes'],base['depths'],base['minimum']
    candidates = [(prod(p**e for p,e in zip(ps,es)),prod(e+1 for e in es),es)
                  for es in product(*(range(E) for E in ds))]
    max_count = prod(ds)
    by_count = [None]*(max_count+2)
    for m,t,_ in candidates:
        if by_count[t] is None or m < by_count[t]:
            by_count[t] = m
    for t in range(max_count,0,-1):
        if by_count[t+1] is not None:
            by_count[t] = min(x for x in [by_count[t],by_count[t+1]] if x is not None)
    capacities = []
    for q in primes_through(max_count+1):
        if q <= max(ps):
            continue
        if q >= minimum:
            m = by_count[q-1]
        else:
            threshold = (minimum+q-1)//q
            admissible = [m for m,t,es in candidates
                          if t-sum(m % n == 0 for n in range(1,threshold)) >= q-1]
            m = min(admissible) if admissible else None
        if m is not None:
            capacities.append({'root_prime':q,'smallest_cylinder_modulus':m})
    capacity = sum((Fraction(1,r['smallest_cylinder_modulus']) for r in capacities),Fraction())
    residual = Fraction(base['exact_uncovered_states'],base['kernel_states'])
    return {'family':'one independent fixed-divisor block per fresh prime; all original prefix primes allowed',
            'prefix_shapes_checked':len(candidates), 'maximum_divisor_count':max_count,
            'roots_checked':len(capacities),'capacity_upper_bound':str(capacity),
            'base_residual_fraction':str(residual),
            'unavoidable_projected_residual_lower_bound':str(max(Fraction(),residual-capacity)),
            'status':'RESTRICTED_REPAIR_EXCLUDED' if capacity < residual else 'NOT_EXCLUDED',
            'root_capacities':capacities}


def check_patch_schema(result, base):
    """Check saved rows against divisor-block identities using integer arithmetic."""
    ps,ds = base['primes'],base['depths']
    qs = result['primes'][len(ps):]
    assert result['primes'][:len(ps)] == ps
    assert qs == sorted(set(qs)) and all(prime(q) and q > max(ps) for q in qs)
    assert result['depths'] == ds+[1]*len(qs)
    assert result['minimum'] == base['minimum']
    assert len({p['root_prime'] for p in result['patches']}) == len(result['patches'])
    k = 0
    for row in base['rows']:
        assert result['rows'][k] == {'exponents':row['exponents']+[0]*len(qs),
                                     'residues':row['residues']+[0]*len(qs)}
        k += 1
    for patch in result['patches']:
        q,a,m,ns = (patch[t] for t in ['root_prime','residue','modulus','divisors'])
        assert q in qs and len(ns) == q-1 and len(set(ns)) == q-1
        assert all(type(n) is int and n >= 1 and n*q >= base['minimum'] for n in ns)
        assert m == lcm(*ns) and 0 <= a < m
        for digit,n in enumerate(ns,1):
            row = result['rows'][k]
            es,rs = row['exponents'][:len(ps)],row['residues'][:len(ps)]
            assert len(es) == len(ps) and all(0 <= e < E for e,E in zip(es,ds))
            assert prod(p**e for p,e in zip(ps,es)) == n
            assert rs == [a % p**e for p,e in zip(ps,es)]
            assert row['exponents'][len(ps):] == [int(t == q) for t in qs]
            assert row['residues'][len(ps):] == [digit-1 if t == q else 0 for t in qs]
            k += 1
    assert k == len(result['rows'])
    return {'status':'DIVISOR_BLOCK_SCHEMAS_VERIFIED','blocks':len(result['patches']),
            'base_rows':len(base['rows']),'patch_rows':k-len(base['rows'])}


def attach_negative_witness(result, base, point):
    """Choose a missed divisor at each independent new prime, then check CRT.

    The independent modular checker checks every emitted row, including its
    finite completion. No decision-diagram count is trusted by this check.
    """
    result['schema_check'] = check_patch_schema(result,base)
    ps,ds = base['primes'],base['depths']
    original = uncovered_progression(point,ps,ds)['residue']
    qs = result['primes'][len(ps):]
    extra = [0]*len(qs)
    for patch in result['patches']:
        a = patch['residue']
        missed = [j for j,n in enumerate(patch['divisors']) if original % n != a % n]
        assert missed, 'Base point lies in a fully patched cylinder'
        extra[qs.index(patch['root_prime'])] = missed[0]
    full_point = list(point)+extra
    result['verification'] = {
        'status':'UNCOVERED', 'point':full_point,
        'uncovered_progression':uncovered_progression(full_point,result['primes'],result['depths']),
        'method':'projected_witness_lift_then_independent_modular_check'}
    result['independent_witness_check'] = verify_uncovered_progression(result)
    result['status'] = 'UNCOVERED'


def search(base, maximum_root=113, seconds=120, node_limit=200000):
    ps,ds = base['primes'],base['depths']
    qs = [q for q in primes_through(maximum_root) if q > max(ps)]
    start = time.monotonic()
    domains,blocks = coordinates(ps,ds,base.get('layout','prime-first'))
    diagram = Diagram(domains,node_limit,start+seconds)
    root = ONE
    for row in base['rows']:
        root = diagram.erase(root,condition_from_row(row,ps,ds,blocks))
        if len(diagram.nodes) > node_limit//4:
            old = diagram
            diagram,root = old.compact(root)
            old.count.cache_clear()
    before = diagram.count(root)
    assert before == base['exact_uncovered_states']
    total = prod(domains)
    result = {'minimum':base['minimum'],'primes':ps+qs,'depths':ds+[1]*len(qs),
              'method':'independent_fresh_prime_divisor_blocks',
              'base_rows':len(base['rows']), 'base_kernel_states':total,
              'base_uncovered_fraction':str(Fraction(before,total)),
              'patches':[], 'rows':[
                  {'exponents':r['exponents']+[0]*len(qs),
                   'residues':r['residues']+[0]*len(qs)} for r in base['rows']],
              'node_limit':node_limit, 'seconds_budget':seconds,
              'shape_policy':'exact best placement among first four distinct admissible moduli'}
    outcome = 'ROOTS_EXHAUSTED'
    try:
        for j,q in enumerate(qs):
            if time.monotonic()-start >= seconds:
                raise DiagramLimit('time limit')
            best = None
            for m,es,chosen in shapes(q,base['minimum'],ps,ds)[:4]:
                selected = [v for block,e in zip(blocks,es) for v in block[:e]]
                gain,assignment,_ = diagram.maximum_cofactor(root,selected)
                if best is None or gain > best[0]:
                    best = gain,assignment,m,es,chosen
            if best is None or best[0] == 0:
                continue
            gain,assignment,m,es,chosen = best
            mask = row_from_condition(es,assignment,ps,ds,blocks)
            a,check_m = crt_prefix(es,mask['residues'],ps)
            assert m == check_m == lcm(*chosen)
            new_root = diagram.erase(root,assignment)
            after = diagram.count(new_root)
            assert before-after == gain
            patch = {'root_prime':q, 'residue':a, 'modulus':m,
                     'divisors':chosen,'newly_removed_base_states':gain,
                     'remaining_base_fraction':str(Fraction(after,total))}
            result['patches'].append(patch)
            for digit,n in enumerate(chosen,1):
                ex = factor_over(n,ps)+[0]*len(qs)
                rs = [a % p**e for p,e in zip(ps,ex)]+[0]*len(qs)
                ex[len(ps)+j] = 1
                rs[len(ps)+j] = digit-1
                result['rows'].append({'exponents':ex,'residues':rs})
            root,before = new_root,after
            print(json.dumps({k:v for k,v in patch.items() if k != 'divisors'}),flush=True)
            if root == ZERO:
                outcome = 'EXACT_PROJECTED_COVER_FOUND'
                break
            old = diagram
            diagram,root = old.compact(root)
            old.count.cache_clear()
    except DiagramLimit as exc:
        outcome = str(exc).upper().replace(' ','_')
    result['search_outcome'] = outcome
    result['base_states_with_uncovered_extension'] = before
    result['projected_uncovered_fraction'] = str(Fraction(before,total))
    if root == ZERO:
        result['verification'] = verify(result,node_limit=300000)
        result['status'] = result['verification']['status']
        assert result['status'] == 'VERIFIED'
    else:
        assignment = {v:0 for v in range(len(domains))}
        cursor = root
        while cursor >= 0:
            v,cs = diagram.nodes[cursor]
            digit = next(j for j,c in enumerate(cs) if diagram.count(c,v+1))
            assignment[v] = digit
            cursor = cs[digit]
        assert cursor == ONE
        point = row_from_condition(ds,assignment,ps,ds,blocks)['residues']
        attach_negative_witness(result,base,point)
    result['elapsed_seconds'] = round(time.monotonic()-start,3)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base')
    parser.add_argument('--maximum-root',type=int,default=113)
    parser.add_argument('--seconds',type=int,default=120)
    parser.add_argument('--nodes',type=int,default=200000)
    parser.add_argument('--capacity-audit',action='store_true')
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    base = json.loads((ROOT/args.base).read_text())
    result = (capacity_audit(base) if args.capacity_audit else
              search(base,args.maximum_root,args.seconds,args.nodes))
    result['base_file'] = args.base
    (ROOT/args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','patches','root_capacities')},indent=2))
