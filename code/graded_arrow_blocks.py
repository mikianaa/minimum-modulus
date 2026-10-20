"""Resource-disjoint arrow supplies and exact joint progression blocks.

Only a complete union of guaranteed progressions certifies this method.
A projected gap does not disprove incidental coverage by expanded rows.
"""
import argparse
from fractions import Fraction
from itertools import product
import json
from math import prod
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tmp'/'construction-deps'))
import numpy as np
from construction_kernel import prime
from owens_branch_search import primes_through


def supply(order, fixed, with_nodes=False):
    """Grades and optional DAG: fixed p-powers and disjoint tail groups."""
    nodes = [{'kind':'unit','grade':1}] if with_nodes else None
    current = [(1,0)]
    for p in order:
        e = fixed.get(p,0)
        current.sort(key=lambda t:(-t[0],t[1]))
        following = current[:]
        for j in range(1,e+1):
            for grade,child in current:
                idx = len(nodes) if with_nodes else len(following)
                if with_nodes:
                    nodes.append(dict(kind='fixed',prime=p,exponent=j,
                                      children=[child],grade=grade*p**j))
                following.append((grade*p**j,idx))
        for start in range(0,len(current)-p+2,p-1):
            group = current[start:start+p-1]
            grade = min(g for g,_ in group)*p**(e+1)
            idx = len(nodes) if with_nodes else len(following)
            if with_nodes:
                nodes.append(dict(kind='tail',prime=p,exponent=e+1,
                                  children=[n for _,n in group],grade=grade))
            following.append((grade,idx))
        current = following
    current.sort(key=lambda t:(-t[0],t[1]))
    return current,nodes


def stage_shapes(q, exponent, tail, minimum, core_primes, core_exponents):
    earlier = primes_through(q-1)
    eligible = [(p,e) for p,e in zip(core_primes,core_exponents) if p<q]
    factor = q**(exponent-int(tail))
    threshold = (minimum+q**exponent-1)//q**exponent
    out = []
    for es in product(*(range(e+1) for _,e in eligible)):
        fixed = {p:e for (p,_),e in zip(eligible,es) if e}
        m = prod(p**e for p,e in fixed.items())
        special = [p for p in earlier if p in fixed]
        ordinary = [p for p in earlier if p not in fixed]
        orders = [earlier,special+ordinary,list(reversed(special))+ordinary]
        best = None
        seen = set()
        for order in orders:
            if tuple(order) in seen:
                continue
            seen.add(tuple(order))
            values,_ = supply(order,fixed)
            good = sum(g>=threshold for g,_ in values)
            count = good//(q-1) if tail else good
            if best is None or count>best['available_blocks']:
                best = dict(root_prime=q,exponent=exponent,tail=tail,
                            order=order,fixed_exponents={str(p):e for p,e in fixed.items()},
                            input_threshold=threshold,available_inputs=good,
                            available_blocks=count,modulus=m*factor)
        if best['available_blocks']:
            out.append(best)
    return out


def best_block(uncovered, shapes):
    best = None
    period = len(uncovered)
    for shape in shapes:
        m = shape['modulus']
        assert period%m==0
        counts = uncovered.reshape(period//m,m).sum(axis=0,dtype=np.int64)
        k = min(m,shape['available_blocks'])
        indices = np.argsort(-counts,kind='stable')[:k]
        residues = [int(a) for a in indices if counts[a]>0]
        gain = int(counts[indices].sum())
        if best is None or gain>best[0]:
            best = gain,dict(shape,residues=residues)
    return best


def add_block(counts, block, change):
    if block is not None:
        for a in block['residues']:
            counts[a::block['modulus']] += change


def check(result):
    """Check grade induction, disjoint stage spectra, and literal coverage."""
    minimum = result['minimum']
    assert minimum>=2
    ps,es = result['core_primes'],result['core_exponents']
    assert len(set(ps))==len(ps) and all(prime(p) for p in ps)
    assert all(type(e) is int and e>=0 for e in es)
    L = prod(p**e for p,e in zip(ps,es))
    assert L==result['period']
    used = {}
    covered = bytearray(L)
    selected_blocks = 0
    for block in result['blocks']:
        if block is None:
            continue
        q,j,tail = block['root_prime'],block['exponent'],block['tail']
        assert prime(q) and type(j) is int and j>=1 and type(tail) is bool
        previous = used.setdefault(q,[])
        for old_j,old_tail in previous:
            assert not (j==old_j or (tail and old_j>=j) or (old_tail and j>=old_j))
        previous.append((j,tail))
        fixed = {int(p):e for p,e in block['fixed_exponents'].items()}
        order = block['order']
        assert len(set(order))==len(order) and all(prime(p) and p<q for p in order)
        assert set(fixed)<=set(order) and all(type(e) is int and e>=1 for e in fixed.values())
        grades = [1]
        for p in order:
            e = fixed.get(p,0)
            ordered = sorted(grades,reverse=True)
            copies = [g*p**a for a in range(e+1) for g in grades]
            arrows = [ordered[k+p-2]*p**(e+1)
                      for k in range(0,len(ordered)-p+2,p-1)]
            grades = copies+arrows
        threshold = (minimum+q**j-1)//q**j
        good = sum(g>=threshold for g in grades)
        size = q-1 if tail else 1
        assert good==block['available_inputs']
        assert good//size==block['available_blocks']
        assert threshold==block['input_threshold']
        modulus = prod(p**e for p,e in fixed.items())*q**(j-int(tail))
        assert modulus==block['modulus'] and L%modulus==0
        residues = block['residues']
        assert len(residues)==len(set(residues)) and len(residues)*size<=good
        assert all(type(a) is int and 0<=a<modulus for a in residues)
        for a in residues:
            covered[a::modulus] = b'\x01'*len(range(a,L,modulus))
        selected_blocks += len(residues)
    missed = [i for i,v in enumerate(covered) if not v]
    assert len(missed)==result['guaranteed_uncovered']
    return dict(status='VERIFIED_CONSTRUCTION_SCHEMA' if not missed else 'INCOMPLETE_GUARANTEED_COVER',
                period=L,guaranteed_uncovered=len(missed),progressions=selected_blocks,
                first_projected_gap=missed[0] if missed else None,
                gap_scope='Gap in guaranteed progressions; not a disproof of incidental coverage by expanded rows')


def capacity_audit(minimum,largest,exponents):
    """Upper bound for precisely the shape/order menus searched here."""
    ps=[2,3,5,7]
    stages=[]
    total=Fraction()
    for q in primes_through(largest):
        depth=exponents[ps.index(q)]+1 if q in ps else 1
        for j in range(1,depth+1):
            candidates=stage_shapes(q,j,j==depth,minimum,ps,exponents)
            if not candidates:
                continue
            bound=lambda s:Fraction(min(s['modulus'],s['available_blocks']),s['modulus'])
            best=max(candidates,key=bound)
            total+=bound(best)
            stages.append(dict(root_prime=q,exponent=j,tail=j==depth,
                               best_shape=best,capacity=str(bound(best))))
    return dict(status='RESTRICTED_MENU_EXCLUDED' if total<1 else 'NOT_EXCLUDED',
                minimum=minimum,largest_prime=largest,core_primes=ps,
                core_exponents=exponents,total_capacity=str(total),
                stages=stages,scope='Only the recorded three order policies, one shape per largest-prime/exponent stage, and guaranteed progressions. Does not bound incidental arrow coverage or all possible constructions.')


def search(minimum,largest,exponents,sweeps=2):
    core_primes=[2,3,5,7]
    L=prod(p**e for p,e in zip(core_primes,exponents))
    stages=[]
    shapes=[]
    start=time.monotonic()
    for q in primes_through(largest):
        depth=exponents[core_primes.index(q)]+1 if q in core_primes else 1
        for j in range(1,depth+1):
            stage=(q,j,j==depth)
            candidates=stage_shapes(*stage,minimum,core_primes,exponents)
            if candidates:
                stages.append(stage)
                shapes.append(candidates)
    counts=np.zeros(L,dtype=np.int16)
    blocks=[None]*len(stages)
    trace=[]
    for sweep in range(sweeps):
        improved=0
        for i,stage in enumerate(stages):
            before=int(np.count_nonzero(counts==0))
            add_block(counts,blocks[i],-1)
            baseline=int(np.count_nonzero(counts==0))
            gain,proposal=best_block(counts==0,shapes[i])
            after=baseline-gain
            assert after<=before
            if after<before:
                blocks[i]=proposal
                improved+=1
            add_block(counts,blocks[i],1)
            assert int(np.count_nonzero(counts==0))==after
            if after<before:
                info=dict(sweep=sweep+1,root=stage[0],exponent=stage[1],tail=stage[2],
                          before=before,after=after,progression_modulus=proposal['modulus'],
                          progressions=len(proposal['residues']))
                trace.append(info)
                print(json.dumps(info),flush=True)
            if after==0:
                break
        if not improved or not np.any(counts==0):
            break
    result=dict(minimum=minimum,largest_prime=largest,core_primes=core_primes,
                core_exponents=exponents,period=L,method='graded_resource_disjoint_arrow_blocks',
                blocks=blocks,trace=trace,guaranteed_uncovered=int(np.count_nonzero(counts==0)),
                elapsed_seconds=time.monotonic()-start,
                recipe='Reconstruct supply(order,fixed,with_nodes=True), retain grades >= input_threshold in descending order, and allocate consecutive disjoint groups to recorded residues. A fixed root uses one input; a tail root uses q-1 in digit order 1,...,q-1. Fixed coordinates match the target residue; inner arrows have that target prefix. Finitize every tail occurrence with a distinct auxiliary prime > largest_prime, allocated by stage, residue, then preorder traversal. Gaps in guaranteed progressions are not asserted to miss expanded rows.')
    result['verification']=check(result)
    result['capacity_audit']=capacity_audit(minimum,largest,exponents)
    result['status']=result['verification']['status']
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum',type=int,default=50)
    parser.add_argument('--largest-prime',type=int,default=113)
    parser.add_argument('--core-exponents',default='5,3,2,1')
    parser.add_argument('--sweeps',type=int,default=2)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    exponents=list(map(int,args.core_exponents.split(',')))
    if len(exponents)!=4 or any(e<0 for e in exponents):
        parser.error('Specify four nonnegative core exponents')
    result=search(args.minimum,args.largest_prime,exponents,args.sweeps)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','period','guaranteed_uncovered','elapsed_seconds','verification']},indent=2))
