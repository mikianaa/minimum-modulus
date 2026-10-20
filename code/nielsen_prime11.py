"""Explicit interpretation of Nielsen's ten prime-11 inputs.

This is a transcription candidate until its exhaustive relative audit passes.
Owens's stated permutation moves the first prime-5 input to digit 4.
Tail centers are allowed above their lowest prefix: e.g. the ternary
27-tail centered at 21 expands into 12 mod 27 and a tail at prefix 21.
"""
from itertools import product
from math import prod
import json
from pathlib import Path

from graded_arrow_blocks import np
from construction_kernel import sizes,verify
from owens_seed_core import seed
from joint_root_repair import validate,target_points

PS=[2,3,5,7]
DS=[6,4,3,2]


def fixed(p,e,r):
    return [{p:('fixed',e,r,0)}]


def tail(p,e,a,d):
    return [{p:('tail',e,a,d)}]


def multiply(a,b):
    result=[]
    for x,y in product(a,b):
        assert not set(x)&set(y)
        result.append(x|y)
    return result


def arrow(p,e,a,inputs,permutation=None):
    result=[]
    for digit,body in enumerate(inputs,1):
        d=permutation.get(digit,digit) if permutation else digit
        result+=multiply(tail(p,e,a,d),body)
    return result


def expand(family,depths=DS):
    options=[]
    for p,E in zip(PS,depths):
        if p not in family:
            options.append([(0,0)])
            continue
        kind,e,a,d=family[p]
        if kind=='fixed':
            assert e<E and 0<=a<p**e
            options.append([(e,a)])
        else:
            assert e<=E and 0<=a<p**(E-1) and 1<=d<p
            options.append([(k,(a+d*p**(k-1))%p**k) for k in range(e,E)]
                           +[(E,a*(p-1)+d-1)])
    return [dict(exponents=[e for e,r in row],residues=[r for e,r in row]) for row in product(*options)]


def inputs(binary_residue=1,permutation=True,depths=DS,raw=False):
    one=[{}]
    two=fixed(2,1,1)
    four=fixed(2,2,binary_residue)
    eight=tail(2,3,binary_residue,1)
    three=fixed(3,1,1)
    nine=fixed(3,2,3)
    twentyseven=fixed(3,3,3)
    t27=tail(3,3,21,2)
    times=multiply
    five_map=permutation if isinstance(permutation,dict) else ({1:4,4:1} if permutation else None)
    fives=lambda xs:arrow(5,1,0,xs,five_map)
    sevens=lambda xs:arrow(7,1,0,xs)
    ternary=lambda xs:sum((times(fixed(3,2,1+3*d),body) for d,body in enumerate(xs)),[])
    a1=four
    a2=eight
    a3=times(three,two)+twentyseven+times(t27,two)
    a4=times(three,four)+times(twentyseven,four)+times(t27,eight)
    a5=times(three,eight)+times(nine,eight)
    a6=fives([one,two,three,four])
    a7=fives([eight,times(three,two)+times(nine,two),times(three,four),a5])
    a8=ternary([one,two,four])+arrow(3,4,3,[one,four])+fives([
        t27,times(t27,two),times(t27,four),times(t27,eight)])
    a9=sevens([one,two,three,fives([one,two,[],four]),four,eight])
    part1=fives([ternary([one,four,[]]),[],[],[]])
    part2=sevens([a3,a4,times(three,eight),
        fives([ternary([[],[],one])+times(nine,two),eight,[],three+times(nine,four)]),
        a8,
        fives([ternary([[],[],eight]),times(three,two),times(three,four),times(three,eight)])
        +times(nine,eight)])
    families=[a1,a2,a3,a4,a5,a6,a7,a8,a9,part1+part2]
    return families if raw else [[row for family in item for row in expand(family,depths)] for item in families]


def seed_certificate(minimum):
    result=seed(minimum,True)
    dims=sizes(PS,DS)
    states=np.indices(dims,sparse=True)
    residual=np.ones(dims,dtype=bool)
    for row in result['rows']:
        hit=True
        for s,p,E,e,r in zip(states,PS,DS,row['exponents'],row['residues']):
            hit=hit & (s==r if e==E else (s//(p-1))%p**e==r)
        residual&=~hit
    result['exact_uncovered_states']=int(residual.sum())
    verify(result)
    return result


def audit(minimum=42):
    core=seed_certificate(minimum)
    rows=inputs()
    result=dict(core=core,region=dict(exponents=[2,2,0,0],residues=[1,0,0,0],
                unions={'1':[1,3,4,7]}),levels=DS,root_prime=11,tail_mode=True,inputs=rows,
                method='literal_nielsen_prime11_with_owens_five_permutation')
    _,points,needed=target_points(core,result['region'],DS,True)
    seen={}
    duplicates=[]
    gaps=[]
    below=[]
    for d,classes in enumerate(rows,1):
        covered=np.zeros(len(points),dtype=bool)
        for row in classes:
            if 11*prod(p**e for p,e in zip(PS,row['exponents']))<minimum:
                below.append(dict(input=d,row=row))
            key=tuple(row['exponents'])
            if key in seen:
                duplicates.append([seen[key],d,list(key)])
            seen[key]=d
            hit=np.ones(len(points),dtype=bool)
            for i,(p,E,e,r) in enumerate(zip(PS,DS,row['exponents'],row['residues'])):
                hit&=(points[:,i]==r if e==E else points[:,i]//(p-1)%p**e==r)
            covered|=hit
        missed=np.flatnonzero(~covered)
        gaps.append(dict(input=d,rows=len(classes),missing=len(missed),
                         first_missing=points[int(missed[0])].tolist() if len(missed) else None))
    result['audit']=dict(resource_duplicates=duplicates,input_gaps=gaps,target_states=needed,
                         below_minimum_rows=below)
    result['status']='UNVERIFIED_TRANSCRIPTION'
    if not duplicates and not below and not any(g['missing'] for g in gaps):
        result['verification']=validate(result)
        result['status']=result['verification']['status']
    return result


if __name__=='__main__':
    result=audit()
    Path('data/nielsen-prime11-transcription.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],audit=result['audit']),indent=2))
