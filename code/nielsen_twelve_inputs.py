"""Twelve disjoint elementary inputs on 1 mod 4 intersect 2 mod 3.

These are the first twelve inputs of Nielsen's section 4.8, with centers
chosen explicitly. They use only core primes, so the additional two
coordinates below are unconstrained.
"""
from nielsen_prime11 import fixed,tail,multiply,arrow,expand


def inputs(levels=(3,4,2,1)):
    one=[{}]
    two=fixed(2,1,1)
    four=fixed(2,2,1)
    eight=tail(2,3,1,1)
    three=fixed(3,1,2)
    binary=[one,two,four,eight]
    threes=[multiply(three,b) for b in binary]
    nines=[arrow(3,2,2,binary[:2]),arrow(3,2,2,binary[2:])]
    fives=[arrow(5,1,0,binary),arrow(5,1,0,threes)]
    first=[four,eight]+threes+nines+fives
    all_twelve=[one,two]+first
    families=first+[arrow(7,1,0,all_twelve[:6]),arrow(7,1,0,all_twelve[6:])]
    result=[]
    seen=set()
    for group in families:
        rows=[]
        for family in group:
            for row in expand(family,levels):
                row=dict(exponents=row['exponents']+[0,0],residues=row['residues']+[0,0])
                assert tuple(row['exponents']) not in seen
                seen.add(tuple(row['exponents']))
                rows.append(row)
        result.append(rows)
    return result
