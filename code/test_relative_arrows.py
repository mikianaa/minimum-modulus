"""Exhaustive state tests of relative-arrow borrowing and own moduli."""
from itertools import product
import json
from math import prod
from pathlib import Path

from relative_arrow_search import Catalog,make_supply,validate


def main():
    catalog=Catalog()
    for q,a,m in [(3,0,2),(5,1,4),(7,3,4)]:
        fixed={2:m.bit_length()-1}
        borrowed=catalog.borrowed(a,m)
        order=[p for p in [2,3,5] if p<q]
        values,nodes=make_supply(order,fixed,borrowed)
        selected=[i for g,i in sorted(values)[:q-1]]
        assert len(selected)==q-1
        catalog.append(dict(root_prime=q,residue=a,modulus=m,
                            nodes=nodes,root_inputs=selected,
                            order=order,borrowed={str(p):ds for p,ds in borrowed.items()}))
    independent=Catalog(catalog.sources,independent=True)
    ps=[2,3,5,7]
    depths=[4,1,1,1]
    coordinates=[[a+d*p**(E-1) for a in range(p**(E-1)) for d in range(1,p)]
                 for p,E in zip(ps,depths)]
    points=[dict(zip(ps,x)) for x in product(*coordinates)]

    def node_value(s,idx,x,allow_borrow):
        source=catalog.sources[s]
        n=source['nodes'][idx]
        if n['kind']=='unit':
            return True
        p,e=n['prime'],n['exponent']
        a=source['residue']
        if n['kind']=='fixed':
            return x[p]%p**e==a%p**e and node_value(s,n['children'][0],x,allow_borrow)
        prefix=a%p**(e-1)
        if x[p]%p**(e-1)!=prefix:
            return False
        z=(x[p]-prefix)//p**(e-1)
        assert z>0
        while z%p==0:
            z//=p
        d=z%p
        if d in n['digits']:
            return node_value(s,n['children'][n['digits'].index(d)],x,allow_borrow)
        if not allow_borrow:
            return False
        old=catalog.by_prime[p]
        return node_value(old,catalog.sources[old]['root_inputs'][d-1],x,True)

    def own_source(s,x):
        source=catalog.sources[s]
        p=source['root_prime']
        z=x[p]
        while z%p==0:
            z//=p
        return node_value(s,source['root_inputs'][z%p-1],x,False)

    covers=0
    dependency_implications=0
    for x in points:
        prior=False
        for s,source in enumerate(catalog.sources):
            owned=own_source(s,x)
            if x[2]%source['modulus']==source['residue']:
                assert prior or owned
                covers+=1
            for idx in range(len(source['nodes'])):
                if node_value(s,idx,x,True):
                    assert node_value(s,idx,x,False) or prior
                    dependency_implications+=1
            prior=prior or owned
        assert prior
    tests=0
    for s,source in enumerate(catalog.sources):
        for idx in range(len(source['nodes'])):
            for m in [1,2,4,8]:
                for a in range(m):
                    fast=catalog.contains(s,idx,a,m)
                    other=independent.contains(s,idx,a,m)
                    assert fast==other
                    if fast:
                        assert all(node_value(s,idx,x,True) for x in points if x[2]%m==a)
                    tests+=1

    def leaves(s,idx):
        source=catalog.sources[s]
        n=source['nodes'][idx]
        if n['kind']=='unit':
            return [{}]
        result=[]
        for child in n['children']:
            for row in leaves(s,child):
                assert n['prime'] not in row
                row[n['prime']]=(n['exponent'],n['kind']=='tail')
                result.append(row)
        return result
    resources=[]
    for s,source in enumerate(catalog.sources):
        p=source['root_prime']
        for child in source['root_inputs']:
            for row in leaves(s,child):
                row[p]=(1,True)
                assert prod(q**e for q,(e,_) in row.items())>=3
                resources.append(row)
    def intersect(a,b):
        for p in ps:
            ea,ta=a.get(p,(0,False))
            eb,tb=b.get(p,(0,False))
            if not ta and not tb and ea!=eb:
                return False
            if ta and not tb and eb<ea or tb and not ta and ea<eb:
                return False
        return True
    for i,a in enumerate(resources):
        assert not any(intersect(a,b) for b in resources[i+1:])
    saved=json.loads(Path('data/construction-50-relative-113.json').read_text())
    assert validate(saved)['status']=='VERIFIED_PARTIAL_RELATIVE_SCHEMA'
    result=dict(status='PASS',complete_small_kernel_states=len(points),
                containment_comparisons=tests,relative_cover_instances=covers,
                borrowed_implications=dependency_implications,
                pairwise_disjoint_exponent_families=len(resources),
                saved_minimum_50_schema='VERIFIED_PARTIAL_RELATIVE_SCHEMA')
    Path('data/relative-arrow-tests.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
