"""Independent full-state audits of core marginals and relative regions."""
from itertools import product
import json
from math import prod
from pathlib import Path
import random

from construction_kernel import matches,sizes,resources,labels
from core_region_arrows import Core,Catalog,search,validate,assignments


def main():
    rng=random.Random(503)
    comparisons=0
    ps,ds=[2,3,5],[3,2,1]
    dims=sizes(ps,ds)
    points=list(product(*(range(n) for n in dims)))
    exps=list(resources(ps,ds,3))
    for _ in range(12):
        rows=[dict(exponents=list(es),residues=[rng.randrange(labels(p,E,e))
              for p,E,e in zip(ps,ds,es)]) for es in rng.sample(exps,8)]
        missed=[x for x in points if not any(matches(x,r,ps,ds) for r in rows)]
        cert=dict(minimum=3,primes=ps,depths=ds,rows=rows,exact_uncovered_states=len(missed))
        core=Core(cert)
        assert all(bool(core.residual[x])==(x in missed) for x in points)
        for es in exps:
            marginal=core.marginal(core.residual,es)
            for rs in product(*(range(labels(p,E,e)) for p,E,e in zip(ps,ds,es))):
                row=dict(exponents=es,residues=rs)
                expected=[x for x in missed if matches(x,row,ps,ds)]
                assert int(marginal[rs])==len(expected)
                if expected:
                    pid=core.profile(es,rs)
                    for i,p in enumerate(ps):
                        actual=set(core.profiles[pid][i])
                        wanted={int(core.values[i][x[i]]) for x in expected}
                        assert actual==wanted
                comparisons+=1
    # A core with real gaps, followed by two noncore arrow stages.
    rows=[dict(exponents=[2,0,0],residues=[0,0,0]),
          dict(exponents=[0,1,0],residues=[0,0,0]),
          dict(exponents=[0,0,1],residues=[0,0,0]),
          dict(exponents=[1,1,0],residues=[0,1,0])]
    missed=[x for x in points if not any(matches(x,r,ps,ds) for r in rows)]
    cert=dict(minimum=3,primes=ps,depths=ds,rows=rows,exact_uncovered_states=len(missed))
    forest_mode='--forest' in __import__('sys').argv
    union_mode='--unions' in __import__('sys').argv
    result=search(cert,largest=11,max_volume=96,per_shape=4,seconds=30,
                  branch_sensitive=True,forest_supply=forest_mode,union_targets=union_mode)
    assert result['guaranteed_uncovered_states']==0
    Path('data/core-region-test-control.json').write_text(json.dumps(result,indent=2)+'\n')
    catalog=Catalog(Core(cert))
    for source in result['sources']:
        catalog.append(source)
    allps=[2,3,5,7,11]
    def representative(p,E,state):
        a,d=divmod(state,p-1)
        return a+(d+1)*p**(E-1)
    count=0
    implications=0
    def in_region(state,region):
        choices=[region.get('unions',{}).get(str(i),[r]) for i,r in enumerate(region['residues'])]
        return any(matches(state,dict(exponents=region['exponents'],residues=rs),ps,ds)
                   for rs in product(*choices))
    union_comparisons=0
    for es in product(*(range(E+1) for E in ds)):
        marginal=core.marginal(core.residual,es)
        for gain,rs,unions in assignments(core,marginal,es,2,True):
            if unions:
                region=dict(exponents=es,residues=rs,unions=unions)
                assert gain==sum(in_region(x,region) for x in points if core.residual[x])
                union_comparisons+=1
    for state in product(*(range(d) for d in dims),range(6),range(10)):
        xs={p:representative(p,E,s) for p,E,s in zip(allps,ds+[1,1],state)}
        def own(s,idx):
            node=result['sources'][s]['nodes'][idx]
            if node['kind']=='unit':
                return True
            if node['kind']=='forest':
                p=node['prime']
                def branch_hit(b):
                    e,a=b['exponent'],b['prefix']
                    hit=((xs[p]-a)%p**e==0 if b['kind']=='fixed' else
                         (xs[p]-a-b['digit']*p**(e-1))%p**e==0)
                    return hit and own(s,b['child'])
                return any(branch_hit(b) for b in node['branches'])
            p,e,a=node['prime'],node['exponent'],node['prefix']
            if node['kind']=='fixed':
                return xs[p]%p**e==a and own(s,node['children'][0])
            if xs[p]%p**(e-1)!=a:
                return False
            y=(xs[p]-a)//p**(e-1)
            while y%p==0:
                y//=p
            d=y%p
            return d in node['digits'] and own(s,node['children'][node['digits'].index(d)])
        covered=any(matches(state[:3],r,ps,ds) for r in rows)
        for s,source in enumerate(result['sources']):
            q=source['root_prime']
            digit=xs[q]%q
            hit=own(s,source['root_inputs'][digit-1])
            if in_region(state[:3],source['region']):
                assert covered or hit
                implications+=1
            covered=covered or hit
        count+=1
    # A false root minimum and a false omitted-core input must be rejected.
    bad=json.loads(json.dumps(result))
    bad['minimum']=50
    try:
        validate(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('Inconsistent minimum accepted')
    bad=json.loads(json.dumps(result))
    source=bad['sources'][0]
    key=next(iter(source['configs']))
    source['configs'][key]['digits']=[]
    try:
        validate(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('False core omission accepted')
    if forest_mode:
        bad=json.loads(json.dumps(result))
        node=next(n for s in bad['sources'] for n in s['nodes'] if n['kind']=='forest')
        node['branches'].append(node['branches'][0].copy())
        try:
            validate(bad)
        except AssertionError:
            pass
        else:
            raise AssertionError('Duplicated exponent resource accepted')
    report=dict(status='PASS',exact_marginal_comparisons=comparisons,
                forest_supply=forest_mode,branch_sensitive=True,
                union_targets=union_mode,union_target_comparisons=union_comparisons,
                expanded_relative_states=count,relative_coverage_implications=implications,
                false_minimum_rejected=True,false_core_omission_rejected=True)
    name='core-region-union-tests' if union_mode else ('core-region-forest-tests' if forest_mode else 'core-region-arrow-tests')
    Path('data/'+name+'.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
