"""Relative arrow construction on exact core-kernel regions, not samples.

Every target can fix a prime-tail digit. Missing core inputs are justified
against the entire saved core residual. Other missing inputs refer to
earlier largest-prime stages. A finite complete cover requires zero gap.
"""
import argparse
from functools import lru_cache
from fractions import Fraction
from itertools import product
import json
from math import prod
from pathlib import Path
import time

from graded_arrow_blocks import np
from construction_kernel import sizes,labels,verify,prime
from owens_branch_search import primes_through
from prime_forest_supply import pack


class Core:
    def __init__(self,certificate):
        self.certificate=certificate
        self.ps=certificate['primes']
        self.ds=certificate['depths']
        self.dims=sizes(self.ps,self.ds)
        self.total=prod(self.dims)
        self.values=[np.array([a+d*p**(E-1) for a in range(p**(E-1))
                              for d in range(1,p)],dtype=np.int64)
                     for p,E in zip(self.ps,self.ds)]
        self.residual=np.ones(self.dims,dtype=bool)
        for row in certificate['rows']:
            self.residual[self.index(row['exponents'],row['residues'])]=False
        assert int(self.residual.sum())==certificate['exact_uncovered_states']
        self.profile_ids={}
        self.profiles=[]
        self.contexts=[]

    def index(self,es,rs,unions=None):
        indices=[]
        unions=unions or {}
        for i,(p,E,e,r,size) in enumerate(zip(self.ps,self.ds,es,rs,self.dims)):
            choices=unions.get(str(i),[r])
            if e==E:
                indices.append(np.array(choices))
            else:
                indices.append(np.flatnonzero(np.isin((np.arange(size)//(p-1))%p**e,choices)))
        return np.ix_(*indices)

    def region_index(self,region):
        if 'pieces' in region:
            mask=np.zeros(self.dims,dtype=bool)
            for piece in region['pieces']:
                mask[self.region_index(piece)]=True
            return mask
        return self.index(region['exponents'],region['residues'],region.get('unions'))

    def marginal(self,residual,es):
        shape=[]
        summed=[]
        for p,E,e,size in zip(self.ps,self.ds,es,self.dims):
            if e==E:
                shape.append(size)
            else:
                k=len(shape)
                shape.extend([p**(E-1-e),p**e,p-1])
                summed.extend([k,k+2])
        x=residual.reshape(shape)
        return x.sum(axis=tuple(summed),dtype=np.int64) if summed else x.astype(np.int64)

    def profile(self,es,rs,residual=None,epoch=0,unions=None):
        key=epoch,tuple(es),tuple(rs),tuple((k,tuple(v)) for k,v in sorted((unions or {}).items()))
        if key in self.profile_ids:
            return self.profile_ids[key]
        index=self.index(es,rs,unions)
        return self._profile(index,self.residual if residual is None else residual,epoch,key)

    def _profile(self,index,residual,epoch,key):
        if key in self.profile_ids:
            return self.profile_ids[key]
        region=residual[index]
        assert region.any()
        profile=[]
        for i in range(len(self.ps)):
            used=region.any(axis=tuple(j for j in range(len(self.ps)) if j!=i))
            indices=index[i].reshape(-1)[used]
            profile.append(tuple(map(int,self.values[i][indices])))
        idx=len(self.profiles)
        self.profiles.append(tuple(profile))
        self.contexts.append((epoch,index,residual))
        self.profile_ids[key]=idx
        return idx

    def condition(self,pid,p,values):
        epoch,index,residual=self.contexts[pid]
        i=self.ps.index(p)
        indices=[v.reshape(-1) for v in index]
        allowed=set(values)
        indices[i]=np.array([j for j in indices[i] if int(self.values[i][j]) in allowed],dtype=np.int64)
        key='condition',epoch,tuple(tuple(map(int,x)) for x in indices)
        return self._profile(np.ix_(*indices),residual,epoch,key)

    def configs(self,profile_id):
        result={}
        for p,E,values in zip(self.ps,self.ds,self.profiles[profile_id]):
            f=E-1
            while len({x%p**f for x in values})>1:
                f-=1
            prefix=values[0]%p**f
            digits=set()
            for x in values:
                y=(x-prefix)//p**f
                assert y>0
                while y%p==0:
                    y//=p
                digits.add(y%p)
            result[p]=(f,prefix,tuple(sorted(digits)))
        return result

    def forest_states(self,pid):
        return {p:(E,tuple(sorted(int(x%p**(E-1))*(p-1)+int(x//p**(E-1))-1
                                 for x in values)))
                for p,E,values in zip(self.ps,self.ds,self.profiles[pid])}


def supply(order,configs,borrowed,materialize=True,forest_states=None):
    nodes=[dict(kind='unit',grade=1)]
    current=[(1,0)]
    for p in order:
        if forest_states is not None and p in forest_states:
            E,states=forest_states[p]
            old=sorted(current,key=lambda x:(-x[0],x[1]))
            nxt=current[:]
            for forest in pack(p,E,states,len(old)):
                branches=[dict(kind=b['kind'],exponent=b['exponent'],prefix=b['prefix'],
                               digit=b['digit'],child=old[b['old_input']][1]) for b in forest]
                grade=min(p**b['exponent']*old[b['old_input']][0] for b in forest)
                idx=len(nodes) if materialize else len(nxt)
                if materialize:
                    nodes.append(dict(kind='forest',prime=p,branches=branches,grade=grade))
                nxt.append((grade,idx))
            current=nxt
            continue
        f,a,required=configs.get(p,(0,0,tuple(range(1,p))))
        free=borrowed.get(p,())
        digits=tuple(d for d in required if d not in free)
        if not digits:
            return None,None
        old=sorted(current,key=lambda x:(-x[0],x[1]))
        nxt=current[:]
        for j in range(1,f+1):
            for grade,child in old:
                idx=len(nodes) if materialize else len(nxt)
                if materialize:
                    nodes.append(dict(kind='fixed',prime=p,exponent=j,prefix=a%p**j,
                                      children=[child],grade=grade*p**j))
                nxt.append((grade*p**j,idx))
        k=len(digits)
        for start in range(0,len(old)-k+1,k):
            group=old[start:start+k]
            grade=min(g for g,_ in group)*p**(f+1)
            idx=len(nodes) if materialize else len(nxt)
            if materialize:
                nodes.append(dict(kind='tail',prime=p,exponent=f+1,prefix=a,
                                  digits=list(digits),borrowed=list(free),
                                  children=[i for _,i in group],grade=grade))
            nxt.append((grade,idx))
        current=nxt
    return sorted(current),nodes if materialize else None


class Catalog:
    def __init__(self,core,branch_sensitive=False):
        self.core=core
        self.sources=[]
        self.by_prime={}
        self.branch_sensitive=branch_sensitive

    def append(self,source):
        p=source['root_prime']
        assert p not in self.by_prime
        self.by_prime[p]=len(self.sources)
        self.sources.append(source)

    @lru_cache(maxsize=300000)
    def contains(self,s,idx,profile_id):
        node=self.sources[s]['nodes'][idx]
        if node['kind']=='relative_core_union':
            old=self.by_prime[node['reference_prime']]
            assert old<s
            prior=self.sources[old]['nodes'][node['reference_input']]
            assert prior['kind']=='row_union' and prior['levels']==node['levels'] and prior['tail_mode']==node['tail_mode']
            combined=dict(node,kind='row_union',rows=node['rows']+prior['rows'])
            return row_union_contains(self.core,profile_id,combined)
        if node['kind']=='row_union':
            return row_union_contains(self.core,profile_id,node)
        if node['kind']=='unit':
            return True
        if node['kind']=='forest':
            p=node['prime']
            i=self.core.ps.index(p)
            values=self.core.profiles[profile_id][i]
            remaining=set(values)
            for b in node['branches']:
                e,a=b['exponent'],b['prefix']
                hit=[x for x in values if (x%p**e==a if b['kind']=='fixed' else
                     x%p**(e-1)==a and x//p**(e-1)==b['digit'])]
                if hit:
                    child_pid=self.core.condition(profile_id,p,hit) if self.branch_sensitive else profile_id
                    if self.contains(s,b['child'],child_pid):
                        remaining.difference_update(hit)
            return not remaining
        p,e,a=node['prime'],node['exponent'],node['prefix']
        if p in self.core.ps:
            values=self.core.profiles[profile_id][self.core.ps.index(p)]
            if node['kind']=='fixed':
                return all(x%p**e==a for x in values) and self.contains(s,node['children'][0],profile_id)
            if any(x%p**(e-1)!=a for x in values):
                return False
            possible=set()
            digit_values={}
            for x in values:
                y=(x-a)//p**(e-1)
                assert y>0
                while y%p==0:
                    y//=p
                possible.add(y%p)
                digit_values.setdefault(y%p,[]).append(x)
        else:
            assert node['kind']=='tail' and e==1 and a==0
            possible=range(1,p)
        own=dict(zip(node['digits'],node['children']))
        for d in possible:
            if d in own:
                child_profile=(self.core.condition(profile_id,p,digit_values[d])
                               if self.branch_sensitive and p in self.core.ps else profile_id)
                if not self.contains(s,own[d],child_profile):
                    return False
            elif d in node['borrowed']:
                old=self.by_prime[p]
                assert old<s
                if not self.contains(old,self.sources[old]['root_inputs'][d-1],profile_id):
                    return False
            else:
                return False
        return True

    def borrowed(self,profile_id):
        result={}
        for s,source in enumerate(self.sources):
            ds=[d for d,i in enumerate(source['root_inputs'],1) if self.contains(s,i,profile_id)]
            if ds:
                result[source['root_prime']]=tuple(ds)
        return result


@lru_cache(maxsize=128)
def row_union_table(ps,levels,tail_mode,rows,independent):
    dims=tuple(p**(E-1)*(p-1) if tail_mode else p**E for p,E in zip(ps,levels))
    table=np.zeros(dims,dtype=bool)
    for es,rs in rows:
        selections=[]
        for p,E,e,r,n in zip(ps,levels,es,rs,dims):
            states=np.arange(n)
            if independent:
                values=(states//(p-1)+(states%(p-1)+1)*p**(E-1)) if tail_mode else states
                residue=(r//(p-1)+(r%(p-1)+1)*p**(E-1)) if tail_mode and e==E else r
                mask=(values-residue)%p**e==0
            else:
                projection=states if tail_mode and e==E else ((states//(p-1) if tail_mode else states)%p**e)
                mask=projection==r
            selections.append(np.flatnonzero(mask))
        table[np.ix_(*selections)]=True
    return table


@lru_cache(maxsize=128)
def lower_coordinate_map(p,G,E,tail_mode,independent):
    result=[]
    for s in range(p**(G-1)*(p-1)):
        x=s//(p-1)+(s%(p-1)+1)*p**(G-1)
        if not tail_mode:
            result.append(x%p**E)
            continue
        a=x%p**(E-1)
        if independent:
            digits=[d for d in range(1,p) if any((x-a-d*p**(k-1))%p**k==0 for k in range(E,G+1))]
            assert len(digits)==1
            digit=digits[0]
        else:
            y=(x-a)//p**(E-1)
            while y%p==0:
                y//=p
            digit=y%p
        result.append(a*(p-1)+digit-1)
    return np.array(result,dtype=np.int64)


@lru_cache(maxsize=256)
def compact_residual_profile(core,pid,levels,tail_mode,independent):
    _,index,residual=core.contexts[pid]
    current=residual[index]
    for i,(p,G,E,v) in enumerate(zip(core.ps,core.ds,levels,index)):
        mapping=lower_coordinate_map(p,G,E,tail_mode,independent)[v.reshape(-1)]
        size=p**(E-1)*(p-1) if tail_mode else p**E
        current=np.stack([np.take(current,np.flatnonzero(mapping==s),axis=i).any(axis=i)
                          for s in range(size)],axis=i)
    return current


def row_union_contains(core,pid,node,independent=False,use_cache=True):
    """Check a mixed-prime union on the actual conditioned residual."""
    _,index,residual=core.contexts[pid]
    local_size=prod(p**(E-1)*(p-1) if node['tail_mode'] else p**E for p,E in zip(core.ps,node['levels']))
    if use_cache and local_size<=1_000_000:
        key=tuple((tuple(r['exponents']),tuple(r['residues'])) for r in node['rows'])
        table=row_union_table(tuple(core.ps),tuple(node['levels']),node['tail_mode'],key,independent)
        supported=[]
        for p,G,E,values in zip(core.ps,core.ds,node['levels'],core.profiles[pid]):
            states=np.array([(x%p**(G-1))*(p-1)+x//p**(G-1)-1 for x in values],dtype=np.int64)
            supported.append(np.unique(lower_coordinate_map(p,G,E,node['tail_mode'],independent)[states]))
        possibility=table[np.ix_(*supported)]
        if possibility.all():
            return True
        if not possibility.any():
            return False
        target=compact_residual_profile(core,pid,tuple(node['levels']),node['tail_mode'],independent)
        return not (target&~table).any()
    coordinates=[core.values[i][v.reshape(-1)] for i,v in enumerate(index)]
    remaining=residual[index].copy()
    for row in node['rows']:
        selections=[]
        for values,p,global_E,E,e,r in zip(coordinates,core.ps,core.ds,node['levels'],
                                          row['exponents'],row['residues']):
            if not node['tail_mode'] or e<E:
                mask=(values-r)%p**e==0
            else:
                a,d0=divmod(r,p-1)
                if independent:
                    mask=np.zeros(len(values),dtype=bool)
                    for k in range(E,global_E+1):
                        mask|=(values-a-(d0+1)*p**(k-1))%p**k==0
                else:
                    def hit(x):
                        if int(x)%p**(E-1)!=a:
                            return False
                        y=(int(x)-a)//p**(E-1)
                        while y%p==0:
                            y//=p
                        return y%p==d0+1
                    mask=np.array([hit(x) for x in values])
            selections.append(np.flatnonzero(mask))
        if any(len(s)==0 for s in selections):
            continue
        remaining[np.ix_(*selections)]=False
        if not remaining.any():
            return True
    return False


def joint_source(certificate):
    ps=certificate['core']['primes']
    nodes=[dict(kind='row_union',levels=certificate['levels'],tail_mode=certificate.get('tail_mode',False),
                rows=rows,grade=min(prod(p**e for p,e in zip(ps,r['exponents'])) for r in rows))
           for rows in certificate['inputs']]
    return dict(root_prime=certificate['root_prime'],region=certificate['region'],nodes=nodes,
                root_inputs=list(range(len(nodes))),joint_certificate=certificate)


def partial_joint_source(certificate,digits):
    q=certificate['root_prime']
    assert len(digits)==len(certificate['inputs'])==len(set(digits))
    assert all(type(d) is int and 1<=d<q for d in digits)
    by_digit=dict(zip(digits,certificate['inputs']))
    ps=certificate['core']['primes']
    nodes=[]
    for d in range(1,q):
        rows=by_digit.get(d,[])
        grade=min((prod(p**e for p,e in zip(ps,r['exponents'])) for r in rows),default=None)
        nodes.append(dict(kind='row_union',levels=certificate['levels'],tail_mode=certificate.get('tail_mode',False),
                          rows=rows,grade=grade))
    return dict(root_prime=q,nodes=nodes,root_inputs=list(range(q-1)),
                partial_joint_certificate=certificate,assigned_digits=list(digits))


def best_order(q,configs,borrowed,minimum,forest_states=None):
    ps=primes_through(q-1)
    special=[p for p in ps if p in configs]
    ordinary=sorted((p for p in ps if p not in special),key=lambda p:(p-1-len(borrowed.get(p,())),p))
    best=None
    factors={}
    for p in ps:
        f,_,required=configs.get(p,(0,0,tuple(range(1,p))))
        k=len([d for d in required if d not in borrowed.get(p,())])
        if k==0:
            return None,ps
        factors[p]=Fraction(f+1)+Fraction(1,k)
    interleaved=sorted(ps,key=lambda p:(-factors[p],p))
    for order in [ps,special+ordinary,list(reversed(special))+ordinary,interleaved]:
        values,_=supply(order,configs,borrowed,False,forest_states)
        if values is None:
            return None,order
        good=sum(q*g>=minimum for g,_ in values)
        if best is None or good>best[0]:
            best=good,order
    return best


def validate(result):
    certificate=result['core']
    assert result['minimum']==certificate['minimum']
    verify(certificate)
    core=Core(certificate)
    residual=core.residual.copy()
    sources=[]
    by_prime={}
    dependencies=0
    branch_sensitive=result.get('branch_sensitive',False)

    @lru_cache(maxsize=300000)
    def contains(s,idx,pid):
        node=sources[s]['nodes'][idx]
        if node['kind']=='relative_core_union':
            old=by_prime[node['reference_prime']]
            assert old<s
            prior=sources[old]['nodes'][node['reference_input']]
            assert prior['kind']=='row_union' and prior['levels']==node['levels'] and prior['tail_mode']==node['tail_mode']
            combined=dict(node,kind='row_union',rows=node['rows']+prior['rows'])
            return row_union_contains(core,pid,combined,True)
        if node['kind']=='row_union':
            return row_union_contains(core,pid,node,True)
        if node['kind']=='unit':
            return True
        if node['kind']=='forest':
            p=node['prime']
            i=core.ps.index(p)
            values=core.profiles[pid][i]
            covered=set()
            for b in node['branches']:
                e,a=b['exponent'],b['prefix']
                hit=[x for x in values if ((x-a)%p**e==0 if b['kind']=='fixed' else
                     (x-a-b['digit']*p**(e-1))%p**e==0)]
                if hit:
                    child_pid=core.condition(pid,p,hit) if branch_sensitive else pid
                    if contains(s,b['child'],child_pid):
                        covered.update(hit)
            return all(x in covered for x in values)
        p,e,a=node['prime'],node['exponent'],node['prefix']
        if p in core.ps:
            i=core.ps.index(p)
            values=core.profiles[pid][i]
            if node['kind']=='fixed':
                return all((x-a)%p**e==0 for x in values) and contains(s,node['children'][0],pid)
            if any((x-a)%p**(e-1) for x in values):
                return False
            # Independent arithmetic test for every possible tail digit.
            possible=[d for d in range(1,p) if any(
                (x-a-d*p**(k-1))%p**k==0
                for x in values for k in range(e,core.ds[i]+1))]
            digit_values={d:[x for x in values if any((x-a-d*p**(k-1))%p**k==0
                          for k in range(e,core.ds[i]+1))] for d in possible}
        else:
            assert node['kind']=='tail' and e==1 and a==0
            possible=range(1,p)
        own=dict(zip(node['digits'],node['children']))
        for d in possible:
            if d in own:
                child_profile=core.condition(pid,p,digit_values[d]) if branch_sensitive and p in core.ps else pid
                if not contains(s,own[d],child_profile):
                    return False
            elif d in node['borrowed']:
                old=by_prime[p]
                assert old<s
                if not contains(old,sources[old]['root_inputs'][d-1],pid):
                    return False
            else:
                return False
        return True

    epoch=0
    conditioning=result.get('profile_conditioning','core')
    assert conditioning in ('core','current')

    def check_region(region):
        if 'pieces' in region:
            assert region['pieces']
            # Validate every rectangular component before building the union.
            for piece in region['pieces']:
                es,rs=piece['exponents'],piece['residues']
                assert len(es)==len(rs)==len(core.ps)
                assert all(type(e) is int and 0<=e<=E for e,E in zip(es,core.ds))
                assert all(type(r) is int and 0<=r<labels(p,E,e) for p,E,e,r in zip(core.ps,core.ds,es,rs))
                assert 'unions' not in piece and 'pieces' not in piece
            selected=(residual if conditioning=='current' else core.residual)&core.region_index(region)
            index=core.index([0]*len(core.ps),[0]*len(core.ps))
            return core._profile(index,selected,epoch,('pieces',epoch,json.dumps(region,sort_keys=True)))
        es,rs=region['exponents'],region['residues']
        assert len(es)==len(rs)==len(core.ps)
        assert all(type(e) is int and 0<=e<=E for e,E in zip(es,core.ds))
        assert all(type(r) is int and 0<=r<labels(p,E,e) for p,E,e,r in zip(core.ps,core.ds,es,rs))
        unions=region.get('unions',{})
        for k,vs in unions.items():
            i=int(k)
            assert str(i)==k and 0<=i<len(core.ps) and vs and len(vs)==len(set(vs))
            assert all(type(v) is int and 0<=v<labels(core.ps[i],core.ds[i],es[i]) for v in vs)
        return core.profile(es,rs,residual if conditioning=='current' else None,epoch,unions)

    pending=list(result['derived_regions'])

    def apply_derived(item):
        nonlocal epoch
        pid=check_region(item['region'])
        s=by_prime[item['source_prime']]
        assert all(contains(s,i,pid) for i in sources[s]['root_inputs'])
        region=item['region']
        residual[core.region_index(region)]=False
        epoch+=1

    for source in result['sources']:
        while pending and pending[0].get('after_sources',len(result['sources']))==len(sources):
            apply_derived(pending.pop(0))
        q=source['root_prime']
        assert prime(q) and max(core.ps)<q<=result['largest_prime'] and q not in by_prime
        assert not sources or q>sources[-1]['root_prime']
        if 'partial_joint_certificate' in source:
            from joint_root_repair import validate_packing
            packed=source['partial_joint_certificate']
            assert packed['core']==certificate
            assert partial_joint_source(packed,source['assigned_digits'])==source
            validate_packing(packed)
            by_prime[q]=len(sources)
            sources.append(source)
            continue
        pid=check_region(source['region'])
        if 'partial_borrowing_certificate' in source:
            from partial_borrow_transport import validate as validate_transport,source as transport_source
            record=source['partial_borrowing_certificate']
            assert record['core']==certificate and transport_source(record)==source
            prior=sources[by_prime[11]]
            assert prior['root_inputs']==list(range(10))
            assert all(n['kind']=='row_union' and n['levels']==record['levels'] for n in prior['nodes'])
            assert [n['rows'] for n in prior['nodes']]==record['old_prime11_inputs']
            validate_transport(record)
            by_prime[q]=len(sources)
            sources.append(source)
            assert all(contains(len(sources)-1,i,pid) for i in source['root_inputs'])
            residual[core.region_index(source['region'])]=False
            epoch+=1
            continue
        if 'joint_certificate' in source:
            from joint_root_repair import validate as validate_joint
            joint=source['joint_certificate']
            assert joint['core']==certificate and joint_source(joint)==source
            validate_joint(joint)
            by_prime[q]=len(sources)
            sources.append(source)
            assert all(contains(len(sources)-1,i,pid) for i in source['root_inputs'])
            residual[core.region_index(source['region'])]=False
            epoch+=1
            continue
        configs={int(p):(v['fixed_exponent'],v['prefix'],tuple(v['digits']))
                 for p,v in source['configs'].items()}
        assert set(configs)==set(core.ps)
        # Each omitted core branch must have no base-uncovered point in
        # the entire target region, not merely no sampled point.
        for p,E,values in zip(core.ps,core.ds,core.profiles[pid]):
            f,a,ds=configs[p]
            assert type(f) is int and 0<=f<E and 0<=a<p**f
            assert all((x-a)%p**f==0 for x in values)
            actual=tuple(d for d in range(1,p) if any(
                (x-a-d*p**(k-1))%p**k==0
                for x in values for k in range(f+1,E+1)))
            assert actual==ds and len(ds)>0
        borrowed={int(p):tuple(ds) for p,ds in source['borrowed'].items()}
        for p,ds in borrowed.items():
            assert p in by_prime and p not in configs and len(ds)==len(set(ds))
            old=by_prime[p]
            for d in ds:
                assert 1<=d<p and contains(old,sources[old]['root_inputs'][d-1],pid)
                dependencies+=1
        order=source['order']
        assert len(order)==len(set(order)) and all(prime(p) and p<q for p in order)
        assert set(core.ps)<=set(order) and set(borrowed)<=set(order)
        nodes=source['nodes']
        assert nodes[0]==dict(kind='unit',grade=1)
        active=[0]
        cursor=1
        for p in order:
            if result.get('forest_supply',False) and p in core.ps:
                i=core.ps.index(p)
                E=core.ds[i]
                values=core.profiles[pid][i]
                following=active[:]
                used=set()
                while cursor<len(nodes) and nodes[cursor].get('prime')==p:
                    node=nodes[cursor]
                    assert node['kind']=='forest' and node['branches']
                    covered=set()
                    grades=[]
                    for b in node['branches']:
                        e,a,child=b['exponent'],b['prefix'],b['child']
                        assert type(e) is int and 1<=e<=E and child in active
                        assert (e,child) not in used
                        used.add((e,child))
                        if b['kind']=='fixed':
                            assert e<E and 0<=a<p**e and b['digit']==0
                            covered.update(x for x in values if (x-a)%p**e==0)
                        else:
                            assert b['kind']=='tail' and e==E and 0<=a<p**(E-1)
                            assert 1<=b['digit']<p
                            covered.update(x for x in values if (x-a-b['digit']*p**(E-1))%p**E==0)
                        grades.append(p**e*nodes[child]['grade'])
                    assert all(x in covered for x in values)
                    assert node['grade']==min(grades)
                    following.append(cursor)
                    cursor+=1
                active=following
                continue
            f,a,required=configs.get(p,(0,0,tuple(range(1,p))))
            free=borrowed.get(p,())
            ranked=sorted(active,key=lambda i:(-nodes[i]['grade'],i))
            following=active[:]
            for e in range(1,f+1):
                for child in ranked:
                    expected=dict(kind='fixed',prime=p,exponent=e,prefix=a%p**e,
                                  children=[child],grade=nodes[child]['grade']*p**e)
                    assert nodes[cursor]==expected
                    following.append(cursor)
                    cursor+=1
            digits=[d for d in required if d not in free]
            k=len(digits)
            assert k>0
            for j in range(len(ranked)//k):
                children=ranked[j*k:(j+1)*k]
                expected=dict(kind='tail',prime=p,exponent=f+1,prefix=a,digits=digits,
                              borrowed=list(free),children=children,
                              grade=min(nodes[i]['grade'] for i in children)*p**(f+1))
                assert nodes[cursor]==expected
                following.append(cursor)
                cursor+=1
            active=following
        assert cursor==len(nodes)
        selected=source['root_inputs']
        assert len(selected)==len(set(selected))==q-1 and set(selected)<=set(active)
        assert all(q*nodes[i]['grade']>=result['minimum'] for i in selected)
        by_prime[q]=len(sources)
        sources.append(source)
        region=source['region']
        residual[core.region_index(region)]=False
        epoch+=1
    for item in pending:
        assert item.get('after_sources',len(sources))==len(sources)
        apply_derived(item)
    gap=int(residual.sum())
    assert gap==result['guaranteed_uncovered_states']
    contains.cache_clear()
    return dict(status='VERIFIED_COMPLETE_REGION_SCHEMA' if gap==0 else 'VERIFIED_PARTIAL_REGION_SCHEMA',
                kernel_states=core.total,guaranteed_uncovered_states=gap,
                root_stages=len(sources),checked_borrowed_inputs=dependencies,
                checked_core_branch_profiles=sum(not any(k in s for k in ('joint_certificate','partial_joint_certificate','partial_borrowing_certificate')) for s in sources)*len(core.ps),
                checked_joint_roots=sum('joint_certificate' in s for s in sources),
                checked_partial_joint_roots=sum('partial_joint_certificate' in s for s in sources),
                checked_partial_borrowing_stages=sum('partial_borrowing_certificate' in s for s in sources),
                scope='Zero gap certifies the finite completed distinct cover. Nonzero gap concerns only guaranteed coverage.')


def assignments(core,marginal,es,per_shape,union_targets):
    """Single cylinders, plus unions of a p-class and a different p^2-class."""
    result=[]
    flat=marginal.reshape(-1)
    for index in np.argsort(-flat,kind='stable')[:per_shape]:
        rs=tuple(map(int,np.unravel_index(int(index),marginal.shape)))
        result.append((int(flat[index]),rs,{}))
    if union_targets and marginal.size<=252:
        for i,(p,E,e) in enumerate(zip(core.ps,core.ds,es)):
            if p not in (2,3) or e!=2 or e>=E:
                continue
            for a in range(p):
                for b in range(p*p):
                    if b%p==a:
                        continue
                    vs=sorted(list(range(a,p*p,p))+[b])
                    merged=np.take(marginal,vs,axis=i).sum(axis=i,keepdims=True)
                    vector=merged.reshape(-1)
                    for index in np.argsort(-vector,kind='stable')[:per_shape]:
                        rs=tuple(map(int,np.unravel_index(int(index),merged.shape)))
                        result.append((int(vector[index]),rs,{str(i):vs}))
    return sorted(result,key=lambda x:-x[0])


def search(certificate,largest=113,max_volume=2520,per_shape=4,seconds=180,conditioning='current',branch_sensitive=False,forest_supply=False,union_targets=False,initial_joint=None,initial_transport=None,resume_result=None,initial_packing=None,packing_digits=None):
    start=time.monotonic()
    core=Core(certificate)
    residual=core.residual.copy()
    shapes=[]
    for es in product(*(range(E+1) for E in core.ds)):
        volume=prod(labels(p,E,e) for p,E,e in zip(core.ps,core.ds,es))
        if volume<=max_volume:
            shapes.append((volume,es))
    shapes.sort()
    catalog=Catalog(core,branch_sensitive)
    derived=[]
    trace=[]
    stopped='ROOTS_EXHAUSTED'
    epoch=0
    resumed_after=None
    if resume_result is not None:
        assert initial_joint is None and initial_transport is None and initial_packing is None
        assert resume_result['core']==certificate
        assert resume_result.get('forest_supply',False)==forest_supply
        assert resume_result.get('branch_sensitive',False)==branch_sensitive
        assert resume_result.get('profile_conditioning','core')==conditioning
        validate(resume_result)
        for source in resume_result['sources']:
            catalog.append(source)
            if 'region' in source:
                residual[core.region_index(source['region'])]=False
        derived=list(resume_result['derived_regions'])
        for item in derived:
            residual[core.region_index(item['region'])]=False
        assert int(residual.sum())==resume_result['guaranteed_uncovered_states']
        epoch=len(catalog.sources)+len(derived)
        resumed_after=max(catalog.by_prime,default=max(core.ps))
    if initial_joint is not None:
        from joint_root_repair import validate as validate_joint
        assert initial_joint['core']==certificate
        validate_joint(initial_joint)
        q=initial_joint['root_prime']
        assert q==next(p for p in primes_through(largest) if p>max(core.ps))
        catalog.append(joint_source(initial_joint))
        residual[core.region_index(initial_joint['region'])]=False
        epoch+=1
    if initial_packing is not None:
        from joint_root_repair import validate_packing
        assert initial_joint is None and initial_packing['core']==certificate
        assert initial_packing['root_prime']==next(p for p in primes_through(largest) if p>max(core.ps))
        validate_packing(initial_packing)
        catalog.append(partial_joint_source(initial_packing,packing_digits))
        epoch+=1
    if initial_transport is not None:
        from partial_borrow_transport import validate as validate_transport,source as transport_source
        assert 11 in catalog.by_prime
        assert initial_transport['core']==certificate
        assert initial_transport['old_prime11_inputs']==[n['rows'] for n in catalog.sources[catalog.by_prime[11]]['nodes']]
        validate_transport(initial_transport)
        catalog.append(transport_source(initial_transport))
        residual[core.region_index(initial_transport['region'])]=False
        epoch+=1
    for q in primes_through(largest):
        if q<=max(core.ps) or q in catalog.by_prime or resumed_after is not None and q<=resumed_after:
            continue
        if time.monotonic()-start>seconds:
            stopped='TIME_LIMIT'
            break
        while True:
            best=None
            for volume,es in shapes:
                bound=core.total//volume*(4 if union_targets and volume<=252 else 1)
                if best is not None and bound<=best['gain']:
                    continue
                marginal=core.marginal(residual,es)
                for gain,rs,unions in assignments(core,marginal,es,per_shape,union_targets):
                    if gain==0 or best is not None and gain<=best['gain']:
                        continue
                    pid=core.profile(es,rs,residual if conditioning=='current' else None,epoch,unions)
                    configs=core.configs(pid)
                    forests=core.forest_states(pid) if forest_supply else None
                    borrowed=catalog.borrowed(pid)
                    already=next((p for p,ds in borrowed.items() if len(ds)==p-1),None)
                    good,order=best_order(q,configs,borrowed,certificate['minimum'],forests) if already is None else (None,[])
                    if already is None and good<q-1:
                        continue
                    best=dict(gain=gain,region=dict(exponents=list(es),residues=list(rs)),
                              configs=configs,borrowed=borrowed,order=order,
                              available=good,already=already,forests=forests)
                    if unions:
                        best['region']['unions']=unions
            if best is None:
                break
            region=best['region']
            if best['already'] is not None:
                derived.append(dict(source_prime=best['already'],region=region,after_sources=len(catalog.sources)))
                residual[core.region_index(region)]=False
                epoch+=1
                if not residual.any():
                    break
                continue
            values,nodes=supply(best['order'],best['configs'],best['borrowed'],forest_states=best['forests'])
            good=[(g,i) for g,i in values if q*g>=certificate['minimum']]
            assert len(good)==best['available'] and len(good)>=q-1
            source=dict(root_prime=q,region=region,order=best['order'],
                        configs={str(p):dict(fixed_exponent=f,prefix=a,digits=list(ds))
                                 for p,(f,a,ds) in best['configs'].items()},
                        borrowed={str(p):list(ds) for p,ds in best['borrowed'].items()},
                        nodes=nodes,root_inputs=[i for _,i in good[:q-1]],available_inputs=len(good))
            catalog.append(source)
            before=int(residual.sum())
            residual[core.region_index(region)]=False
            after=int(residual.sum())
            epoch+=1
            assert before-after==best['gain']
            info=dict(root_prime=q,region=region,gain=best['gain'],available_inputs=len(good),
                      guaranteed_uncovered_states=after,uncovered_fraction=after/core.total,
                      borrowed_inputs=sum(len(ds) for ds in best['borrowed'].values()),
                      elapsed_seconds=round(time.monotonic()-start,3))
            trace.append(info)
            print(json.dumps(info),flush=True)
            break
        if not residual.any():
            stopped='CANDIDATE_COMPLETE'
            break
    result=dict(minimum=certificate['minimum'],core=certificate,sources=catalog.sources,
                derived_regions=derived,guaranteed_uncovered_states=int(residual.sum()),
                kernel_states=core.total,method='exact_core_region_relative_arrows',
                largest_prime=largest,maximum_region_volume=max_volume,targets_per_shape=per_shape,
                profile_conditioning=conditioning,
                branch_sensitive=branch_sensitive,
                forest_supply=forest_supply,
                union_targets=union_targets,
                resumed_after_prime=resumed_after,
                search_outcome=stopped,trace=trace,elapsed_seconds=time.monotonic()-start,
                recipe='Core rows and each source DAG specify the own modulus families. A source root uses its q-1 inputs in digit order at q-tail exponent 1. Internal fixed/tail nodes use their explicit prefixes and digits. Missing core digits are justified by the exhaustive core residual; borrowed noncore digits reference earlier sources. Close all own tail occurrences with distinct auxiliary primes above largest_prime in core-row, then source/preorder order. These finite completions preserve the checked resource distinctions and minimum.')
    catalog.contains.cache_clear()
    result['verification']=validate(result)
    result['status']=result['verification']['status']
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('core')
    parser.add_argument('--largest-prime',type=int,default=113)
    parser.add_argument('--max-volume',type=int,default=2520)
    parser.add_argument('--per-shape',type=int,default=4)
    parser.add_argument('--seconds',type=float,default=180)
    parser.add_argument('--conditioning',choices=['current','core'],default='current')
    parser.add_argument('--branch-sensitive',action='store_true')
    parser.add_argument('--forest-supply',action='store_true')
    parser.add_argument('--union-targets',action='store_true')
    parser.add_argument('--initial-joint')
    parser.add_argument('--initial-transport')
    parser.add_argument('--resume')
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    initial=json.loads(Path(args.initial_joint).read_text()) if args.initial_joint else None
    transport=json.loads(Path(args.initial_transport).read_text()) if args.initial_transport else None
    resumed=json.loads(Path(args.resume).read_text()) if args.resume else None
    result=search(json.loads(Path(args.core).read_text()),args.largest_prime,args.max_volume,args.per_shape,args.seconds,args.conditioning,args.branch_sensitive,args.forest_supply,args.union_targets,initial,transport,resumed)
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['verification'],indent=2))
