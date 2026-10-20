"""Exact partition checks of the stratified moment/overflow assembly."""
from fractions import Fraction as F
from itertools import product
from stratified_moment_tail import assemble
from refine_tail_parameter import prepare
from verify_prefix_refinement import scaled
from verify_minimum_modulus import I,S,ceildiv
from moment_product_bound import MomentProductBound
from verify_stoploss import L
import bisect


if __name__=='__main__':
    checks=0
    for cutoff,rounding in product((5,10,30),(0,1,37)):
        atoms=[(1,F(1,5)),(2,F(1,4)),(5,F(1,5)),(9,F(1,5)),(25,F(1,10)),(100,F(1,20))]
        class Bins:pass
        b=Bins();b.bins=[1,4,8,16];b.upper=[3,7,15,31]
        masses=[F(0)]*4;moments=[F(0)]*4;retained=[[] for _ in range(4)];outside=[]
        for y,p in atoms:
            if y>=cutoff:outside.append((y,p));continue
            j=max(j for j,lo in enumerate(b.bins) if lo<=y)
            masses[j]+=p;moments[j]+=y*p;retained[j].append((y,p))
        b.mass=[ceildiv(p.numerator*S,p.denominator)+rounding for p in masses]
        b.moment=[max(0,p.numerator*S//p.denominator-rounding) for p in moments]
        orders=[(12,4),(16,4)];coef=[F(1,1000),F(1,10000)]
        coeff=[ceildiv(c.numerator*S,c.denominator) for c in coef]
        U=[sum(p*y**s for y,p in atoms) for s in (3,4)]
        got,details=assemble(b,orders,coeff,[I.rat(u.numerator,u.denominator) for u in U],rounding+1,rounding+1)
        table=prepare(details)
        replay=sum(min(ceildiv(c*u,S) for c,u in zip(coeff,cell)) for cell in zip(*table))
        assert replay+scaled(details['stratified_overflow_tail_upper'])==got
        exact=sum(min(c*sum(p*y**s for y,p in cell) for c,s in zip(coef,(3,4))) for cell in retained+[outside])
        assert F(got,S)>=exact,(cutoff,rounding,got,exact)
        for j,group in enumerate(details['stratified_selected_groups']):
            true=sum(p*y**(orders[j][0]//orders[j][1])
                     for i,cell in enumerate(retained)
                     if details['stratified_bin_order_indices'][i]==j for y,p in cell)
            assert F(group['moment_upper'])>=true
        for new in ((F(1,2000),F(1,30000)),(F(1,100),F(1,20000))):
            bound=sum(c*F(g['moment_upper']) for c,g in zip(new,details['stratified_selected_groups']))
            true=sum(new[j]*sum(p*y**(orders[j][0]//orders[j][1]) for y,p in cell)
                     for cell,j in zip(retained,details['stratified_bin_order_indices']))
            assert bound>=true
            cc=[ceildiv(c.numerator*S,c.denominator) for c in new]
            switched=sum(min(ceildiv(c*u,S) for c,u in zip(cc,cell)) for cell in zip(*table))
            exact_min=sum(min(c*sum(p*y**s for y,p in cell) for c,s in zip(new,(3,4))) for cell in retained)
            assert F(switched,S)>=exact_min
        partial,_=assemble(b,orders,coeff,[I.rat(U[0].numerator,U[0].denominator),None],rounding+1,rounding+1)
        assert partial>=got
        checks+=1
    print('PASS:',checks,'exact stratified moment partitions, selected-group enclosures, 18 regrouped-coefficient and 18 reselected-order comparisons; exact state replay.')
    checks=0
    for tilt in (0,L//3):
        full=MomentProductBound(40,8,'reachable');N=len(full.bins)
        mass=[F(0)]*N;mass[1]=F(1);moment=mass.copy();C=F(1)
        for n,q in enumerate((2,3,5,7,11,13),1):
            nu=F(L,L-tilt);zero=1-nu/q
            mm=[x*zero for x in mass];ww=[x*zero for x in moment]
            H=1
            while q**H<2*10**32:H+=1
            for h in range(2,H+1):
                prob=nu*F(q-1,q**h)
                for j,b in enumerate(full.bins):
                    if j and b*h<=full.bins[-1]:
                        k=bisect.bisect_right(full.bins,b*h)-1
                        mm[k]+=mass[j]*prob;ww[k]+=moment[j]*h*prob
            mass,moment=mm,ww;C*=1+nu/(q-1);full.update(q,tilt)
            ep=n*N*109;em=ceildiv(ep*C.numerator,C.denominator)
            for a,b,c,d in zip(full.mass,full.moment,mass,moment):
                assert 0<=F(a,S)-c<=F(ep,S)
                assert 0<=d-F(b,S)<=F(em,S)
            checks+=1
    print('PASS:',checks,'exact retained-bin mass and moment rounding-error enclosures.')
