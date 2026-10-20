"""Exact retained path checks with deliberately enlarged exponent remainders."""
from fractions import Fraction as F
import bisect
from moment_product_bound import MomentProductBound
from verify_minimum_modulus import S
from verify_stoploss import L


if __name__=='__main__':
    count=0
    for cutoff in (20,2000,2*10**20):
        b=MomentProductBound(100,8,'reachable');b.tail_cutoff=cutoff
        mass=[F(0)]*len(b.bins);mass[1]=F(1);moment=mass.copy();omitted=F(0)
        for q,D in ((2,L//100),(3,L//5),(5,L//4),(11,L//3)):
            v=F(L,L-D);H=1
            while q**H<cutoff:H+=1
            omitted+=v/F(q**H)
            new=[x*(1-v/q) for x in mass];mom=[x*(1-v/q) for x in moment]
            for h in range(2,H+1):
                prob=v*F(q-1,q**h)
                for j,low in enumerate(b.bins):
                    if j and low*h<=b.bins[-1]:
                        k=bisect.bisect_right(b.bins,low*h)-1
                        new[k]+=mass[j]*prob;mom[k]+=moment[j]*h*prob
            mass,moment=new,mom;b.update(q,D)
            assert all(F(x,S)>=y for x,y in zip(b.mass,mass))
            assert all(F(x,S)<=y for x,y in zip(b.moment,moment))
            assert F(b.error,S)>=omitted
            count+=1
    print('Product truncation:',count,'exact retained mass/moment recursions and omitted-exponent bounds')
