"""Independent exact checks for the reachable-product support refinement."""
from fractions import Fraction as F
from bisect import bisect_right
from moment_product_bound import MomentProductBound
from verify_minimum_modulus import S,ceildiv
from verify_stoploss import L
from test_stoploss import exact_product,exact_excess


def check():
    comparisons=0;paths=0
    for R in (2,4,16):
        b=MomentProductBound(100,R,'reachable')
        states={(1,1)}
        for _ in range(12):
            new=set()
            for rounded,true in states:
                for h in range(1,min(108,100//rounded)+1):
                    j=bisect_right(b.bins,rounded*h)-1
                    assert b.bins[j]<=true*h<=b.upper[j]
                    new.add((b.bins[j],true*h));paths+=1
            states=new
        ps=[2,3,5,7];ds=[0,L//4,L//5,L//3]
        dist,mean=exact_product(ps,ds,100)
        for q,D in zip(ps,ds):b.update(q,D)
        pref=b.prefixes();mean_up=ceildiv(mean.numerator*S,mean.denominator)
        for j in range(301):
            t=F(j,3);got=F(b.excess_upper(t.numerator*S//t.denominator,mean_up,pref),S)
            assert got>=exact_excess(t,dist,mean)
            comparisons+=1
    print('PASS:',comparisons,'exact excess bounds;',paths,'retained-path support checks.')


if __name__=='__main__':check()
