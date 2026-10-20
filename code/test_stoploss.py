"""Independent rational tests of the new product and conditional bounds."""
from fractions import Fraction as F
from itertools import product as cartesian
from random import Random
from verify_stoploss import ProductBound,SmallBound,L
from verify_minimum_modulus import S,ceildiv


def exact_product(primes,deltas,cap):
    dist={1:F(1)};mean=F(1)
    for q,D in zip(primes,deltas):
        nu=F(L,L-D);new={}
        for n,u in dist.items():
            for h in range(1,cap//n+1):
                prob=1-nu/q if h==1 else nu*F(q-1,q**h)
                new[n*h]=new.get(n*h,F(0))+u*prob
        dist=new;mean*=1+nu/(q-1)
    return dist,mean


def exact_excess(t,dist,mean):
    return mean-t+sum((t-n)*u for n,u in dist.items() if n<t)


def check():
    cases=0
    ps=[5,7,11];ds=[0,L//4,L//5]
    dist,mean=exact_product(ps,ds,40)
    for resolution in (2,8,100):
        b=ProductBound(40,resolution)
        for q,D in zip(ps,ds):b.update(q,D)
        pref=b.prefixes();mu=ceildiv(mean.numerator*S,mean.denominator)
        for t in [F(n,3) for n in range(1,120)]:
            tlo=t.numerator*S//t.denominator
            upper=b.excess_upper(tlo,mu,pref)
            assert F(upper,S)>=exact_excess(t,dist,mean)
            cases+=1
    # Compare the entire finite-state conditional bound with an exact sum.
    dist,mean=exact_product([5,7],ds[:2],40)
    b=ProductBound(40,100)
    for q,D in zip([5,7],ds[:2]):b.update(q,D)
    for K,p,D in ((80,13,L//4),(100,17,L//3),(20,19,L//2),(80,13,L//10)):
        sm=SmallBound(3,2);d=F(D,L);total=F(0);keptmean=F(0)
        for a,bb in cartesian(range(4),range(3)):
            tau=(a+1)*(bb+1);prob=F(2,2**(a+1)*3**(bb+1));keptmean+=tau*prob
            excluded=F(0);power=p
            while power<K:
                count=sum(2**aa*3**bbb*power<K for aa in range(a+1) for bbb in range(bb+1))
                excluded+=F((p-1)*count,power);power*=p
            threshold=((p-1)*d+excluded)/tau
            total+=prob*tau*exact_excess(threshold,dist,mean)
        total+=(3-keptmean)*mean
        exact=total/((p-1)*(1-d))
        got=sm.loss_upper(K,p,D,ceildiv(mean.numerator*S,mean.denominator),b)
        assert F(got,S)>=exact
        assert F(got,S)-exact<F(1,10**25)
        cases+=1
    # Dependent coordinates, conflicting residues, and repeated cylinders.
    rng=Random(40000)
    comparison_cases=0
    for _ in range(128):
        terms=[]
        for j in range(8):
            a=rng.randrange(3);bb=rng.randrange(3)
            terms.append((a,bb,rng.randrange(2**a),rng.randrange(3**bb),F(rng.randrange(1,5),rng.randrange(1,6))))
        for threshold in (F(1,2),F(1),F(2)):
            actual=F(0)
            for x in range(3):
                for offset in range(6):
                    y=(offset+2*x)%9
                    z=sum(c for a,bb,r,s,c in terms if x%(2**a)==r and y%(3**bb)==s)
                    actual+=max(F(0),z-threshold)/18
            dominant=F(0)
            for a in range(3):
                for bb,prob in enumerate((F(1,2),F(1,3),F(1,6))):
                    z=sum(c for u,v,r,s,c in terms if a>=u and bb>=v)
                    dominant+=prob*max(F(0),z-threshold)/3
            assert actual<=dominant
            comparison_cases+=1
    print(f'PASS: {cases} independent exact-rational stop-loss checks; {comparison_cases} dependent-coordinate convex comparisons.')


if __name__=='__main__':check()
