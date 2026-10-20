"""Independent rational recursion and actual finite distorted-measure checks."""
import math
import random
from fractions import Fraction as Q
from recursive_atom_bound import RecursiveAtoms
from verify_minimum_modulus import S, ONE, I, moments, factorizations
from verify_stoploss import L


def run():
    count=0; comparisons=0
    for seed in range(8):
        rng=random.Random(seed);K=12+(seed%4)*3;cap=160;degree=1+seed//4
        atoms=RecursiveAtoms(K,cap,degree);measure=[Q(1)];mod=1;nus={};C=T=ONE
        B={1:Q(1)};F={1:Q(1)}
        for p,e,D in [(2,2,L//8 if seed%2 else 0),(3,1,L//5),(5,1,L//4),(7,1,L//3)]:
            v=Q(L,L-D);delta=Q(D,L)
            # Independently form all caps simultaneously from the old map.
            new={};oldF=dict(F)
            def h(m):
                power=p
                while m*power<K:power*=p
                return Q(p,power*(p-1))
            for d in oldF:
                if d*p>cap:continue
                product=math.prod(1+u/Q(q-1) for q,u in nus.items())*B[d]
                for q,u in nus.items():
                    a=0;t=d
                    while t%q==0:a+=1;t//=q
                    if a:product*=(Q(a+1)+Q(1,q-1))/(1+u/Q(q-1))
                H=product/Q(p-1)
                for m,b in B.items():
                    if m*p<K:H-=(Q(1,p-1)-h(m))*B[d]*b/B[math.gcd(d,m)]
                for m in oldF:
                    n=math.lcm(d,m)
                    if n<=cap:H-=h(m)*(B[n]-oldF[n])
                a=min(v*oldF[d],oldF[d]+v*H)
                if degree==2:
                    H2=math.prod(1+u*Q(3*q-1,(q-1)**2) for q,u in nus.items())*B[d]/(p-1)**2
                    for q,u in nus.items():
                        exponent=0;t=d
                        while t%q==0:exponent+=1;t//=q
                        if exponent:
                            H2*=(Q((exponent+1)**2)+Q(2*exponent+3,q-1)+Q(2,(q-1)**2))/(1+u*Q(3*q-1,(q-1)**2))
                    a=min(a,oldF[d]+H+v*H2)
                power=p
                while d*power<=cap:new[d*power]=a/power;power*=p
            nus[p]=v
            for n in range(1,atoms.limit+1):
                t=n;value=Q(1)
                for q,u in nus.items():
                    if t%q==0:
                        value*=u
                        while t%q==0:t//=q;value/=q
                if t==1:B[n]=value
            F.update(new)
            m1,m2=moments(K,p,C,T,{q:I.rat(u.numerator,u.denominator) for q,u in nus.items() if q<p},factorizations(K))
            gains=atoms.moment_gains(p)
            newmod=mod*p**e
            allowed=[d for d in range(K,newmod+1) if newmod%d==0 and d%p==0]
            residues={d:rng.randrange(d) for d in allowed}
            bad=[any(x%d==r for d,r in residues.items()) for x in range(newmod)]
            nextmeasure=[Q(0)]*newmod;actual=[Q(0),Q(0)]
            for x in range(mod):
                fibre=list(range(x,newmod,mod));alpha=Q(sum(bad[y] for y in fibre),len(fibre))
                actual[0]+=measure[x]*alpha;actual[1]+=measure[x]*alpha**2
                for y in fibre:
                    if bad[y]:density=(alpha-delta)/(alpha*(1-delta)) if alpha>delta else Q(0)
                    else:density=1/(1-min(alpha,delta))
                    nextmeasure[y]=measure[x]*density/len(fibre)
            assert sum(nextmeasure)==1
            assert actual[0]<=Q(m1.hi-gains[1],S)
            assert actual[1]<=Q(m2.hi-gains[2],S)
            atoms.update(p,D)
            for n,f in F.items():
                assert Q(atoms.F[n],S)>=f
                assert Q(atoms.F[n],S)-f<Q(1,10**32)
                count+=1
                if newmod%n==0:
                    for a in range(n):
                        assert sum(nextmeasure[a::n])<=Q(atoms.F[n],S)
                        comparisons+=1
            measure=nextmeasure;mod=newmod
            vi=I.rat(v.numerator,v.denominator)
            C=C*(ONE+vi*I.rat(1,p-1));T=T*(ONE+vi*I.rat(3*p-1,(p-1)**2))
    print('Recursive atoms:',count,'rational caps,',comparisons,'actual residue comparisons, 64 moment checks')


if __name__=='__main__':run()
