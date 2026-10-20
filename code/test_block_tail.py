"""Independent rational checks of the accelerated future-product envelope."""
from fractions import Fraction as Q
from itertools import product
from block_tail_bound import GeometricEnvelope,block_values,chord_sum,exp_interval
from verify_minimum_modulus import I,S
from verify_stoploss import L


def enclosed(interval,value):
    assert Q(interval.lo,S)<=value<=Q(interval.hi,S),(interval.lo,interval.hi,value)


def run():
    probabilities=losses=chords=exponentials=0
    for D in (0,L//5,L//3,L//2):
        qs=[5,7,11];nu=Q(L,L-D);H=8;env=GeometricEnvelope(H)
        dist=[Q(1)]+[Q(0)]*H;mean=Q(1)
        for q in qs:
            rho=nu/q;env.update(q,D)
            for a in range(1,12):assert nu/q**a<=rho**a
            dist=[sum(dist[j]*(1-rho)*rho**(k-j) for j in range(k+1)) for k in range(H+1)]
            mean*=(1-rho)/(1-2*rho)
            got,remainder=env.probabilities()
            for g,v in zip(got,dist):enclosed(g,v);probabilities+=1
            enclosed(env.mean,mean)
            assert Q(remainder,S)>=mean-sum((1<<j)*v for j,v in enumerate(dist))
        points=[1,2,5,11]
        bounds=block_values(env,points,qs[0],len(qs),D)
        delta=Q(D,L)
        # Sum the actual per-prime losses over a finite part of the ORIGINAL
        # exponent law; all discarded terms in this comparison are positive.
        for y,bound in zip(points,bounds):
            actual=Q(0)
            for index,q in enumerate(qs):
                for exponents in product(range(6),repeat=index):
                    z=y;prob=Q(1)
                    for old,e in zip(qs[:index],exponents):
                        prob*=1-nu/old if e==0 else nu*(old-1)/old**(e+1)
                        z*=e+1
                    actual+=prob*max(0,Q(z,q-1)-delta)/(1-delta)
            assert Q(bound,S)>=actual;losses+=1
    for b,u in ((1,3),(3,7),(10,19)):
        for d in (Q(1,3),Q(2),Q(5)):
            masses=[Q(1,7),Q(2,7),Q(1,7)]
            xs=[Q(b),Q(b+u,2),Q(u)];p=sum(masses);m=sum(x*v for x,v in zip(xs,masses))
            p_hi=(p.numerator*S+p.denominator-1)//p.denominator
            m_lo=m.numerator*S//m.denominator
            state=dict(bins=[b],upper=[u],mass=[str(p_hi)],moment=[str(m_lo)],mass_error='1',moment_error='1')
            def f(x):return max(0,x-d)+max(0,2*x-3*d)
            lo=I.rat(f(Q(b)).numerator,f(Q(b)).denominator).hi
            hi=I.rat(f(Q(u)).numerator,f(Q(u)).denominator).hi
            got=chord_sum(state,[lo],[hi]);actual=sum(v*f(x) for x,v in zip(xs,masses))
            assert Q(got,S)>=actual;chords+=1
    for x in (Q(0),Q(1,10),Q(1),Q(5)):
        got=exp_interval(I.rat(x.numerator,x.denominator));term=total=Q(1)
        for j in range(1,201):term*=x/j;total+=term
        following=term*x/201;tail=following/(1-x/202)
        assert Q(got.lo,S)<=total and Q(got.hi,S)>=total+tail;exponentials+=1
    # A pure second-power moment equals the exact geometric generating
    # function and is dominated by the exponential growth bound.
    env=GeometricEnvelope(8);exact=Q(1)
    for q in (11,13,17):
        env.update(q,L//3);rho=Q(3,2*q);exact*=(1-rho)/(1-4*rho)
    assert Q(env.growth(8,4).hi,S)>=exact
    result=dict(status='PASSED',rational_probability_enclosures=probabilities,
                actual_loss_comparisons=losses,chord_comparisons=chords,
                exponential_enclosures=exponentials,growth_comparisons=1)
    return result


if __name__=='__main__':
    import argparse,hashlib,json
    from pathlib import Path
    ap=argparse.ArgumentParser();ap.add_argument('--output');a=ap.parse_args();result=run()
    result['source_sha256']={str(p).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in (Path('code/block_tail_bound.py'),Path(__file__))}
    if a.output:Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
