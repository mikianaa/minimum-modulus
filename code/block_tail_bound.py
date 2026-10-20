"""Integer block tail bounds conditional on a completed terminal report.

For future primes q, E_q is dominated by a geometric variable of parameter
rho_q=nu/q; E_q+1 <= 2**E_q. Positive power-sum recurrences enclose the
product envelope, without propagating every bin at every prime.
"""
import argparse,hashlib,json
from pathlib import Path
from verify_minimum_modulus import I,ONE,S,ceildiv,primes_to,fixed_text
from verify_stoploss import L
from verify_prefix_refinement import scaled
from fractional_moment_tail import fractional_tail_bound,power
from refine_tail_parameter import prepare


def exp_interval(x):
    assert x.lo>=0
    term=total=ONE
    for j in range(1,10000):
        term=term*x/I.rat(j);total=total+term
        if (j+2)*S>2*x.hi:
            following=term*x/I.rat(j+1)
            tail=following/(ONE-x/I.rat(j+2))
            if tail.hi<=10:return I(total.lo,total.hi+tail.hi)
    raise AssertionError('exponential series did not converge')


class GeometricEnvelope:
    def __init__(self,H=16):
        assert H>=1
        self.H=H;self.powers=[I(0) for _ in range(H+1)]
        self.zero=ONE;self.mean=ONE;self.rhomax=I(0);self.count=0

    def update(self,q,D):
        assert 0<=D<=L//2
        rho=I.rat(L,q*(L-D));assert 2*rho.hi<S
        self.rhomax=I(max(self.rhomax.lo,rho.lo),max(self.rhomax.hi,rho.hi))
        term=rho
        for k in range(1,self.H+1):
            self.powers[k]=self.powers[k]+term;term=term*rho
        self.zero=self.zero*(ONE-rho)
        self.mean=self.mean*(ONE-rho)/(ONE-I.rat(2)*rho);self.count+=1

    def probabilities(self):
        h=[ONE]
        for j in range(1,self.H+1):
            value=I(0)
            for k in range(1,j+1):value=value+self.powers[k]*h[j-k]
            h.append(value/I.rat(j))
        probabilities=[self.zero*v for v in h]
        remainder=max(0,self.mean.hi-sum((1<<j)*v.lo for j,v in enumerate(probabilities)))
        return probabilities,remainder

    def growth(self,num,den):
        a=power(2,num,den)
        assert (a*self.rhomax).hi<S
        return exp_interval((a-ONE)*self.powers[1]/(ONE-a*self.rhomax))


def block_values(envelope,points,first_prime,count,D):
    probabilities,remainder=envelope.probabilities()
    threshold=D*(first_prime-1)*S//L;values=[]
    for y in points:
        numerator=y*remainder
        for j,prob in enumerate(probabilities):
            excess=max(0,y*(1<<j)*S-threshold)
            numerator+=ceildiv(prob.hi*excess,S)
        values.append(ceildiv(count*numerator*L,(first_prime-1)*(L-D)))
    return values


def chord_sum(state,lower_values,upper_values):
    total=0;ep=int(state['mass_error']);em=int(state['moment_error'])
    for b,u,p,m,gb,gu in zip(state['bins'],state['upper'],state['mass'],state['moment'],lower_values,upper_values):
        p=int(p);m=int(m)
        if u==b:total+=ceildiv(gu*p,S);continue
        aa=max(0,u*p-m);bb=max(0,m+em-b*max(0,p-ep))
        total+=ceildiv(gb*aa+gu*bb,(u-b)*S)
    return total


def verify(witness,report,terminal,block_size,D,H=16):
    wr=Path(witness).read_bytes();rr=Path(report).read_bytes();w=json.loads(wr);r=json.loads(rr)
    assert r['status']=='VERIFIED' and r['witness_sha256']==hashlib.sha256(wr).hexdigest()
    assert r['cutoff']==w['K'] and r['terminal_index']==w['start']+len(w['deltas'])
    n0=r['terminal_index'];assert terminal>n0 and block_size>=1 and 0<D<=L//2
    ps=primes_to(20*terminal);assert len(ps)>=terminal and ps[n0-1]==r['terminal_prime']>=w['K']
    state=r['stratified_state'];points=state['bins']+state['upper'];cost=[0]*len(points)
    envelope=GeometricEnvelope(H);blocks=0
    for lo in range(n0,terminal,block_size):
        hi=min(terminal,lo+block_size)
        for q in ps[lo:hi]:envelope.update(q,D)
        values=block_values(envelope,points,ps[lo],hi-lo,D)
        cost=[a+b for a,b in zip(cost,values)];blocks+=1
        if blocks%100==0:print('Certified block',blocks,'through prime index',hi,flush=True)
    nb=len(state['bins']);finite=chord_sum(state,cost[:nb],cost[nb:]);table=prepare(r)
    assert table is not None
    coefficients=[]
    for row in r['stratified_orders']:
        num,den=map(int,row['order'].split('/'))
        coefficient=fractional_tail_bound(ps[terminal-1],terminal,ONE,D,num,den)[0]
        coefficients.append(ceildiv(coefficient*envelope.growth(num,den).hi,S))
    remainder=sum(min(ceildiv(b*m,S) for b,m in zip(coefficients,cell)) for cell in zip(*table))
    overflow=[]
    for row in r['stratified_orders']:
        if row['moment_upper'] is None:continue
        num,den=map(int,row['order'].split('/'))
        outside=max(0,scaled(row['moment_upper'])-scaled(row['retained_moment_lower']))
        coefficient=fractional_tail_bound(ps[n0-1],n0,ONE,D,num,den)[0]
        overflow.append(ceildiv(coefficient*outside,S))
    assert overflow
    tail=finite+remainder+min(overflow);margin=scaled(r['mu_lower'])-tail
    assert margin>0,('block tail failed',fixed_text(margin))
    return dict(status='VERIFIED_BLOCK_TAIL_FROM_COMPLETED_REPORT',cutoff=w['K'],minimum_modulus_at_most=w['K']-1,
                base_witness=witness,base_report=report,base_witness_sha256=hashlib.sha256(wr).hexdigest(),
                base_report_sha256=hashlib.sha256(rr).hexdigest(),base_index=n0,terminal_index=terminal,
                terminal_prime=ps[terminal-1],block_size=block_size,blocks=blocks,exponent_truncation=H,
                future_delta_numerator=D,decimal_places=40,mu_lower=r['mu_lower'],
                retained_finite_tail_upper=fixed_text(finite),retained_infinite_tail_upper=fixed_text(remainder),
                omitted_initial_paths_tail_upper=fixed_text(min(overflow)),tail_loss_upper=fixed_text(tail),
                infinite_survivor_lower=fixed_text(margin),finite_computation_replayed=False,
                note='Conditional on the separately completed base certificate. One common future distortion; all geometric-envelope remainders included.')


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled.')
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('witness');ap.add_argument('report')
    ap.add_argument('--terminal',type=int,default=1000000);ap.add_argument('--block-size',type=int,default=1000)
    ap.add_argument('--delta',type=int,default=400000000000);ap.add_argument('--H',type=int,default=16)
    ap.add_argument('--output',required=True);a=ap.parse_args()
    result=verify(a.witness,a.report,a.terminal,a.block_size,a.delta,a.H)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
