"""Float search recurrence checks against independent exact/sparse laws."""
from itertools import product
from fractions import Fraction as F
from math import prod
from joint_product_experiment import product_distribution,reference_product_distribution,conditioned_loss,complementary_values
import numpy as np
from conditioned_complement_bound import loss_upper
from verify_sparse_conditioning import DynamicSmallBound
from sparse_conditioning_experiment import DynamicSparseStates
from bernoulli_conditioning_experiment import JointStates
from verify_minimum_modulus import I,S
from verify_stoploss import L


if __name__=='__main__':
    count=0
    for H,cs in product((1,3,15,60,160),((0,0,0),(1,0,2),(2,3,1))):
        ps=[5,7,11];nu=[1.2,1.5,1.8]
        a=product_distribution(ps,nu,cs,H);b=reference_product_distribution(ps,nu,cs,H)
        assert set(a)==set(b)
        assert max(abs(a[k]-b[k]) for k in a)<1e-14
        ts=np.asarray([1.,2.,5.,20.]);bs=np.asarray([.5,1.,2.,7.]);A=H/3
        for subtract in (False,True):
            got=complementary_values(a,A,ts,bs,.9,subtract)
            ref=np.asarray([sum(m*(max(0.,A+bb+.9*j-t*z)-(max(0.,A-t*z) if subtract else 0.)) for (z,j),m in b.items()) for t,bb in zip(ts,bs)])
            assert np.max(np.abs(got-ref))<1e-12
        count+=1
    for K,D in product((70,100,200),(0,L//7,2*L//5)):
        ps=[2,3,5,7,11];p=11;N=400;params={2:L//5,3:L//4,5:L//3,7:L//5}
        sm=DynamicSmallBound(ps[:2],N,(K-1)//p+1)
        sf=DynamicSparseStates(K,ps,0,2,N)
        for q in ps[:2]:sm.update(q,params[q]);sf.update(q,L/(L-params[q]))
        sf.prob=sf.normalizer*sf.coeff
        joint=JointStates(sf,K,5,5)
        # Any common M1 suffices for this recurrence comparison.
        m1=F(2,3)
        exact,_=loss_upper(sm,K,p,D,ps[2:4],params,5,N,I.rat(2,3))
        got=conditioned_loss(joint,p,D/L,ps[2:4],[L/(L-params[q]) for q in ps[2:4]],ps[:2],[L/(L-params[q]) for q in ps[:2]],N,float(m1))
        assert abs(got-exact/S)<1e-12,(K,D,got,exact/S)
        count+=1
    print('PASS:',count,'dense product/reference and float/integer complement comparisons.')
