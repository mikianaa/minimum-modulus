"""Floating stratified terminal cost and derivatives; discovery only."""
import numpy as np
from fractional_moment_tail import fractional_tail_bound
from verify_minimum_modulus import ONE,S


class StratifiedObjective:
    def __init__(self,lower,upper,P,k,D,orders,overflow_order):
        self.lower=np.asarray(lower,dtype=float);upper=np.asarray(upper,dtype=float)
        width=upper-self.lower;slopes=[];intercepts=[];self.excluded_orders=[]
        for num,den in orders:
            s=num/den
            try:b=fractional_tail_bound(P,k,ONE,D,num,den)[0]/S
            except AssertionError:
                # High orders may fail the analytic stopping conditions at
                # an earlier diagnostic checkpoint. The exact certificate
                # does not catch these failures or omit its listed orders.
                self.excluded_orders.append((num,den));continue
            a=np.divide(upper**s-self.lower**s,width,out=np.zeros_like(width),where=width>0)
            c=np.where(width>0,self.lower**s-a*self.lower,self.lower**s)
            slopes.append(b*a);intercepts.append(b*c)
        assert slopes
        self.slopes=np.asarray(slopes);self.intercepts=np.asarray(intercepts)
        num,den=overflow_order;self.order=num/den
        self.coefficient=fractional_tail_bound(P,k,ONE,D,num,den)[0]/S

    def evaluate(self,P,M,U):
        raw=self.intercepts*P+self.slopes*M
        values=np.maximum(0.,raw);choice=np.argmin(values,axis=0);ix=np.arange(len(P))
        AP=self.intercepts[choice,ix].copy();AM=self.slopes[choice,ix].copy()
        negative=raw[choice,ix]<0;AP[negative]=0;AM[negative]=0
        mean=np.divide(M,P,out=np.zeros_like(M),where=P>0);t=self.order
        use_jensen=mean>self.lower
        retained=np.maximum(mean,self.lower)**(t-1)*M
        outside=U-float(retained.sum());tail=float(values[choice,ix].sum())
        moment_cost=0.
        if outside>0:
            b=self.coefficient;tail+=b*outside;moment_cost=b*U
            AP+=np.where(use_jensen,b*(t-1)*mean**t,0.)
            AM-=b*np.where(use_jensen,t*mean**(t-1),self.lower**(t-1))
        return tail,AP,AM,moment_cost
