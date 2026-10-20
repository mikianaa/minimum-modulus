"""Downward-rounded product distribution, for exploratory calculations."""
import numpy as np
from functools import lru_cache


def bin_values(cap,resolution):
    values=[0,1]
    while values[-1]<cap:
        values.append(min(cap,values[-1]+max(1,values[-1]//resolution)))
    return values


class Binned:
    def __init__(self,cap,resolution=256):
        self.values=np.asarray(bin_values(cap,resolution),dtype=np.int64)
        self.maps={}

    def update(self,dist,q,nu):
        new=(1-nu/q)*dist
        weight=nu*(q-1)/q**2;h=2
        while weight>1e-30:
            if h not in self.maps:
                product=self.values*h
                keep=product<=self.values[-1]
                self.maps[h]=(keep,np.searchsorted(self.values,product[keep],side='right')-1)
            keep,indices=self.maps[h]
            new+=np.bincount(indices,weights=weight*dist[keep],minlength=len(dist))
            weight/=q;h+=1
        return new

    def update_pair(self,dist,moment,q,nu):
        new=(1-nu/q)*dist;newmoment=(1-nu/q)*moment
        weight=nu*(q-1)/q**2;h=2
        while weight>1e-30:
            if h not in self.maps:
                product=self.values*h;keep=product<=self.values[-1]
                self.maps[h]=(keep,np.searchsorted(self.values,product[keep],side='right')-1)
            keep,indices=self.maps[h]
            new+=np.bincount(indices,weights=weight*dist[keep],minlength=len(dist))
            newmoment+=np.bincount(indices,weights=h*weight*moment[keep],minlength=len(dist))
            weight/=q;h+=1
        return new,newmoment


@lru_cache(maxsize=8)
def upper_support(cap,resolution):
    values=np.asarray(bin_values(cap,resolution))
    # Each nontrivial update grows the rounded value by at least
    # 2R/(R+1); hence there are at most N such updates on a retained path.
    R=resolution;assert R>=2
    n=0;num=den=1
    while num<=int(values[-1])*den:n+=1;num*=2*R;den*=R+1
    rnum=(R+1)**(n-1);rden=R**(n-1)
    return np.asarray([int(b) if b<R else (int(b)*rnum+rden-1)//rden for b in values])


class MomentPrefixes:
    def __init__(self,values,dist,moment,resolution,support_mode='uniform'):
        self.upper=upper_support(int(values[-1]),resolution)
        if support_mode=='reachable':
            from reachable_support import reachable_upper_support
            self.upper=np.asarray(reachable_upper_support(tuple(map(int,values))))
        self.values=values
        width=self.upper-values
        a=np.divide(self.upper*dist-moment,width,out=np.zeros_like(dist),where=width>0)
        self.P=np.cumsum(dist);self.M=np.cumsum(moment)
        self.A=np.cumsum(a);self.B=np.cumsum(a*values)

    def complement(self,t):
        j=np.searchsorted(self.values,t,side='right')-1
        k=np.searchsorted(self.upper,t,side='right')-1
        return t*self.P[k]-self.M[k]+t*(self.A[j]-self.A[k])-(self.B[j]-self.B[k])
