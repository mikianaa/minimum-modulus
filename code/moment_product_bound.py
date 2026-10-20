"""Outward integer bounds with true first moments inside downward bins."""
import bisect
from verify_stoploss import ProductBound,L
from verify_minimum_modulus import S,ceildiv


class MomentProductBound(ProductBound):
    def __init__(self,cap,resolution,support_mode='uniform'):
        super().__init__(cap,resolution)
        R=resolution;assert R>=2
        n=0;num=den=1
        while num<=cap*den:n+=1;num*=2*R;den*=R+1
        self.max_updates=n-1
        rnum=(R+1)**(n-1);rden=R**(n-1)
        self.upper=[b if b<R else ceildiv(b*rnum,rden) for b in self.bins]
        assert support_mode in ('uniform','reachable')
        self.support_mode=support_mode
        if support_mode=='reachable':
            from reachable_support import reachable_upper_support
            self.upper=list(reachable_upper_support(tuple(self.bins)))
        self.moment=[0]*len(self.bins);self.moment[1]=S

    def update(self,q,D):
        super().update(q,D)
        den0=q*(L-D);num0=den0-L
        new=[m*num0//den0 for m in self.moment]
        H=1;qH=q
        while qH<self.tail_cutoff:H+=1;qH*=q
        if self.support_mode=='reachable':assert H<=108
        power=q*q;num=L*(q-1)
        for h in range(2,H+1):
            den=(L-D)*power
            for j,target in self.maps[h]:
                if self.moment[j]:new[target]+=self.moment[j]*num*h//den
            power*=q
        self.moment=new

    def prefixes(self):
        P=[];M=[];A=[];B=[];p=m=a=b=0
        for lower,upper,u,v in zip(self.bins,self.upper,self.mass,self.moment):
            p+=u;m+=v
            coefficient=ceildiv(upper*u-v,upper-lower) if upper>lower else 0
            assert coefficient>=0
            a+=coefficient;b+=lower*coefficient
            P.append(p);M.append(m);A.append(a);B.append(b)
        return P,M,A,B

    def excess_upper(self,tlo,mean_upper,prefixes):
        tlo=min(tlo,self.bins[-1]*S);assert tlo>=0
        j=bisect.bisect_right(self.bins,tlo//S)-1
        k=bisect.bisect_right(self.upper,tlo//S)-1
        P,M,A,B=prefixes
        numerator=tlo*P[k]-S*M[k]+tlo*(A[j]-A[k])-S*(B[j]-B[k])
        comp=ceildiv(numerator,S);assert comp>=0
        bound=mean_upper-tlo+comp+ceildiv(tlo*self.error,S)
        assert bound>=0
        return bound
