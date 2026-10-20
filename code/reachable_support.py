"""Integer support bounds for the directed graph of rounded products.

This file uses only integer arithmetic, including in exploratory searches.
"""
import bisect
from functools import lru_cache


@lru_cache(maxsize=16)
def reachable_upper_support(bins, max_multiplier=108):
    assert bins[0]==0 and bins[1]==1 and max_multiplier>=2
    upper=[0]*len(bins);upper[1]=1
    for j,b in enumerate(bins[1:],1):
        if not upper[j]:continue
        for h in range(2,min(max_multiplier,bins[-1]//b)+1):
            k=bisect.bisect_right(bins,b*h)-1
            assert k>j, 'Nontrivial transitions must strictly increase bins.'
            upper[k]=max(upper[k],h*upper[j])
    # Prefix queries require monotone upper endpoints. Empty bins can be
    # assigned any support containing their lower endpoint.
    for j,b in enumerate(bins):
        upper[j]=max(b,upper[j],upper[j-1] if j else 0)
    return tuple(upper)
