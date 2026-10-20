"""Finite recursive residue caps with an analytic infinite-sum fallback.

All proof arithmetic is outward fixed point. F[n] bounds EVERY residue
modulo n, not the probability of a chosen nested comparison event.
Only the old multiplicative cap is used outside the finite table.
"""
import math
from verify_minimum_modulus import I, ONE, S, ceildiv, factorizations
from verify_stoploss import L


class RecursiveAtoms:
    def __init__(self, K, cap, density_degree=1):
        assert isinstance(K, int) and K > 2 and isinstance(cap, int) and cap >= 2
        assert density_degree in (1,2)
        self.K, self.cap, self.density_degree = K, cap, density_degree
        self.limit = max(cap, (K-1)//2)
        self.fs = factorizations(self.limit)
        self.largest = [0]+[f[-1][0] if f else 1 for f in self.fs[1:]]
        self.B = {1: ONE}
        self.F = {1: S}
        self.nu = {}; self.primes = []; self.C = ONE; self.T = ONE
        self.count = 0

    def coefficient(self, p, m):
        """h_p(m)=sum_{a>=1, m*p**a>=K} p**(-a), exactly rational."""
        power = p
        while m*power < self.K:
            power *= p
        return I.rat(p, power*(p-1))

    def unrestricted_lcm_sum(self, d):
        ans = self.C*self.B[d]
        for q, a in self.fs[d]:
            ans = ans*(I.rat(a+1)+I.rat(1,q-1))/(ONE+self.nu[q]*I.rat(1,q-1))
        return ans

    def H_upper(self, p, d):
        """Infinite restricted LCM sum, minus finite certified improvements."""
        total = (self.unrestricted_lcm_sum(d)*I.rat(1,p-1)).hi
        for m in range(1, (self.K-1)//p+1):
            if self.largest[m] >= p:
                continue
            b = self.B[d]*self.B[m]/self.B[math.gcd(d,m)]
            excluded = I.rat(1,p-1)-self.coefficient(p,m)
            total -= (excluded*b).lo
        for m in range(1, self.cap+1):
            if self.largest[m] >= p:
                continue
            n = math.lcm(d,m)
            if n > self.cap:
                continue
            gain = max(0, self.B[n].lo-self.F[n])
            total -= (self.coefficient(p,m)*I(gain)).lo
        assert total >= 0
        return total

    def H2_upper(self,p,d):
        # An unrestricted conditional second moment is still a valid bound
        # for the restricted one. No finite improvement is subtracted here.
        ans=self.T*self.B[d]*I.rat(1,(p-1)**2)
        for q,a in self.fs[d]:
            local=I.rat((a+1)**2)+I.rat(2*a+3,q-1)+I.rat(2,(q-1)**2)
            ans=ans*local/(ONE+self.nu[q]*I.rat(3*q-1,(q-1)**2))
        return ans.hi

    def update(self, p, D):
        assert 0 <= D <= L//2 and (not self.primes or p > self.primes[-1])
        v = I.rat(L,L-D)
        # New atoms have a unique earlier part; all old atoms are preserved.
        updates = {}
        for d in range(1, self.cap//p+1):
            if self.largest[d] >= p:
                continue
            old = self.F[d]
            a = (v*I(old)).hi
            if D:
                H=self.H_upper(p,d)
                a = min(a, old+(v*I(H)).hi)
                if self.density_degree==2:
                    a=min(a,old+H+(v*I(self.H2_upper(p,d))).hi)
            power = p
            while d*power <= self.cap:
                updates[d*power] = ceildiv(a,power)
                power *= p
        self.nu[p] = v; self.primes.append(p)
        self.C = self.C*(ONE+v*I.rat(1,p-1))
        self.T = self.T*(ONE+v*I.rat(3*p-1,(p-1)**2))
        for n in range(p, self.limit+1, p):
            if self.largest[n] != p:
                continue
            d = n
            while d % p == 0:
                d //= p
            self.B[n] = self.B[d]*v*I.rat(d,n)
        for n, value in updates.items():
            self.F[n] = min(self.B[n].hi,value)
            self.count += self.F[n] < self.B[n].lo

    def moment_gains(self, p, orders=(1,2)):
        """Lower gains for the restricted integer moments of alpha_p.

        For each n, the LCM coefficient is the divisor Mobius inverse of
        (sum_{d|n} h_p(d))**r. It is nonnegative by its tuple expansion.
        Subtract only gains from n in the finite table; the omitted sum is
        retained in the old moment bound, never dropped.
        """
        coeff = {r:{} for r in orders}; gains = {r:0 for r in orders}
        for n in sorted(self.F):
            divisors = [1]
            for q,a in self.fs[n]:
                old = divisors[:]
                for b in range(1,a+1):
                    divisors.extend(d*q**b for d in old)
            h = I(0)
            for d in divisors:
                h = h+self.coefficient(p,d)
            improvement = max(0,self.B[n].lo-self.F[n])
            for r in orders:
                v = ONE
                for _ in range(r): v = v*h
                for d in divisors:
                    if d != n: v = v-coeff[r][d]
                # Intersection weights are nonnegative; intersection of
                # our enclosure with [0,infinity) tightens it safely.
                assert v.hi >= 0
                v = I(max(0,v.lo),v.hi); coeff[r][n] = v
                gains[r] += improvement*v.lo//S
        return gains
