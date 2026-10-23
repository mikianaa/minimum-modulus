"""Fable 21 Oct 2026: exact small-instance adversaries for the stage loss (see notes/fable-notes-2026-10-21-upper.md). """
"""Exact cap-only adversary in a small instance.

Coordinates: primes Q with exponents G (space X = prod Z_{q^G}).  Stage prime p with
fibre Z_{p^gamma}.  Moduli m p^j (m | prod q^G, m p^j >= K) each get ONE residue r_m mod m
and ONE class c mod p^j.  For a point x, alpha(x) = |union of classes of active moduli| / p^gamma.
Measure P on X is any measure obeying BBMST's conditional atom caps
    P(X_j = t | X_<j) <= nu_j / q_j^G_j .
For a fixed family the worst P is an LP; the family is then improved by local search.
Compare with the convex comparison value (independent model, tower) and its truncation.
"""
import itertools, math, random, sys, json, time
import numpy as np
from scipy.optimize import linprog

def model_values(Q, G, nu, p, gamma, K, delta):
    # independent model: V_q in 0..G_q with P(V_q >= a) = nu_q/q^a (a>=1)
    dists = []
    for q, g, n in zip(Q, G, nu):
        pr = []
        for a in range(g + 1):
            up = n / q**a if a >= 1 else 1.0
            lo = n / q**(a + 1) if a + 1 <= g else 0.0
            pr.append(up - lo)
        dists.append(pr)
    def F(m):
        return sum(p**-j for j in range(1, gamma + 1) if m * p**j >= K)
    tot_untr = tot_tr = 0.0
    for es in itertools.product(*[range(g + 1) for g in G]):
        pr = math.prod(d[e] for d, e in zip(dists, es))
        S = 0.0
        for fs in itertools.product(*[range(e + 1) for e in es]):
            m = math.prod(q**f for q, f in zip(Q, fs))
            S += F(m)
        tot_untr += pr * max(S - delta, 0)
        tot_tr += pr * max(min(S, 1) - delta, 0)
    return tot_untr / (1 - delta), tot_tr / (1 - delta)

class Instance:
    def __init__(self, Q, G, nu, p, gamma, K, delta, seed=0):
        self.Q, self.G, self.nu, self.p, self.gamma, self.K, self.delta = Q, G, nu, p, gamma, K, delta
        self.mods = [q**g for q, g in zip(Q, G)]
        self.points = list(itertools.product(*[range(M) for M in self.mods]))
        self.index = {x: i for i, x in enumerate(self.points)}
        self.n = len(self.points)
        # moduli
        self.moduli = []
        for fs in itertools.product(*[range(g + 1) for g in G]):
            m = math.prod(q**f for q, f in zip(Q, fs))
            for j in range(1, gamma + 1):
                if m * p**j >= K:
                    self.moduli.append((m, j, fs))
        self.fibre = p**gamma
        self.rng = random.Random(seed)
        self.build_caps()

    def build_caps(self):
        # conditional caps: for coordinate k, prefix s (coords < k), value t:
        #   sum_{later} P(s,t,later) <= nu_k/mod_k * sum_{t',later} P(s,t',later)
        rows = []
        Q, mods = self.Q, self.mods
        r = len(Q)
        for k in range(r):
            prefixes = itertools.product(*[range(mods[i]) for i in range(k)])
            for s in prefixes:
                for t in range(mods[k]):
                    row = np.zeros(self.n)
                    for i, x in enumerate(self.points):
                        if x[:k] == s:
                            row[i] = (1.0 if x[k] == t else 0.0) - self.nu[k] / mods[k]
                    rows.append(row)
        self.A_ub = np.array(rows)
        self.b_ub = np.zeros(len(rows))

    def tower_family(self):
        # all residues 0, classes: assign distinct classes greedily by level
        fam = {}
        c1 = 0; c2 = 0
        for (m, j, fs) in self.moduli:
            if j == 1:
                fam[(m, j)] = (0, c1 % self.p); c1 += 1
            else:
                fam[(m, j)] = (0, c2 % self.fibre); c2 += 1
        return fam

    def alpha(self, fam):
        p, gamma, fib = self.p, self.gamma, self.fibre
        al = np.zeros(self.n)
        S = np.zeros(self.n)
        for i, x in enumerate(self.points):
            covered = set()
            s = 0.0
            for (m, j, fs) in self.moduli:
                r, c = fam[(m, j)]
                # x mod m: coordinates
                ok = True
                for q, f, xi in zip(self.Q, fs, x):
                    if f and (xi - r) % q**f != 0:
                        ok = False; break
                if ok:
                    s += p**-j
                    pj = p**j
                    for y in range(c % pj, fib, pj):
                        covered.add(y)
            al[i] = len(covered) / fib
            S[i] = s
        return al, S

    def lp(self, al):
        obj = -np.maximum(al - self.delta, 0)
        res = linprog(obj, A_ub=self.A_ub, b_ub=self.b_ub, A_eq=np.ones((1, self.n)), b_eq=[1.0],
                      bounds=(0, None), method='highs')
        assert res.status == 0, res.message
        return -res.fun / (1 - self.delta), res.x

    def value(self, fam):
        al, S = self.alpha(fam)
        v, P = self.lp(al)
        return v, P, al, S

    def local_search(self, fam, iters=2000, log=None):
        best, P, al, S = self.value(fam)
        keys = list(fam.keys())
        for it in range(iters):
            k = self.rng.choice(keys)
            m, j = k
            r, c = fam[k]
            new = dict(fam)
            if self.rng.random() < 0.5:
                new[k] = (self.rng.randrange(m), c)
            else:
                new[k] = (r, self.rng.randrange(self.p**j))
            v, P2, al2, S2 = self.value(new)
            if v >= best - 1e-12:
                if v > best + 1e-9 and log is not None:
                    log.append((it, v))
                fam, best, P, al, S = new, v, P2, al2, S2
        return fam, best, P, al, S

def run(Q, G, nu, p, gamma, K, delta, iters, seed):
    inst = Instance(Q, G, nu, p, gamma, K, delta, seed)
    untr, tr = model_values(Q, G, nu, p, gamma, K, delta)
    fam = inst.tower_family()
    v0, P0, al0, S0 = inst.value(fam)
    # value of tower family under independent model measure is not needed; LP sup is what matters
    log = []
    fam, best, P, al, S = inst.local_search(fam, iters, log)
    out = dict(Q=Q, G=G, nu=nu, p=p, gamma=gamma, K=K, delta=delta, points=inst.n,
               moduli=len(inst.moduli), model_untruncated=untr, model_truncated=tr,
               tower_lp=v0, search_lp=best, improvements=len(log),
               sat_mass_final=float(P[S > 1].sum()), max_S=float(S.max()))
    return out, inst, fam, P, al, S

if __name__ == '__main__':
    cfg = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    Q = cfg.get('Q', [2, 3, 5]); G = cfg.get('G', [3, 2, 1]); nu = cfg.get('nu', [1.25, 1.25, 1.25])
    p = cfg.get('p', 7); gamma = cfg.get('gamma', 2); K = cfg.get('K', 28); delta = cfg.get('delta', 0.25)
    iters = cfg.get('iters', 1500); seed = cfg.get('seed', 0)
    t = time.time()
    out, inst, fam, P, al, S = run(Q, G, nu, p, gamma, K, delta, iters, seed)
    out['seconds'] = time.time() - t
    print(json.dumps(out, indent=1))
