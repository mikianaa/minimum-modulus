"""Uniform-measure adversary: coordinates with nu=1 force P uniform, so only the family matters.
Exact evaluation of E_unif (alpha_F - delta)_+/(1-delta) and simulated annealing over residues/classes."""
import json, sys, math, time, random, itertools
import numpy as np
from fable_adversary_toy import model_values

class UInst:
    def __init__(self, Q, G, p, gamma, K, delta, seed):
        self.Q, self.G, self.p, self.gamma, self.K, self.delta = Q, G, p, gamma, K, delta
        self.mods = [q**g for q, g in zip(Q, G)]
        self.points = np.array(list(itertools.product(*[range(M) for M in self.mods])))
        self.n = len(self.points); self.fibre = p**gamma
        self.moduli = []
        for fs in itertools.product(*[range(g + 1) for g in G]):
            m = math.prod(q**f for q, f in zip(Q, fs))
            for j in range(1, gamma + 1):
                if m * p**j >= K: self.moduli.append((m, j, fs))
        self.rng = random.Random(seed)
        # precompute residue vectors: x mod m for each modulus (as integer via CRT-free: compare coordinates)
        self.coordmod = []
        for (m, j, fs) in self.moduli:
            self.coordmod.append([q**f for q, f in zip(Q, fs)])
        self.classmask = {}  # (j, c) -> boolean fibre mask
        for j in range(1, gamma + 1):
            pj = p**j
            for c in range(pj):
                mk = np.zeros(self.fibre, dtype=bool); mk[c::pj] = True; self.classmask[(j, c)] = mk

    def active(self, idx, r):
        cm = self.coordmod[idx]
        ok = np.ones(self.n, dtype=bool)
        for k, M in enumerate(cm):
            if M > 1: ok &= (self.points[:, k] % M) == (r % M)
        return ok

    def value(self, fam):
        # fam: list of (r, c) per modulus index
        cov = np.zeros((self.n, self.fibre), dtype=bool)
        for idx, (m, j, fs) in enumerate(self.moduli):
            r, c = fam[idx]
            act = self.active(idx, r)
            cov[act] |= self.classmask[(j, c)]
        al = cov.sum(1) / self.fibre
        return float(np.maximum(al - self.delta, 0).mean() / (1 - self.delta)), al

    def tower(self):
        fam = []; c1 = c2 = 0
        for (m, j, fs) in self.moduli:
            if j == 1: fam.append((0, c1 % self.p)); c1 += 1
            else: fam.append((0, c2 % self.fibre)); c2 += 1
        return fam

    def anneal(self, fam, iters, T0, T1):
        cur, _ = self.value(fam); best = cur; bestfam = list(fam); log = []
        nm = len(self.moduli)
        for it in range(iters):
            T = T0 * (T1 / T0) ** (it / max(1, iters - 1))
            idx = self.rng.randrange(nm); m, j, fs = self.moduli[idx]; r, c = fam[idx]
            new = list(fam); u = self.rng.random()
            if u < 0.45: new[idx] = (self.rng.randrange(m), c)
            elif u < 0.9: new[idx] = (r, self.rng.randrange(self.p**j))
            else:
                i2 = self.rng.randrange(nm)
                if self.moduli[i2][1] == j and i2 != idx:
                    r2, c2 = fam[i2]; new[idx] = (r, c2); new[i2] = (r2, c)
            v, _ = self.value(new)
            if v >= cur or self.rng.random() < math.exp((v - cur) / T):
                fam, cur = new, v
                if v > best + 1e-12: best, bestfam = v, list(fam); log.append((it, v))
        return bestfam, best, log

if __name__ == '__main__':
    cfg = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    Q = cfg.get('Q', [2, 3, 5]); G = cfg.get('G', [3, 2, 1]); p = cfg.get('p', 7); gamma = cfg.get('gamma', 2)
    K = cfg.get('K', 28); delta = cfg.get('delta', 0.25); iters = cfg.get('iters', 20000); seed = cfg.get('seed', 0)
    T0 = cfg.get('T0', 0.01); T1 = cfg.get('T1', 1e-5)
    t = time.time()
    inst = UInst(Q, G, p, gamma, K, delta, seed)
    untr, tr = model_values(Q, G, [1.0] * len(Q), p, gamma, K, delta)
    fam = inst.tower(); v0, al0 = inst.value(fam)
    fam, best, log = inst.anneal(fam, iters, T0, T1)
    vb, al = inst.value(fam)
    out = dict(cfg=cfg, points=inst.n, moduli=len(inst.moduli), model_untruncated=untr, model_truncated=tr,
               tower_value=v0, anneal_best=best, ratio_best_to_truncated=best / tr, ratio_best_to_untruncated=best / untr,
               improvements=len(log), last_improvement=log[-1][0] if log else None,
               frac_saturated_points=float((al >= 0.999).mean()), seconds=time.time() - t,
               family=[[m, j, r, c] for (m, j, fs), (r, c) in zip(inst.moduli, fam)])
    json.dump(out, open(f"{cfg.get('name', 'urun')}.json", 'w'), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'family'}, indent=1))
