"""toy.py + simulated annealing, warm restarts and measure diagnostics."""
import json, sys, math, time, random
import numpy as np
from fable_adversary_toy import Instance, model_values

def anneal(inst, fam, iters, T0, T1, log):
    best, P, al, S = inst.value(fam)
    cur = best; bestfam = dict(fam); bestP = P
    keys = list(fam.keys())
    for it in range(iters):
        T = T0 * (T1 / T0) ** (it / max(1, iters - 1))
        k = inst.rng.choice(keys); m, j = k; r, c = fam[k]
        new = dict(fam)
        u = inst.rng.random()
        if u < 0.45: new[k] = (inst.rng.randrange(m), c)
        elif u < 0.9: new[k] = (r, inst.rng.randrange(inst.p**j))
        else:  # swap classes of two moduli at same level
            k2 = inst.rng.choice(keys)
            if k2[1] == j and k2 != k:
                r2, c2 = fam[k2]; new[k] = (r, c2); new[k2] = (r2, c)
        v, P2, al2, S2 = inst.value(new)
        if v >= cur or inst.rng.random() < math.exp((v - cur) / T):
            fam, cur, P = new, v, P2
            if v > best + 1e-10:
                best, bestfam, bestP = v, dict(fam), P2; log.append((it, v))
    return bestfam, best, bestP

def diagnostics(inst, P):
    # per coordinate: distribution of conditional atom / cap
    out = {}
    for k in range(len(inst.Q)):
        cap = inst.nu[k] / inst.mods[k]
        prefix_mass = {}; atom = {}
        for i, x in enumerate(inst.points):
            s = x[:k]; prefix_mass[s] = prefix_mass.get(s, 0) + P[i]
            atom[(s, x[k])] = atom.get((s, x[k]), 0) + P[i]
        ratios = []
        for (s, t), v in atom.items():
            if prefix_mass[s] > 1e-9:
                ratios.append(v / prefix_mass[s] / cap)
        ratios = np.array(ratios)
        out[inst.Q[k]] = dict(at_cap=float((ratios > 0.999).mean()), zero=float((ratios < 1e-6).mean()),
                              middle=float(((ratios >= 1e-6) & (ratios <= 0.999)).mean()))
    out['support_points'] = int((P > 1e-9).sum())
    return out

if __name__ == '__main__':
    cfg = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    Q = cfg.get('Q', [2, 3, 5]); G = cfg.get('G', [3, 2, 1]); nu = cfg.get('nu', [1.25, 1.25, 1.25])
    p = cfg.get('p', 7); gamma = cfg.get('gamma', 2); K = cfg.get('K', 28); delta = cfg.get('delta', 0.25)
    iters = cfg.get('iters', 6000); seed = cfg.get('seed', 0); T0 = cfg.get('T0', 0.01); T1 = cfg.get('T1', 1e-5)
    t = time.time()
    inst = Instance(Q, G, nu, p, gamma, K, delta, seed)
    untr, tr = model_values(Q, G, nu, p, gamma, K, delta)
    fam = inst.tower_family(); log = []
    fam, best, P = anneal(inst, fam, iters, T0, T1, log)
    al, S = inst.alpha(fam)
    out = dict(cfg=cfg, model_untruncated=untr, model_truncated=tr, anneal_best=best,
               ratio_to_untruncated=best / untr, improvements=len(log), last_improvement=log[-1][0] if log else None,
               sat_mass=float(P[S > 1].sum()), diagnostics=diagnostics(inst, P), seconds=time.time() - t,
               family={f"{m}*{p}^{j}": [r, c] for (m, j), (r, c) in fam.items()})
    name = cfg.get('name', 'run')
    json.dump(out, open(f'{name}.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'family'}, indent=1))
