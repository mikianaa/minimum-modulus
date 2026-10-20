"""Adversarial optimisation of the global class assignment at one stage.

Companion to fable_class_collision_mc.py.  The comparison value of the 14,484
certificate at stage p is E(S_p(D)-delta)_+ with S_p(D)=sum_{m|D}F_p(m); it
treats the p-adic classes of all active moduli as disjoint on every fibre.  A
real family picks one class per modulus m p^j globally.  Here the adversary
(tower residues, all moduli used) chooses the level-1 class ("colour" in Z_p)
of every modulus m p with m >= K/p to MAXIMISE the importance-sampled stage loss
on a training sample of fibres D, by coordinate ascent; the result is then
evaluated on an independent test sample (moduli unseen in training keep a
random colour).  The test value is an honest lower bound on what this adversary
can realise; the certificate's comparison value is the upper bound.  Deeper
levels (j0 >= 2) and the j > j0 part of each chain are counted collision-free,
which favours the adversary.

Also reported: E(min(S,1)-delta)_+ (pure truncation at alpha <= 1) to separate the
truncation effect from class collisions.

Usage: python code/fable_class_collision_opt.py --stages 101,199,307 --train 20000 --test 20000 --sweeps 4
Requires numpy (tmp/fable-venv).  Fable-owned file; see memory fable-loop-plan.
"""
import argparse, json, math, time
import numpy as np
from fable_class_collision_mc import load_params, j0, sample_stage


def build(p, K, divlists, weights, universe):
    """Per-sample arrays of level-1 moduli indices; per-modulus sample index."""
    cp = p / (p - 1)
    S = np.zeros(len(divlists)); free = np.zeros(len(divlists))
    sample_mods = []
    mod_samples = {}
    for n, divs in enumerate(divlists):
        mods = []
        s = 0.0; f = 0.0
        for m in divs:
            j = universe[m]
            s += cp * p ** (-j)
            if j == 1:
                mods.append(m); f += cp / p - 1.0 / p
                mod_samples.setdefault(m, []).append(n)
            else:
                f += cp * p ** (-j)
        S[n] = s; free[n] = f
        sample_mods.append(mods)
    return S, free, sample_mods, mod_samples


def run(p, K, delta, ntrain, ntest, theta, rng, sweeps, seed):
    dp = delta[p]
    t0 = time.time()
    wtr, dtr, nu = sample_stage(p, K, delta, ntrain, theta, rng)
    wte, dte, _ = sample_stage(p, K, delta, ntest, theta, rng)
    universe = {}
    for divs in dtr + dte:
        for m in divs:
            if m not in universe:
                universe[m] = j0(m, p, K)
    Str, ftr, smods_tr, msamp_tr = build(p, K, dtr, wtr, universe)
    Ste, fte, smods_te, msamp_te = build(p, K, dte, wte, universe)

    lvl1 = [m for m, j in universe.items() if j == 1]
    color = {m: int(rng.integers(p)) for m in lvl1}

    def phi(a):
        return np.maximum(np.minimum(a, 1.0) - dp, 0.0)

    # training state: counts[n, c] = number of active level-1 moduli of colour c in sample n
    counts = np.zeros((ntrain, p), dtype=np.int32)
    for n, mods in enumerate(smods_tr):
        for m in mods:
            counts[n, color[m]] += 1
    ncol = (counts > 0).sum(axis=1)

    def objective():
        return float(np.sum(wtr * phi(ncol / p + ftr))) / ntrain / (1 - dp)

    hist = [objective()]
    # moduli ordered by total training weight (heaviest first)
    order = sorted(msamp_tr, key=lambda m: -sum(wtr[i] for i in msamp_tr[m]))
    for sw in range(sweeps):
        moved = 0
        for m in order:
            idx = np.array(msamp_tr[m])
            c0 = color[m]
            w = wtr[idx]
            base_ncol = ncol[idx]
            # removing m from c0 reduces colour count where it was the only one
            lose = (counts[idx, c0] == 1).astype(np.int32)
            n_wo = base_ncol - lose
            # adding to colour c adds 1 where counts[.,c]==0 (after removal for c0)
            cnt = counts[idx].copy()
            cnt[:, c0] -= 1
            gain = (cnt == 0).astype(np.int32)  # shape (k, p)
            vals = (w[:, None] * phi((n_wo[:, None] + gain) / p + ftr[idx][:, None])).sum(axis=0)
            c1 = int(np.argmax(vals))
            if vals[c1] > vals[c0] + 1e-15:
                counts[idx, c0] -= 1
                counts[idx, c1] += 1
                ncol[idx] = n_wo + gain[:, c1]
                color[m] = c1
                moved += 1
        hist.append(objective())
        print(json.dumps(dict(p=p, sweep=sw + 1, moved=moved, train_loss=hist[-1])), flush=True)

    def evaluate(S, free, smods, w, n):
        a = np.array([len({color[m] for m in mods}) / p for mods in smods]) + free
        return dict(loss_S=float(np.sum(w * np.maximum(S - dp, 0))) / n / (1 - dp),
                    loss_truncS=float(np.sum(w * phi(S))) / n / (1 - dp),
                    loss_alpha=float(np.sum(w * phi(a))) / n / (1 - dp),
                    first_moment_S=float(np.sum(w * S)) / n,
                    first_moment_alpha=float(np.sum(w * np.minimum(a, 1))) / n)

    out = dict(p=p, delta=dp, train=ntrain, test=ntest, sweeps=sweeps, seed=seed,
               level1_moduli=len(lvl1), train_history=hist)
    out["train"] = evaluate(Str, ftr, smods_tr, wtr, ntrain)
    out["test"] = evaluate(Ste, fte, smods_te, wte, ntest)
    te = out["test"]
    out["gap_truncation_only"] = 1 - te["loss_truncS"] / te["loss_S"]
    out["gap_optimised_adversary"] = 1 - te["loss_alpha"] / te["loss_S"]
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--witness", default="data/sparse-certificate-14485.json")
    ap.add_argument("--stages", default="101,199,307")
    ap.add_argument("--train", type=int, default=20000)
    ap.add_argument("--test", type=int, default=20000)
    ap.add_argument("--sweeps", type=int, default=4)
    ap.add_argument("--theta", type=float, default=2.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--output")
    a = ap.parse_args()
    K, delta = load_params(a.witness)
    rng = np.random.default_rng(a.seed)
    results = []
    for p in [int(x) for x in a.stages.split(",")]:
        r = run(p, K, delta, a.train, a.test, a.theta, rng, a.sweeps, a.seed)
        results.append(r)
        print(json.dumps(r), flush=True)
    if a.output:
        json.dump(dict(note="Fable diagnostic (20 Oct 2026): adversarially optimised global class "
                            "assignment vs the certificate's comparison value; test values are lower "
                            "bounds on the tower adversary, not proved bounds",
                       witness=a.witness, results=results), open(a.output, "w"), indent=1)


if __name__ == "__main__":
    main()
