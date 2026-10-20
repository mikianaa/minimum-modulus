"""Can a non-tower adversary recover the truncation loss alpha <= 1?

In the comparison model of the 14,484 certificate the stage-p loss is
E(S_p(D)-delta)_+ with S_p(D)=sum_{m|D}F_p(m), attained by the "tower" family
(every residue a_d compatible with one point x*).  Real fibres satisfy
alpha <= 1, and for the tower family the truncated value E(min(S,1)-delta)_+ is
20-25% smaller.  The nested-event lemma does not allow truncation because a
family with incompatible residues can "water-fill": moduli that would overflow
a saturated tower fibre are moved to sibling cones that still have headroom.

This script measures how much of the truncation loss a simple branching
adversary recovers.  Fibres x are sampled (importance sampling) from the
product measure around x* = 0 with the certificate's caps: coordinate q has
agreement depth E_q (P(E_q >= a) = nu_q/q^a) and first disagreeing digit
d_q uniform in 1..q-1.  Each modulus m (all m | L with m p^j >= K, weight
F_p(m)) is placed greedily, heaviest first, either in its tower cone
{x = 0 mod m} or in one of the top-level sibling cones
{x = s q^(b-1) mod q^b, x = 0 mod m/q^b} for a prime q | m with exponent b,
whichever adds the most truncated excess sum_x w_x [phi(S_x+F) - phi(S_x)],
phi(s) = (min(s,1)-delta)_+.  Collisions of p-adic classes are ignored
(alpha = S), so the result isolates the truncation effect.

Reported on the same sample: comparison value E(S_tower-delta)_+, truncated
tower value, and the branching adversary's truncated value.  The greedy is fitted
to the sample (optimistic for the adversary on moduli seen in few fibres).

Usage: python code/fable_branching_adversary.py --stages 101,199 --samples 30000
Fable-owned file; see memory fable-loop-plan.
"""
import argparse, json, math, time
import numpy as np
from fable_class_collision_mc import load_params, j0


def sample_fibres(p, K, delta, n, theta_max, rng):
    qs = [q for q in sorted(delta) if q < p]
    nu = {q: 1 / (1 - delta[q]) for q in qs}
    theta = {q: min(theta_max, 0.9 * q) for q in qs}
    W = np.empty(n); E = []; Dg = []; divl = []
    for i in range(n):
        w = 1.0; e = {}; dg = {}; divs = [1]
        for q in qs:
            th = theta[q]; u = rng.random()
            if u >= th / q:
                a = 0
            else:
                a = 1
                while u < th / q ** (a + 1): a += 1
            nuq = nu[q]
            if a == 0:
                w *= (1 - nuq / q) / (1 - th / q)
            else:
                w *= nuq / th
                divs = divs + [d * q ** k for d in divs for k in range(1, a + 1)]
            e[q] = a
            dg[q] = int(rng.integers(1, q))  # first disagreeing digit (position a+1)
        W[i] = w; E.append(e); Dg.append(dg); divl.append(divs)
    return W, E, Dg, divl, qs


def factor(m, qs):
    f = {}
    for q in qs:
        if m % q == 0:
            b = 0
            while m % q == 0:
                m //= q; b += 1
            f[q] = b
    return f


def index_sample(W, E, Dg, divl, qs, p, K):
    n = len(W); cp = p / (p - 1)
    idx = {}
    for i, divs in enumerate(divl):
        for m in divs:
            idx.setdefault(m, []).append(i)
    idx = {m: np.array(v) for m, v in idx.items()}
    F = {m: cp * p ** (-j0(m, p, K)) for m in idx}
    Earr = {q: np.array([e[q] for e in E]) for q in qs}
    Darr = {q: np.array([d[q] for d in Dg]) for q in qs}
    return idx, F, Earr, Darr


def cone_of(m, choice, idx, Earr, Darr, qs, nall):
    """Sample indices of the cone chosen for m and the activity factor.

    'tower': x = 0 mod m.  (q, a1, s): x agrees with x* in coordinate q exactly to
    depth a1 (< b = exponent of q in m), has digit s at position a1+1 and x*'s digits
    at positions a1+2..b; the other coordinates are tower.  The deeper digits are not
    sampled: the modulus is counted with fractional activity q^-(b-a1-1), the
    conditional probability that they match (exact for a1 = b-1)."""
    if choice == "tower":
        return idx.get(m), 1.0
    q, a1, s = choice
    b = factor(m, qs)[q]
    rest = m // q ** b
    parent = idx.get(rest) if rest > 1 else np.arange(nall)
    if parent is None:
        return None, 0.0
    sel = parent[Earr[q][parent] == a1]
    cone = sel[Darr[q][sel] == s] if len(sel) else sel
    return cone, float(q) ** (-(b - a1 - 1))


def run(p, K, delta, ntrain, ntest, theta, rng, min_support=40):
    dp = delta[p]; t0 = time.time()
    W, E, Dg, divl, qs = sample_fibres(p, K, delta, ntrain, theta, rng)
    idx, F, Earr, Darr = index_sample(W, E, Dg, divl, qs, p, K)
    n = ntrain

    def phi(s):
        return np.maximum(np.minimum(s, 1.0) - dp, 0.0)

    def values(W, idx, F, S):
        n = len(W)
        S_tower = np.zeros(n)
        for m, arr in idx.items():
            S_tower[arr] += F[m]
        return dict(loss_comparison=float(np.sum(W * np.maximum(S_tower - dp, 0))) / n / (1 - dp),
                    loss_tower_truncated=float(np.sum(W * phi(S_tower))) / n / (1 - dp),
                    loss_branching_truncated=float(np.sum(W * phi(S))) / n / (1 - dp))

    # coordinate ascent on the training sample, starting from the all-tower family:
    # each modulus (heaviest first) is moved to the cone that maximises the truncated
    # excess given the current position of all other moduli; two sweeps
    order = sorted(idx, key=lambda m: -float(W[idx[m]].sum()))
    S = np.zeros(n)
    for m, arr in idx.items():
        S[arr] += F[m]
    choice = {m: "tower" for m in idx}
    cones = {}
    for sweep in range(2):
        for m in order:
            f = F[m]
            tower = idx[m]
            if len(tower) < min_support:
                continue
            cur_cone, cur_fac = cones.get(m, (tower, 1.0))
            S[cur_cone] -= f * cur_fac
            best_gain = float(np.sum(W[tower] * (phi(S[tower] + f) - phi(S[tower]))))
            best = ("tower", tower, 1.0)
            for q, b in factor(m, qs).items():
                rest = m // q ** b
                parent = idx.get(rest) if rest > 1 else np.arange(n)
                if parent is None:
                    continue
                for a1 in range(0, b):
                    sel = parent[Earr[q][parent] == a1]
                    if len(sel) < min_support:
                        continue
                    fac = float(q) ** (-(b - a1 - 1))
                    digs = Darr[q][sel]
                    for s in range(1, q):
                        cone = sel[digs == s]
                        if len(cone) < min_support:
                            continue
                        g = float(np.sum(W[cone] * (phi(S[cone] + f * fac) - phi(S[cone]))))
                        if g > best_gain + 1e-18:
                            best_gain = g; best = ((q, a1, s), cone, fac)
            S[best[1]] += f * best[2]
            choice[m] = best[0]; cones[m] = (best[1], best[2])
    n_branch = sum(1 for m in choice if choice[m] != "tower")
    w_total = sum(float(W[idx[m]].sum()) for m in idx)
    w_branch = sum(float(W[idx[m]].sum()) for m in choice if choice[m] != "tower")
    train = values(W, idx, F, S)

    # honest evaluation on a fresh sample; moduli unseen in training stay in the tower
    W2, E2, Dg2, divl2, _ = sample_fibres(p, K, delta, ntest, theta, rng)
    idx2, F2, Earr2, Darr2 = index_sample(W2, E2, Dg2, divl2, qs, p, K)
    S2 = np.zeros(ntest)
    for m in idx2:
        cone, fac = cone_of(m, choice.get(m, "tower"), idx2, Earr2, Darr2, qs, ntest)
        if cone is not None and len(cone):
            S2[cone] += F2[m] * fac
    # moduli chosen in training but absent from every test tower cone can still be
    # active on test sibling cones
    for m, ch in choice.items():
        if ch != "tower" and m not in idx2:
            cone, fac = cone_of(m, ch, idx2, Earr2, Darr2, qs, ntest)
            if cone is not None and len(cone):
                S2[cone] += F[m] * fac
    test = values(W2, idx2, F2, S2)
    for d in (train, test):
        d["gap_truncation"] = 1 - d["loss_tower_truncated"] / d["loss_comparison"]
        d["gap_branching"] = 1 - d["loss_branching_truncated"] / d["loss_comparison"]
        span = d["loss_comparison"] - d["loss_tower_truncated"]
        d["recovered_fraction"] = (d["loss_branching_truncated"] - d["loss_tower_truncated"]) / span if span > 0 else None
    out = dict(p=p, delta=dp, train_samples=ntrain, test_samples=ntest, moduli=len(idx), branched_moduli=n_branch,
               branched_weight_fraction=w_branch / w_total if w_total else None,
               train=train, test=test, seconds=round(time.time() - t0, 1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--witness", default="data/sparse-certificate-14485.json")
    ap.add_argument("--stages", default="101,199")
    ap.add_argument("--samples", type=int, default=30000)
    ap.add_argument("--test", type=int, default=20000)
    ap.add_argument("--min-support", type=int, default=40)
    ap.add_argument("--theta", type=float, default=2.5)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--output")
    a = ap.parse_args()
    K, delta = load_params(a.witness)
    rng = np.random.default_rng(a.seed)
    results = []
    for p in [int(x) for x in a.stages.split(",")]:
        r = run(p, K, delta, a.samples, a.test, a.theta, rng, a.min_support)
        results.append(r); print(json.dumps(r), flush=True)
    if a.output:
        json.dump(dict(note="Fable diagnostic (20 Oct 2026): greedy top-level branching adversary vs tower; "
                            "truncated excess; sample-fitted, not a proved bound",
                       witness=a.witness, results=results), open(a.output, "w"), indent=1)


if __name__ == "__main__":
    main()
