"""Monte Carlo diagnostic: how much of the stage loss in the 14,484 certificate
is an artefact of letting every fibre use disjoint p-adic classes?

At stage p a modulus d = m p^j (m | D, j >= j0(m) = min{j: m p^j >= K}) covers
ONE class modulo p^j of the new coordinate, and that class is the same on every
fibre where m is active.  The certificate's comparison value E(S_p(D)-delta)_+
with S_p(D) = sum_{m|D} F_p(m) implicitly lets the classes of all active moduli
be pairwise disjoint on every fibre simultaneously.  A real family must choose
one class per modulus globally, so on heavy fibres classes collide and the
covered fraction alpha(D) is strictly less than S_p(D).

This script measures the loss realised by the comparison's own extremal
adversary (tower residues, all moduli used) under explicit global class
assignments:
  random  - each modulus gets a uniformly random class of its level;
  greedy  - moduli in decreasing order of activity weight nu(m)/m pick the
            class minimising a pairwise conflict proxy nu(lcm)/lcm with the
            moduli already placed there (only the heaviest --greedy-top moduli
            are optimised, the rest are random).
Levels: j = j0(m) is placed explicitly for j0 in {1,2}; everything deeper (weight
<= p^-3 per modulus) and the j > j0 part of each chain is counted collision-free,
which favours the adversary, so the measured gap is conservative.

Both alpha and S are estimated from one importance-sampled stream (proposal
tails theta/q^a) so their difference has low variance.
Usage: python code/fable_class_collision_mc.py --stages 101,199,307 --samples 50000
Requires numpy (tmp/fable-venv).  Fable-owned file; see memory fable-loop-plan.
"""
import argparse, json, math, time
from math import gcd
import numpy as np

L = 10 ** 12


def primes_to(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def load_params(witness):
    w = json.load(open(witness))
    K = w["K"]; ds = w["deltas"]
    ps = primes_to(max(K, 20 * len(ds)))
    while len(ps) < len(ds): ps = primes_to(2 * len(ps) * 20)
    delta = {ps[i]: ds[i] / L for i in range(len(ds))}
    return K, delta


def j0(m, p, K):
    j = 1; v = m * p
    while v < K: v *= p; j += 1
    return j


def sample_stage(p, K, delta, nsamples, theta_max, rng):
    """Return (weights, list of active-divisor lists, nu_of_m function data)."""
    qs = [q for q in sorted(delta) if q < p]
    nu = {q: 1 / (1 - delta[q]) for q in qs}
    theta = {q: min(theta_max, 0.9 * q) for q in qs}
    weights = np.empty(nsamples)
    divlists = []
    for n in range(nsamples):
        divs = [1]
        w = 1.0
        for q in qs:
            th = theta[q]
            u = rng.random()
            if u >= th / q:
                e = 0
            else:
                e = 1
                while u < th / q ** (e + 1): e += 1
            nuq = nu[q]
            if e == 0:
                w *= (1 - nuq / q) / (1 - th / q)
            else:
                w *= nuq / th
                divs = divs + [d * q ** k for d in divs for k in range(1, e + 1)]
        weights[n] = w
        divlists.append(divs)
    return weights, divlists, nu


def nu_of(m, nu):
    r = 1.0
    for q, v in nu.items():
        if m % q == 0:
            r *= v
    return r


def run(p, K, delta, nsamples, theta_max, rng, greedy_top):
    dp = delta[p]
    t0 = time.time()
    weights, divlists, nu = sample_stage(p, K, delta, nsamples, theta_max, rng)
    # universe of moduli m with their level j0 and activity weight
    universe = {}
    for divs in divlists:
        for m in divs:
            if m not in universe:
                universe[m] = j0(m, p, K)
    cp = p / (p - 1)
    # chain weight F_p(m) = p^{-j0} * p/(p-1); explicit part p^{-j0}; rest collision-free
    lvl1 = [m for m, j in universe.items() if j == 1]
    lvl2 = [m for m, j in universe.items() if j == 2]
    act_w = {m: nu_of(m, nu) / m for m in universe}

    # --- random placement
    rnd_color = {m: int(rng.integers(p)) for m in lvl1}
    rnd_node2 = {m: int(rng.integers(p * p)) for m in lvl2}

    # --- greedy placement (level 1 colours) for the heaviest moduli
    top = sorted(lvl1, key=lambda m: -act_w[m])[:greedy_top]
    members = [[] for _ in range(p)]
    grd_color = {}
    for m in top:
        best, bestc = None, None
        for c in range(p):
            cost = 0.0
            for m2 in members[c]:
                l = m * m2 // gcd(m, m2)
                cost += nu_of(l, nu) / l
            if best is None or cost < best:
                best, bestc = cost, c
        grd_color[m] = bestc
        members[bestc].append(m)
    for m in lvl1:
        if m not in grd_color:
            grd_color[m] = rnd_color[m]
    # level-2 nodes: greedy over the heaviest too (cheap: avoid nodes under heavy colours)
    grd_node2 = dict(rnd_node2)

    def alpha(divs, color, node2):
        cols = set(); nodes = set(); free = 0.0
        for m in divs:
            j = universe[m]
            if j == 1:
                cols.add(color[m]); free += cp / p - 1.0 / p
            elif j == 2:
                nodes.add(node2[m]); free += cp / p ** 2 - 1.0 / p ** 2
            else:
                free += cp * p ** (-j)
        a = len(cols) / p
        for nd in nodes:
            if (nd // p) not in cols:
                a += 1.0 / p ** 2
        return min(1.0, a + free)

    acc = {k: 0.0 for k in ("S", "rnd", "grd", "S2", "rnd2", "grd2", "m1S", "m1rnd", "m1grd")}
    for w, divs in zip(weights, divlists):
        S = sum(cp * p ** (-universe[m]) for m in divs)
        ar = alpha(divs, rnd_color, rnd_node2)
        ag = alpha(divs, grd_color, grd_node2)
        for key, val in (("S", S), ("rnd", ar), ("grd", ag)):
            ex = max(val - dp, 0.0)
            acc[key] += w * ex; acc[key + "2"] += (w * ex) ** 2; acc["m1" + key] += w * val
    n = nsamples
    out = dict(p=p, delta=dp, samples=n, universe=len(universe), level1=len(lvl1), level2=len(lvl2),
               greedy_top=min(greedy_top, len(lvl1)))
    for key in ("S", "rnd", "grd"):
        mean = acc[key] / n
        var = acc[key + "2"] / n - mean ** 2
        out["loss_" + key] = mean / (1 - dp)
        out["se_" + key] = math.sqrt(max(var, 0) / n) / (1 - dp)
        out["first_moment_" + key] = acc["m1" + key] / n
    out["gap_random"] = 1 - out["loss_rnd"] / out["loss_S"]
    out["gap_greedy"] = 1 - out["loss_grd"] / out["loss_S"]
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--witness", default="data/sparse-certificate-14485.json")
    ap.add_argument("--stages", default="101,199,307")
    ap.add_argument("--samples", type=int, default=50000)
    ap.add_argument("--theta", type=float, default=2.5)
    ap.add_argument("--greedy-top", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--output")
    a = ap.parse_args()
    K, delta = load_params(a.witness)
    rng = np.random.default_rng(a.seed)
    results = []
    for p in [int(x) for x in a.stages.split(",")]:
        r = run(p, K, delta, a.samples, a.theta, rng, a.greedy_top)
        results.append(r)
        print(json.dumps(r), flush=True)
    if a.output:
        json.dump(dict(note="Fable diagnostic (20 Oct 2026): stage loss realised by the tower adversary "
                            "under global class assignments vs the certificate's comparison value; "
                            "the assignments are heuristic lower bounds on the adversary, not proved bounds",
                       witness=a.witness, results=results), open(a.output, "w"), indent=1)


if __name__ == "__main__":
    main()
