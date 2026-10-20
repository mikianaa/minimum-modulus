"""Rule-based global class assignment for the tower adversary (honest lower bound).

fable_class_collision_opt.py optimises colours on a finite training sample and is
therefore optimistic for the adversary on the fibres it saw and pessimistic on
unseen moduli.  Here the colour of EVERY p-smooth modulus m <= M (level j0 = 1,
i.e. m >= K/p) is fixed by a deterministic greedy rule, so the evaluation on a
fresh importance sample is an honest estimate of what this adversary realises:

  process m in increasing order; for each colour c accumulate the conflict
  weight sum nu(l)/l over already-coloured m' of colour c with l = lcm(m, m')
  <= k0 * m (i.e. m' | m k for some k <= k0; these are the pairs that co-occur
  on the heaviest fibres containing m); choose the colour of least weight.

Moduli above M (rare on heavy fibres) get a hash colour.  Reported: the stage
loss with this colouring vs the certificate's comparison value E(S-delta)_+ and
the truncated value E(min(S,1)-delta)_+.

Usage: python code/fable_class_collision_rule.py --stages 101,199,307 --M 2000000 --k0 24 --samples 20000
Fable-owned file; see memory fable-loop-plan.
"""
import argparse, json, math, time
from math import gcd
import numpy as np
from fable_class_collision_mc import load_params, j0, sample_stage, nu_of


def smooth_numbers(limit, primes):
    out = [1]
    for q in primes:
        new = []
        for v in out:
            w = v * q
            while w <= limit:
                new.append(w); w *= q
        out += new
    return sorted(out)


def divisors_from_factor(n, primes):
    divs = [1]
    for q in primes:
        if n % q == 0:
            e = 0
            while n % q == 0:
                n //= q; e += 1
            divs = [d * q ** k for d in divs for k in range(e + 1)]
    return divs


def build_rule(p, K, nu, M, k0, rng):
    qs = sorted(nu)
    lo = (K + p - 1) // p  # level-1 moduli: m >= K/p
    ms = [m for m in smooth_numbers(M, qs) if m >= lo]
    color = {}
    nucache = {}

    def nuw(l):
        v = nucache.get(l)
        if v is None:
            v = nu_of(l, nu) / l; nucache[l] = v
        return v

    cost = np.zeros(p)
    for m in ms:
        cost[:] = 0.0
        seen = set()
        for k in range(1, k0 + 1):
            for m2 in divisors_from_factor(m * k, qs):
                if m2 >= lo and m2 < m and m2 not in seen:
                    seen.add(m2)
                    c = color.get(m2)
                    if c is not None:
                        l = m * m2 // gcd(m, m2)
                        cost[c] += nuw(l)
        # ties (in particular all-zero cost) are broken uniformly at random
        best = np.flatnonzero(cost <= cost.min() * (1 + 1e-9) + 1e-300)
        color[m] = int(best[rng.integers(len(best))])
    return color, len(ms)


def run(p, K, delta, nsamples, theta, rng, M, k0):
    dp = delta[p]
    t0 = time.time()
    w, divlists, nu = sample_stage(p, K, delta, nsamples, theta, rng)
    color, nm = build_rule(p, K, nu, M, k0, rng)
    t1 = time.time()
    cp = p / (p - 1)
    accS = accT = accA = 0.0
    for wt, divs in zip(w, divlists):
        S = 0.0; free = 0.0; cols = set()
        for m in divs:
            j = j0(m, p, K)
            S += cp * p ** (-j)
            if j == 1:
                c = color.get(m)
                if c is None:
                    c = hash(m) % p
                cols.add(c); free += 1.0 / (p * (p - 1))
            else:
                free += cp * p ** (-j)
        a = min(1.0, len(cols) / p + free)
        accS += wt * max(S - dp, 0); accT += wt * max(min(S, 1) - dp, 0); accA += wt * max(a - dp, 0)
    n = nsamples
    out = dict(p=p, delta=dp, samples=n, M=M, k0=k0, rule_moduli=nm,
               loss_S=accS / n / (1 - dp), loss_truncS=accT / n / (1 - dp), loss_rule=accA / n / (1 - dp))
    out["gap_truncation_only"] = 1 - out["loss_truncS"] / out["loss_S"]
    out["gap_rule_adversary"] = 1 - out["loss_rule"] / out["loss_S"]
    out["collision_gap_vs_truncated"] = 1 - out["loss_rule"] / out["loss_truncS"]
    out["seconds_rule"] = round(t1 - t0, 1); out["seconds"] = round(time.time() - t0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--witness", default="data/sparse-certificate-14485.json")
    ap.add_argument("--stages", default="101,199,307")
    ap.add_argument("--samples", type=int, default=20000)
    ap.add_argument("--M", type=int, default=2000000)
    ap.add_argument("--k0", type=int, default=24)
    ap.add_argument("--theta", type=float, default=2.5)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--output")
    a = ap.parse_args()
    K, delta = load_params(a.witness)
    rng = np.random.default_rng(a.seed)
    results = []
    for p in [int(x) for x in a.stages.split(",")]:
        r = run(p, K, delta, a.samples, a.theta, rng, a.M, a.k0)
        results.append(r)
        print(json.dumps(r), flush=True)
    if a.output:
        json.dump(dict(note="Fable diagnostic (20 Oct 2026): rule-based global class assignment for the "
                            "tower adversary (honest lower bound on the adversary) vs the certificate's "
                            "comparison value", witness=a.witness, results=results), open(a.output, "w"), indent=1)


if __name__ == "__main__":
    main()
