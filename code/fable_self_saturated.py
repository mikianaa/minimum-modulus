"""How much of the tower's overflow E(S-1)_+ is carried by self-saturated moduli?

A modulus m with S_p(m) > 1 (its own divisors already saturate a fibre) is active
only on saturated fibres of the tower family, so its whole weight is wasted there;
a non-tower family can move such moduli to shallower branches with headroom at no
cost in mass.  If the weight of these moduli on saturated fibres is comparable to
the overflow E(S-1)_+, the truncation gap of the tower is recoverable in principle.

Reports, per stage: comparison value, overflow E(S-1)_+/(1-delta), the weight of
self-saturated moduli (S(m) > 1) and of nearly self-saturated ones (S(m) > 1/2) on
saturated fibres, all importance-sampled.
Usage: python code/fable_self_saturated.py --stages 101,199,307 --samples 20000
Fable-owned file; see memory fable-loop-plan.
"""
import argparse, json
import numpy as np
from fable_class_collision_mc import load_params, j0, sample_stage


def run(p, K, delta, n, theta, rng):
    dp = delta[p]
    W, divl, nu = sample_stage(p, K, delta, n, theta, rng)
    cp = p / (p - 1)
    Fc = {}

    def F(m):
        v = Fc.get(m)
        if v is None:
            v = cp * p ** (-j0(m, p, K)); Fc[m] = v
        return v

    Sm = {}

    def S_of(m):
        v = Sm.get(m)
        if v is None:
            # divisors of m among the sample's divisor list are not available; compute directly
            divs = [1]
            x = m
            for q in sorted(nu):
                if x % q == 0:
                    e = 0
                    while x % q == 0:
                        x //= q; e += 1
                    divs = [d * q ** k for d in divs for k in range(e + 1)]
            v = sum(F(d) for d in divs); Sm[m] = v
        return v

    acc = dict(comp=0.0, over=0.0, selfsat=0.0, halfsat=0.0, satmass=0.0, first=0.0)
    for w, divs in zip(W, divl):
        S = sum(F(m) for m in divs)
        acc["first"] += w * S
        acc["comp"] += w * max(S - dp, 0)
        if S > 1:
            acc["over"] += w * (S - 1)
            acc["satmass"] += w
            acc["selfsat"] += w * sum(F(m) for m in divs if S_of(m) > 1)
            acc["halfsat"] += w * sum(F(m) for m in divs if S_of(m) > 0.5)
    out = dict(p=p, delta=dp, samples=n)
    for k, v in acc.items():
        out[k] = v / n / (1 - dp) if k != "satmass" else v / n
    out["selfsat_over_overflow"] = out["selfsat"] / out["over"] if out["over"] else None
    out["halfsat_over_overflow"] = out["halfsat"] / out["over"] if out["over"] else None
    out["overflow_share_of_comparison"] = out["over"] / out["comp"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--witness", default="data/sparse-certificate-14485.json")
    ap.add_argument("--stages", default="101,199,307")
    ap.add_argument("--samples", type=int, default=20000)
    ap.add_argument("--theta", type=float, default=2.5)
    ap.add_argument("--seed", type=int, default=9)
    a = ap.parse_args()
    K, delta = load_params(a.witness)
    rng = np.random.default_rng(a.seed)
    for p in [int(x) for x in a.stages.split(",")]:
        print(json.dumps(run(p, K, delta, a.samples, a.theta, rng)), flush=True)


if __name__ == "__main__":
    main()
