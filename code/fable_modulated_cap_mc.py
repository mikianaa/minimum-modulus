"""Monte Carlo diagnostic for the fibre-dependent (modulated) cap idea.

Nested comparison model of the 14,484 certificate at stage p (cutoff K):
independent E_q (q<p) with P(E_q>=a)=nu_q/q^a, D=prod q^{E_q},
S_p(D)=sum_{m|D} F_p(m), F_p(m)=sum_{j>=1, m p^j>=K} p^{-j},
loss_p = E(S_p(D)-delta_p)_+/(1-delta_p).

Modulated model (NOT a proved bound; the hypothesis under test): the tail
factor nu_q is replaced by w_q = min(nu_q, 1+(nu_q-1) S_q(D_{<q})/delta_q),
i.e. the nu_q inflation is only available in proportion to how heavy the
fibre already was at stage q.

Both models are estimated from one importance-sampled stream (proposal tails
theta_q/q^a) with likelihood-ratio weights, so the difference has low variance.
Usage: python code/fable_modulated_cap_mc.py --stages 101,199,307 --samples 200000
Requires numpy (tmp/fable-venv). Fable-owned file; see memory fable-loop-plan.
"""
import argparse, json, math, time, sys
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


def S_stage(divs, p, K):
    c = p / (p - 1)
    return sum(c * p ** (-j0(m, p, K)) for m in divs)


def run(p, K, delta, nsamples, theta_max, rng):
    qs = [q for q in sorted(delta) if q < p]
    nu = {q: 1 / (1 - delta[q]) for q in qs}
    dp = delta[p]
    theta = {q: min(theta_max, 0.9 * q) for q in qs}
    acc = dict(ind_loss=0.0, mod_loss=0.0, ind_loss2=0.0, mod_loss2=0.0, ind_m1=0.0, mod_m1=0.0, hits=0)
    for _ in range(nsamples):
        divs = [1]
        w_ind = 1.0; w_mod = 1.0
        for q in qs:
            th = theta[q]
            u = rng.random()
            # proposal: P(E>=a)=th/q^a
            if u >= th / q:
                e = 0
            else:
                e = 1
                while u < th / q ** (e + 1): e += 1
            nuq = nu[q]
            if delta[q] > 0 and nuq > 1.0000001:
                Sq = S_stage(divs, q, K)
                wq = min(nuq, 1 + (nuq - 1) * Sq / delta[q])
            else:
                wq = nuq
            if e == 0:
                w_ind *= (1 - nuq / q) / (1 - th / q)
                w_mod *= (1 - wq / q) / (1 - th / q)
            else:
                w_ind *= nuq / th
                w_mod *= wq / th
                divs = divs + [d * q ** k for d in divs for k in range(1, e + 1)]
        S = S_stage(divs, p, K)
        ex = max(S - dp, 0.0)
        acc["ind_m1"] += w_ind * S; acc["mod_m1"] += w_mod * S
        acc["ind_loss"] += w_ind * ex; acc["mod_loss"] += w_mod * ex
        acc["ind_loss2"] += (w_ind * ex) ** 2; acc["mod_loss2"] += (w_mod * ex) ** 2
        acc["hits"] += ex > 0
    n = nsamples
    out = dict(p=p, delta=dp, samples=n, hits=acc["hits"])
    for k in ("ind", "mod"):
        mean = acc[k + "_loss"] / n
        var = acc[k + "_loss2"] / n - mean ** 2
        out[k + "_loss"] = mean / (1 - dp)
        out[k + "_loss_se"] = math.sqrt(max(var, 0) / n) / (1 - dp)
        out[k + "_first_moment"] = acc[k + "_m1"] / n
    out["gain_relative"] = 1 - out["mod_loss"] / out["ind_loss"] if out["ind_loss"] else None
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--witness", default="data/sparse-certificate-14485.json")
    ap.add_argument("--stages", default="101,199,307")
    ap.add_argument("--samples", type=int, default=100000)
    ap.add_argument("--theta", type=float, default=2.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--output")
    a = ap.parse_args()
    K, delta = load_params(a.witness)
    rng = np.random.default_rng(a.seed)
    results = []
    for p in [int(x) for x in a.stages.split(",")]:
        t = time.time()
        r = run(p, K, delta, a.samples, a.theta, rng)
        r["seconds"] = round(time.time() - t, 1)
        results.append(r)
        print(json.dumps(r), flush=True)
    if a.output:
        json.dump(dict(note="Fable diagnostic; modulated-cap model is a hypothesis, not a proved bound", results=results), open(a.output, "w"), indent=1)


if __name__ == "__main__":
    main()
