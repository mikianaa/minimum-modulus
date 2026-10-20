"""Three-stage distortion (BBMST comparison-model) bound on the covered density of a
finite family of arithmetic progressions with DISTINCT {2,3,5}-smooth moduli >= K.

Stages: p=2 (moduli 2^j>=K), p=3 (2^a 3^b, b>=1, >=K), p=5 (2^a 3^b 5^c, c>=1, >=K).
Stage loss = E(S_p(D)-delta_p)_+/(1-delta_p) with the nested comparison model
P(E_q>=a)=nu_q/q^a, nu_q=1/(1-delta_q), S_p(D)=sum_{m|D} F_p(m),
F_p(m)=sum_{j>=1, m p^j>=K} p^-j  (exactly as in src/conditional-stoploss-improvement.tex).
If the minimal total loss is < 1 no such family covers the integers.

Result (18/19 Oct 2026): K=9 -> 1.060 (not excluded), K=10 -> 0.960 (excluded).
This reproduces Harrington-Klein-Lowrance-Trifonov's m<=9 (arXiv:2605.18644, Thm 1.12)
but does not settle their Problem 1 (m=9). Floating point; exponent tails truncated at
2^60, 3^40 (negligible). Fable-owned file; see memory fable-loop-plan.
Usage: python code/fable_smooth235_bound.py [K ...]
"""
import sys
import numpy as np


def F(p, m, K):
    j = 1; v = m * p
    while v < K: v *= p; j += 1
    return p ** (-j) * p / (p - 1)


def tables(K, A=60, B=40):
    S2 = np.array([sum(F(3, 2 ** i, K) for i in range(a + 1)) for a in range(A)])
    S3 = np.array([[sum(F(5, 2 ** i * 3 ** j, K) for i in range(a + 1) for j in range(b + 1)) for b in range(B)] for a in range(A)])
    return S2, S3


def bound(K, S2, S3, d1, d2, d3):
    A, B = S3.shape
    nu2 = 1 / (1 - d1); nu3 = 1 / (1 - d2)
    P2 = np.array([1 - nu2 / 2] + [nu2 / 2 ** (a + 1) for a in range(1, A)])
    P3 = np.array([1 - nu3 / 3] + [nu3 * 2 / 3 ** (b + 1) for b in range(1, B)])
    loss1 = max(F(2, 1, K) - d1, 0) / (1 - d1)
    loss2 = (P2 * np.maximum(S2 - d2, 0)).sum() / (1 - d2)
    loss3 = (P2[:, None] * P3[None, :] * np.maximum(S3 - d3, 0)).sum() / (1 - d3)
    return loss1 + loss2 + loss3, (loss1, loss2, loss3)


if __name__ == "__main__":
    Ks = [int(x) for x in sys.argv[1:]] or [9, 10]
    grid = [i / 100 for i in range(0, 51)]
    for K in Ks:
        S2, S3 = tables(K)
        best = min((bound(K, S2, S3, a, b, c) + ((a, b, c),) for a in grid for b in grid for c in grid), key=lambda t: t[0])
        print(f"K={K}: minimal total loss {best[0]:.5f} (stage losses {tuple(round(float(x),4) for x in best[1])}, deltas {best[2]}); "
              f"undistorted sum {bound(K, S2, S3, 0, 0, 0)[0]:.4f}; {'EXCLUDED' if best[0] < 1 else 'not excluded'}")
