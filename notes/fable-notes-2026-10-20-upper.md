# Fable notes, 20 October 2026: upper side, target 10,000

Fable-owned working notes (no Codex files edited). Scripts: `code/fable_class_collision_mc.py`,
`code/fable_class_collision_opt.py`, `code/fable_class_collision_rule.py`,
`code/fable_branching_adversary.py`; results `data/fable-class-collision-*.json`,
`data/fable-branching-adversary-*.json`. Environment `tmp/fable-venv`.

## Starting point

The certified bound is 14,484 (`data/latest-certified-bound.json`; the README headline
14,496 is stale). The stage loss used by the certificate is the comparison value

    loss_p = E(S_p(D) - delta_p)_+ / (1 - delta_p),   S_p(D) = sum_{m | D} F_p(m),

with D the independent-model divisor pattern of the earlier primes. The 18 Oct
Monte Carlo (`fable_modulated_cap_mc.py`, exact S_p with all divisors) agrees with the
certificate to < 1 %, so the *evaluation* of the comparison model has no slack; any
improvement must come from the comparison itself. Codex's own routes (cap tuning,
lineage accounting, per-stage truncation inside the lemma) are refuted; see memory
`minimum-modulus-14484-audit`.

Loss budget (profile `data/fable-loss-profile-14485.json`), with the ratio loss /
first moment: 0.93 at p <= 29 (irreducible first moment, 0.12 in total), 0.72 at 43,
0.42 at 101, 0.17 at 307, 0.06 at 1009. The comparison only does real work for
p >= 40; stages p >= 41 carry about 0.83 of the total loss 0.9999. Since cutoff 10,000
needs about 0.09 more survivor mass, a uniform 11 % reduction of the stage losses for
p >= 41 would reach the target.

## The unexploited structure: one class per modulus

At stage p a modulus d = m p^j covers ONE class modulo p^j of the new coordinate, and
that class c_d is the same on every fibre x with x in H_d (x = a_d mod m). The real
covered fraction of the fibre over x is

    alpha(x) = mu( union_{d : x in H_d} c_d )  <=  min(1, S_p(x)),

while the comparison value lets every fibre use pairwise disjoint classes
(alpha = S) and ignores alpha <= 1. Two consequences:

1. **Truncation.** Saturated fibres (S > 1) carry 21-25 % of the stage excess at
   p = 101, 199, 307 (E(min(S,1)-delta)_+ vs E(S-delta)_+ in the independent model).
2. **Collisions.** Level-1 moduli (m >= K/p, weight 1/p each) share only p classes.
   A heavy fibre has N_act = tau_big(D) active moduli; the adversary needs a colouring
   m -> Z_p that is injective on the active set of every heavy fibre simultaneously.
   The conflict graph is complete (any two moduli are co-active on some heavy fibre),
   so collisions are forced; the question is their weight.

For the TOWER adversary (all residues compatible with one point x*) both effects
are rigorous in the independent model: alpha is an increasing function of the
agreement vector, and the caps give stochastic domination for increasing functions.
The comparison lemma needs more (arbitrary residues), and Codex's two-event toy shows
truncation fails there: a family can "water-fill" by spreading moduli with a common
factor over incompatible residues (e.g. all moduli with 2^4 || m over the classes
mod 16); such a family is globally consistent, so no comparison lemma of the present
kind can prove truncation. The real adversary is somewhere between.

## Measurements (stage loss relative to the comparison value)

Tower adversary, independent model, importance sampling (theta = 2.5):

| p | delta | truncation only | collisions: random colouring | sample-fitted optimum (train) | deterministic rule (M = 3e6, k0 = 30) |
|---|---|---|---|---|---|
| 101 | 0.182 | -22.7 % | -47 % | -5 % below truncated | -13.7 % below truncated |
| 199 | 0.269 | -24.4 % | | -4.6 % | -26.7 % (rule covers too few moduli: K/p = 72) |
| 307 | 0.288 | -25.3 % | | -3.4 % | -32.5 % (same caveat, K/p = 47) |

The sample-fitted optimum (coordinate ascent on 40-60k fibres) is optimistic for the
adversary (free colours for rarely seen moduli); the deterministic greedy rule is a
valid family but weak for large p. Honest estimate of the forced collision cost for the
tower adversary: 5-10 % of the truncated stage loss.

Branching (non-tower) adversary: each modulus may move from its tower cone to a
top-level sibling cone {x = s q^(b-1) mod q^b} (coordinate ascent from the all-tower
family, moves only with >= 40 supporting fibres). Train: recovers 27 % (p = 101) and
8 % (p = 199) of the truncation gap; test (fresh fibres, unseen moduli stay tower):
recovers nothing (value below the tower). The water-filling freedom exists in
principle but this restricted adversary cannot exploit it robustly.

Full sweep (test sample; "collisions" = sample-fitted optimum on the training sample,
optimistic for the adversary, measured relative to the truncated value):

| p | delta | truncation only | collisions (train optimum) |
|---|---|---|---|
| 43 | 0.074 | -20.7 % | -2.3 % |
| 61 | 0.119 | -22.2 % | -3.0 % |
| 79 | 0.156 | -21.1 % | -3.5 % |
| 101 | 0.182 | -22.7 % | -4.9 % |
| 199 | 0.269 | -24.4 % | -4.6 % |
| 307 | 0.288 | -25.3 % | -3.4 % |
| 503 | 0.349 | -24.8 % | -4.1 % |
| 1009 | 0.370 | -31.7 % | -3.7 % |

The truncation share of the stage excess is a uniform 21-25 % for 43 <= p <= 503 and
rises to 32 % at p = 1009 (data/fable-class-collision-p43-79.json, -p101, -p199-307,
-p503-1009). Forced class collisions add only a few per cent against a well-chosen
colouring; the large gaps of random colourings (40-70 %) are not forced.

## Is the truncation gap recoverable by water-filling? (`code/fable_self_saturated.py`)

In a single coordinate the answer is yes: the overflow of a saturated branch is carried
by the deepest moduli, which are active only on that branch, and they can be moved to
an earlier-diverging branch with headroom at no cost in mass (the two-event toy). In
the multi-coordinate model the overflow at a saturated fibre is instead carried by
moduli that are shared with many unsaturated fibres. Weight on saturated fibres
(S > 1) of self-saturated moduli (S_p(m) > 1, their whole cone is saturated, so the
tower wastes them and a family may move them freely) and of half-saturated moduli
(S_p(m) > 1/2), as a fraction of the overflow E(S-1)_+:

| p | overflow / comparison | self-saturated / overflow | half-saturated / overflow |
|---|---|---|---|
| 43 | 0.205 | 0.158 | 0.390 |
| 61 | 0.209 | 0.138 | 0.356 |
| 79 | 0.221 | 0.120 | 0.327 |
| 101 | 0.232 | 0.109 | 0.300 |
| 199 | 0.237 | 0.070 | 0.231 |
| 307 | 0.226 | 0.045 | 0.174 |
| 503 | 0.275 | 0.038 | 0.137 |
| 1009 | 0.255 | 0.021 | 0.088 |

A modulus with small S_p(m) has most of its cone on unsaturated fibres, where its
weight produces real excess; moving it to a sibling branch shifts its activity to
lighter fibres (the comparison's convexity works against the adversary there). So
at most roughly the half-saturated share of the overflow is recoverable by any
rearrangement of residues, and 70-90 % of the overflow (15-23 % of the stage loss,
growing with p) appears to be slack that NO family can realise. The sample-fitted
branching experiments are consistent with this (no robust recovery).

## Assessment

* Against the strongest adversary constructed here the comparison value is 25-35 %
  too large at p in [100, 300]; a rigorous capture of 11 % over p >= 41 would give
  10,000. So the target is plausibly TRUE within a refined comparison method, but
* there is no proof route yet. The coordinate induction behind the comparison lemma
  allows the arrangement of classes to depend on the earlier coordinates; the real
  constraint is a GLOBAL assignment (one class and one residue per modulus). The
  union within a colour produces step functions after one coordinate is integrated
  out, and only Markov-type bounds survive, which are weaker than the present bound.
* Needed: a lemma controlling globally consistent assignments. The concrete target
  suggested by the measurements: for any family F and capped Q,

      E_Q (min(S_F,1) - delta)_+  <=  E_model (min(S_T,1) - delta)_+  +  R,

  where R bounds the water-filling gain and is at most the saturated-fibre weight of
  moduli whose cones are mostly saturated (self/half-saturated share above). The
  one-coordinate toy shows R cannot be dropped; the multi-coordinate structure
  (one residue per modulus, moduli shared between fibres) is what should make R
  small. A proof must keep the global consistency of residues across fibres, which
  the coordinate-by-coordinate induction of the present comparison lemma discards.
  Class collisions (2-5 %) are a minor second-order effect by comparison.
