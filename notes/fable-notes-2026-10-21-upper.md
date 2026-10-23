# Fable notes, 21 October 2026: is the convex comparison tight? (upper side)

Fable-owned notes (no Codex files edited). Scripts: `code/fable_adversary_toy.py`
(exact small-instance adversaries); results `data/fable-adversary-toy-results.json`.
Environment `tmp/fable-venv` (scipy added for the LP).

## Question

The 20 Oct notes identified the truncation alpha <= 1 as the main slack of the stage
bound loss_p = E(S_p(D) - delta)_+/(1 - delta): against the TOWER family in the
independent model the truncated value is 21-25 % below the comparison value at
p = 43..503, and the sampled "branching" adversaries recovered none of it robustly.
That measurement fixed the measure to the independent model and restricted the family.
Here both are optimised exactly in small instances.

Instance: coordinates = small primes Q with exponents G (space X = prod Z_{q^G}); stage
prime p with fibre Z_{p^gamma}; moduli m p^j (m | prod q^G, m p^j >= K) each get ONE
residue mod m and ONE class mod p^j; alpha(x) = measure of the union of the classes of
the moduli active at x. Loss = E_P(alpha - delta)_+/(1 - delta).

* cap-only adversary: P ranges over ALL measures with BBMST's conditional atom caps
  P(X_k = t | X_<k) <= nu_k/q_k^G_k (an LP for a fixed family, scipy/HiGHS), the family
  is annealed on top (residues, classes, class swaps);
* uniform adversary: nu = 1 at every coordinate forces P uniform; only the family is
  annealed, the loss is evaluated exactly.

Reference values: the comparison (independent model, V_q = min(E_q, G_q)) with and
without truncation min(S, 1), and the tower family's actual value.

## Results

Cap-only adversary (LP measure + annealed family, 8000 steps; all three coordinates
have nu = 1.25 unless stated, so the measure may zero out 20 % of every fibre):

| instance | comparison | truncated | LP+anneal | ratio to comparison |
|---|---|---|---|---|
| Q=2,3,5 G=3,2,1 p=7 K=28 delta=.25 (seed 11 / 12) | 0.2878 | 0.2134 | 0.2545 / 0.2716 | 0.88 / 0.94 |
| same, K=14 | 0.4365 | 0.3170 | 0.4065 | 0.93 |
| same, nu=1.5, delta=.3 | 0.4147 | 0.2823 | 0.3728 | 0.90 |

The optimal measures are two-level: a fraction 1/nu of every fibre at the cap, the
rest at zero (this is BBMST's measure of a stage-q cover with alpha_q = delta_q
exactly everywhere). The anneals were still improving at the end.

Uniform measure (nu = 1 at every coordinate, so P is forced uniform; this is the real
regime for coordinates q <= 29, where the certificate's nu_q - 1 is 1e-4 .. 8e-3):

| instance | comparison | truncated | tower (actual) | best family | ratio to comparison | above truncated |
|---|---|---|---|---|---|---|
| Q=2,3,5 G=3,2,1 p=7 g=2 K=28 d=.25 (2 seeds) | 0.1828 | 0.1424 | 0.1228 | 0.1749 / 0.1742 | 0.957 / 0.953 | +23 % |
| same, K=14 | 0.2909 | 0.2231 | 0.1840 | 0.2773 | 0.953 | +24 % |
| same, gamma=3 | 0.1901 | 0.1470 | 0.1197 | 0.1828 | 0.962 | +24 % |
| G=4,2,1, K=56 | 0.1205 | 0.0896 | 0.0706 | 0.1129 | 0.937 | +26 % |
| Q=2,3,5,7 G=3,2,1,1 p=11 g=2 K=66 d=.25 | 0.0898 | 0.0730 | 0.0343 | 0.0836 | 0.931 | +14 % |
| same, K=110 d=.2 | 0.0613 | 0.0518 | 0.0307 | 0.0577 | 0.941 | +11 % |
| same, K=66 d=.4 | 0.0739 | 0.0530 | 0.0173 | 0.0605 | 0.819 | +14 % |
| Q=2,3,5,7 G=4,2,1,1 p=31 g=1 K=248 d=.2 | 0.00985 | 0.00956 | 0.00815 | 0.00878 (unconverged) | 0.89 | -8 % |
| same, K=496 d=.15 | 0.00662 | 0.00647 | 0.00604 | 0.00604 (no move accepted) | 0.91 | -7 % |

In every instance with real saturation (p = 7, 11) the best family has NO saturated
fibre (frac_saturated 0 .. 0.6 %): it spreads the residues so that alpha stays in
[delta, 1] on most of the active mass, and reaches 93-96 % of the untruncated
comparison value, 11-26 % ABOVE the truncated tower value. The p = 31 instances have
only 3 % truncation to begin with; the annealer in a 5040-point space barely moves,
so 0.89-0.91 are lower bounds for the adversary there (the gap is class collisions on
31 classes plus sub-threshold waste, not truncation).

## Conclusions

1. Truncation (alpha <= 1) is NOT a usable lead. The 20 Oct estimate that 70-90 % of
   the overflow is unrealisable was an artefact of fixing the measure and restricting
   the family moves; the exact adversary water-fills under the uniform measure and
   recovers nearly all of it. The one-coordinate toy generalises.
2. With measure freedom (nu > 1 at the saturating coordinates) the cap-only adversary
   gets to 88-94 % and was still improving: no cap-only lemma can give the truncated
   value either.
3. The convex comparison is therefore within roughly 5-10 % of the exact single-stage
   adversary in these instances, and the residual is a combinatorial optimum (packing
   of residue classes), not something a clean inequality will capture. The target
   10,000 needs 11 % over p >= 41: it is out of reach of any single-stage refinement
   of the comparison (truncation, collisions, modulated caps, measure structure),
   not just of parameter tuning.
4. What remains: (a) at coordinates q <= 29 the measure is forced uniform, so the
   stages decouple there and the per-stage bound is essentially the truth; (b) the
   only coupling left is at coordinates q >= 31, where BBMST's two-level measure
   (cap on uncovered points of heavy fibres only) is produced by the stage-q cover of
   the SAME family; charging that shaping would need a multi-stage accounting that the
   lineage experiments (17 Oct) showed to be delicate; (c) a different distortion
   rule. None of these is a bounded experiment; the sensible next upper-side step is a
   multi-stage toy (two consecutive stages, one family, exact adversary) to measure
   how much the coupling in (b) is worth before investing in a lemma.
5. For the headline number: the block tail (section "A rigorously accelerated block
   tail") gives 3.7e-5 of extra loss budget at 14,485, and the early-stage sensitivity
   is about 1.9e-5 of loss per unit of cutoff, so re-running the prefix with the block
   tail would move 14,484 to about 14,482. Not worth a 70-minute certificate.

