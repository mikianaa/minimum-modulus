# Fable autonomous loop notes, 18 October 2026

Fable-owned working notes (no Codex files were edited). Everything below is
reproducible from `code/fable_*.py`; the Python environment is `tmp/fable-venv`
(python-sat, numpy, ortools).

## Upper side (certified bound 14,484, target 10,000)

1. **Per-stage loss profile** of the 14,484 certificate for the first 170
   stages (primes 2..1013), from an instrumented copy of the verifier:
   `data/fable-loss-profile-14485.json`. Losses are flat at 2.2-2.9 % of the
   surviving mass for p = 17..113 (the optimizer equalised marginal losses);
   0.857 of the total loss 0.9999 falls below p = 1013. For p <= 29 the
   distortion is tiny (delta_p < 0.01), so the loss equals the first moment
   E[S_p] and cannot be reduced by this kind of argument.
2. **Cap-inflation ceiling.** Replaying with every cap nu_q = 1/(1-delta_q)
   replaced by 1 lowers the early loss 0.857 -> 0.730 (ceiling 0.126, more than
   the ~0.09 needed): `data/fable-nu-cap-diagnostic-14485.json`. The mechanism
   is that in BBMST's distorted measure the weight nu_q/q^gamma is attained only
   on uncovered residues of fibres with alpha_q >= delta_q (heavy fibres), whose
   mass is at most M_1/delta_q (0.17 at p = 101, 0.10 at p = 199).
3. **Refutation.** The realisable gain of fibre-dependent caps
   w_q = min(nu_q, 1 + (nu_q-1) S_q(D_<q)/delta_q), measured by importance-
   sampled Monte Carlo on the nested comparison model
   (`code/fable_modulated_cap_mc.py`; the independent-model estimates agree
   with the certificate to < 1 %), is only 1.0 / 1.7 / 2.0 / 2.2 % of the
   stage loss at p = 101 / 199 / 307 / 503, against ceilings of 9 / 21 / 28 /
   38 %. Deep towers (the configurations that carry the excess) are heavy at
   every earlier stage, so they receive the nu inflation legitimately. The
   first moment drops ~7 % but the excess barely moves. Do not retry
   fibre-dependent caps. The comparison's extremal adversary never covers the
   tower; only coverage or a genuinely global constraint can help.

## Lower side (record 42; constructions)

Literature frontier found (minimal LCM of a distinct covering system with
minimum modulus m): m = 3: 120 (Krukenberg); 4: 360 (Dalton-Trifonov 2022,
JIS 22.9.1); 5: 1440 and 6: 5040 (Klein, arXiv:2508.18062); 7: 10080
(J. Zhang and S. Zhang, arXiv:2607.19029, July 2026, Gurobi filtering).
Open: m = 8 (Harrington-Klein-Lowrance-Trifonov, arXiv:2605.18644, give an
m = 8 system with LCM 2^8 3^3 5^2 = 172800 using only the primes 2, 3, 5 and
prove m <= 9 for {2,3,5}-smooth LCMs).

Tools: `code/fable_sat_cover.py` (python-sat; CaDiCaL ignores interrupts, so
timed runs use Glucose4) and `code/fable_cpsat_cover.py` (OR-Tools CP-SAT,
far faster), both exact: one Boolean per (modulus d | N, residue), at most one
residue per modulus, every class mod N covered, translation symmetry broken by
fixing residue 0 for modulus N. `code/fable_lcm_scan.py` adds the reciprocal
filter and a coarse exact relaxation at modulus g = N/p (infeasible => proof).

Validation against the literature (CP-SAT, 8 workers):

| m | N | result | time |
|---|---|--------|------|
| 3 | 108 / 120 | UNSAT / SAT | 0.0 s |
| 4 | 240 / 336 / 360 | UNSAT / UNSAT / SAT | 473 s / 29 s / 0.1 s |
| 5 | 1080 / 1440 | UNKNOWN (600 s) / SAT, 32 moduli | 51 s |
| 6 | 5040 | SAT, 55 moduli | 3.3 s |
| 7 | 5040 / 10080 | UNKNOWN / UNKNOWN (600 s) | |

SAT side reproduces Klein's values exactly. Beyond N ~ 10^3 the plain model
neither proves UNSAT nor finds the known m = 7 system in 10 minutes;
redundant counting constraints and LP linearisation did not help N = 240.
A first m = 8 scan (7-smooth N from 10080) returned only UNKNOWN for 10080 and
10800 and was stopped. Reaching m = 8 needs the structured filtering of
Zhang-Zhang (or a hierarchical/lifting search), not a bigger time limit.

**Capability test (19 Oct, 00:36).** CP-SAT with redundant counting
constraints, LP linearisation level 2 and 8 workers did not find the known
m = 7, N = 10080 covering system in 5131 s (UNKNOWN). The plain exact model is
therefore not a viable route to the open m = 8 case; the SAT line is closed
for this loop unless a structured (filtering / lifting) formulation is built.

**Coarse-filter power (19 Oct).** For N = 1080, m = 5 (UNSAT by Klein), the
exact relaxation at modulus g is FEASIBLE in 0.1 s for g = 36 and 54 (too
coarse to see the obstruction) and UNKNOWN for g = 72 even with hours of
wall time (the run spanned a system sleep). The coarse relaxations are not
cheap certificates at this size either; the obstruction needs the full
residue structure or a problem-specific argument.

## HKLT Problem 1 via the distortion method (19 Oct, 08:05)

Harrington-Klein-Lowrance-Trifonov (arXiv:2605.18644) ask whether a distinct
covering system with {2,3,5}-smooth moduli and minimum modulus 9 exists; they
prove m <= 9 by truncated inclusion-exclusion and expect m = 9 impossible.
The repo's comparison model needs only three stages here (p = 2, 3, 5), so
the bound can be evaluated exactly (`code/fable_smooth235_bound.py`): the
minimal total loss over the distortion parameters is 1.060 for K = 9 (not
excluded) and 0.960 for K = 10 (excluded). The method therefore reproduces
exactly HKLT's threshold m <= 9 and does not settle m = 9; closing the 0.06
gap would need an argument that uses coverage (the comparison's extremal
adversary is a non-covering tower), not better parameters.

## New task (19 Oct 08:10): construct a covering system with minimum modulus 43

The user redirected the loop to M = 43 and asked to continue until an
example is found. Starting analysis (from Owens' thesis text and Codex's
construction-search.tex):

* Owens' 42 construction contains exactly one class of modulus 42: the "2" in
  the third input of the top-level 7-arrow on the even branch,
  7↑(..., 3(2, 4, 3↑(1,2)), ...), i.e. H = 10 (mod 42). Every other modulus
  in Sections 3.1-3.4 is at least 45, so raising the minimum to 43 removes
  only this class. All later sections use moduli >= 43 automatically.
* The removed class covered, in the first input of that 3 inside the third
  input of the 7, all of the even branch: the hole 4 part (2 mod 4, density
  1/84), hole 8 (1/168), hole 16 (1/336), hole 32 (1/672); 0 mod 32 is
  covered by 64↑. Total new obligation about 1/42.
* Downstream, Owens' later budgets silently rely on this class wherever a
  7↑ is refilled on a "first input of 3" branch and its third input is
  marked x: Section 3.8 (three 7↑ copies on hole 4), 3.10 (four copies),
  3.12 (hole 8, "only three inputs" becomes four), 3.15 ("only four inputs"
  becomes five). Those counts were exact, so each costs one extra set.
* Codex's density audit (construction-search.tex, Prop. "restricted density
  obstruction") shows {2,3,5,7}-smooth spare moduli alone cannot cover H
  (capacity 2759/3600), and the fresh-prime set-building rule stalls
  (rational relaxation B < q-1 for q >= 127). A repair must therefore be
  coupled: rearrange early placements (prime 3 sets, the 81↑, the 36·7^d
  allocation) to create surplus sets.

Tooling to reuse: code/construction_kernel.py (exact finite-kernel
verifier) and code/owens_seed_core.py (literal encoding of Owens'
Sections 3.1-3.4, parametrised by the minimum modulus).

## RESULT (19 Oct, ~11:30): a candidate covering system with minimum modulus 43

`data/fable-m43-construction.json` records the construction:

1. Owens' 42-system with its single modulus-42 class 10 (mod 42) deleted
   (Sections 3.1-3.20 otherwise verbatim; later budgets only need H = 10 mod 42
   to be covered, which steps 2-3 provide).
2. The uncovered remainder H' of H after Sections 3.1-3.4 is the product
   (all 7-adic states of 3 mod 7) x C', with C' explicit (hole 4 x 5-classes
   1-4 and 20 mod 25; hole 8 x classes 1-4 and the 5-tail (20,4); hole 16 x
   classes 2-4; hole 32 x class 4; all in 1 mod 3); density 71/4200.
3. 14 spare classes/arrow families with {2,3,5,7}-smooth moduli unused by
   Owens (49, 98, 147, 343^, 686^, 1029^, 1372^, 1715^, 196, 245, 3360, 6720,
   23520, 47040) cover 7-adic rows and the hole-16/32 class-4 parts.
4. 99 progressions, each covered by an arrow q^ with its own fresh prime
   q in {89, 97, ..., 797} (89 is fresh because all of Owens' 89-moduli contain
   67, which is excluded as a base prime for q = 89). Each arrow's q-1 inputs
   come from Codex's disjoint set-building rule; a base-prime order reaching
   >= q-1 sets is recorded for every piece (99/99 certified). All resulting
   moduli contain q, hence are distinct from Owens', from the spares and from
   each other. Design found by CP-SAT over refinement shapes with one fresh
   prime per piece (scratchpad m43/design_ilp10.py); piece shapes are genuine
   progressions (no tail exponents).
5. Finitization with auxiliary primes > 797 as in Nielsen/Codex.

Verification: all 50,803,200 kernel states of H covered (depths 6,4,3,3);
40,000 random integers of H covered when rows are evaluated as congruences;
kernel/integer semantics agree on 3,000 samples with Codex's matches().

Caveats: inherits the unverified correctness of Owens' thesis (Sections
3.5-3.20 were never machine-checked) and of Codex's transcription of
Sections 3.1-3.4 (owens_seed_core.py). The exclusion of 67 for q = 89 rests on
the thesis text (89 is used only to fill one entry of a 67-arrow).

### Simplification and verification (19 Oct, afternoon)

* Re-solving the design to minimise the number of pieces gives **44 fresh-prime
  arrows with primes 89, 97, ..., 331** (previously 99 up to 797); coverage of
  H re-verified on the kernel (0 states remaining), 44/44 set-building orders
  certified. Shapes used: (4,1,2,2)x12, (3,2,1,2)x9, (3,1,1,1)x7, (3,1,1,2)x6,
  (2,2,1,2)x6, (4,1,1,1)x2, (3,1,2,1)x1, (5,1,1,2)x1. This is optimal for the
  present scheme (one fresh prime per progression, Codex's set-building rule).
* Independent transcription audit: Fable re-read Owens' Sections 3.1-3.4 and
  matched every family in owens_seed_core.py. One cyclic relabelling of the
  mod-25 classes inside 2 mod 5 (seventh-prime input 6) is a branch automorphism
  and does not affect the structure. The residual after the prime-7 stage
  matches the thesis diagrams hole by hole.
* Running: machine check that the symbolic modulus families of each fresh
  arrow are pairwise disjoint (scratchpad m43/check_families.py).

### The 44-piece patch at a glance

Spare classes (moduli unused by Owens, all >= 43):

| modulus | conditions |
|---|---|
| 49 | x = 3 (mod 49) |
| 98 | x = 0 (mod 2); x = 10 (mod 49) |
| 147 | x = 1 (mod 3); x = 17 (mod 49) |
| 343 | arrow family: x = 24 (mod 49) with first nonzero base-7 digit 1; one class at each level >= 3 |
| 686 | x = 0 (mod 2); arrow family: x = 24 (mod 49) with first nonzero base-7 digit 2; one class at each level >= 3 |
| 1029 | x = 1 (mod 3); arrow family: x = 24 (mod 49) with first nonzero base-7 digit 3; one class at each level >= 3 |
| 1372 | x = 2 (mod 4); arrow family: x = 24 (mod 49) with first nonzero base-7 digit 4; one class at each level >= 3 |
| 1715 | x = 4 (mod 5); arrow family: x = 24 (mod 49) with first nonzero base-7 digit 5; one class at each level >= 3 |
| 196 | x = 2 (mod 4); x = 31 (mod 49) |
| 245 | x = 4 (mod 5); x = 38 (mod 49) |
| 3360 | x = 8 (mod 32); x = 1 (mod 3); x = 4 (mod 5); x = 3 (mod 7) |
| 6720 | arrow family: x = 24 (mod 32) with first nonzero base-2 digit 1; one class at each level >= 6; x = 1 (mod 3); x = 4 (mod 5); x = 3 (mod 7) |
| 23520 | x = 16 (mod 32); x = 1 (mod 3); x = 4 (mod 5); x = 24 (mod 49) |
| 47040 | arrow family: x = 16 (mod 32) with first nonzero base-2 digit 1; one class at each level >= 6; x = 1 (mod 3); x = 4 (mod 5); x = 31 (mod 49) |

Fresh-prime arrows, grouped by the region they cover (hole = 2-adic class, class = residue of the prefix mod 5):

| region | shape 2^a 3^b 5^c 7^d | progression modulus | pieces | fresh primes |
|---|---|---|---|---|
| 2 mod 4 x 5-class 1 | (3,1,1,1) | 840 | 2 | 89, 103 |
| 4 mod 8 x 5-class 1 | (3,1,1,1) | 840 | 1 | 97 |
| 4 mod 8 x 5-class 3 | (3,1,1,1) | 840 | 1 | 101 |
| 2 mod 4 x 5-class 3 | (3,1,1,1) | 840 | 2 | 107, 113 |
| 4 mod 8 x 5-class 2 | (3,1,1,1) | 840 | 1 | 109 |
| 8 mod 16 x 5-class 3 | (4,1,1,1) | 1680 | 1 | 127 |
| 8 mod 16 x 5-class 2 | (4,1,1,1) | 1680 | 1 | 131 |
| 2 mod 4 x 5-class 2 | (3,1,1,2) | 5880 | 6 | 137, 139, 149, 151, 157, 163 |
| 4 mod 8 x 20 mod 25 | (3,1,2,1) | 4200 | 1 | 167 |
| 2 mod 4 x 5-class 4 | (2,2,1,2) | 8820 | 6 | 173, 179, 181, 191, 193, 197 |
| 16 mod 32 x 5-class 4 | (5,1,1,2) | 23520 | 1 | 199 |
| 4 mod 8 x 5-class 4 | (3,2,1,2) | 17640 | 9 | 211, 223, 227, 229, 233, 239, 241, 251, 257 |
| 2 mod 4 x 20 mod 25 | (4,1,2,2) | 58800 | 12 | 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331 |
* Modulus distinctness inside each fresh arrow: the explicit symbolic expansion
  of the set-building rule is exponential in the number of base primes (about
  2^60 classes for q ~ 300), so a full enumeration is impossible; the lemma's
  exponent-range argument is the proof. As a machine test, the rule was expanded
  with base primes <= 13, 23 and 31 for all eight shapes used (up to 8,292
  classes per instance): 0 duplicate moduli in every case.

**Verification status summary (19 Oct):** (i) hole H' computed on the kernel
and cross-checked by integer sampling; (ii) Owens 3.1-3.4 transcription audited
independently (structure identical); (iii) coverage of H by seed + 14 spares +
44 fresh arrows verified on all 50.8M kernel states; (iv) set counts certified
for all 44 arrows; (v) distinctness: structural proof + truncated machine tests.
Not verified: Owens' Sections 3.5-3.20 themselves (no transcription exists).

## RESULT (20 Oct): a candidate covering system with minimum modulus 44

`data/fable-m44-construction.json`. Owens' construction uses modulus 43 exactly
once: the set '1' input of the top-level 43-arrow in Section 3.14 on the branch
B = 2 mod 4 x 4 mod 5 (fourth 5-input) x 2 mod 3 (middle 3-input). Deleting it
opens K = B x 1 (mod 43) (ordering the arrow's inputs so that '1' sits in class
1); after Sections 3.1-3.4, 70% of B is uncovered (7-inputs 1,2,3,6 fully,
input 5 only in 24 mod 25, input 4 nothing, input 7 recursively; confirmed on
the kernel). Because every K-piece carries x = 1 (mod 43), the set-building
count doubles, and a single arrow at q = 89 would cover all of K (91 >= 88
sets) - but 89 is also needed by the M=43 patch, which is exactly one prime
short without it. Joint re-design (both holes, one fresh prime per piece,
set-building orders optimised further: +46 feasible (shape, prime) pairs) gives
**43 arrows in total, primes 89 and 97..317**: 31 for H' and 12 for K
(B x {3,6} mod 7 at 181, 191; B x 0 mod 7 at 193; B x 1 mod 7 split mod 9 at
271, 277, 281; B x 2 mod 7 split mod 16 at 283..311; B x 5 mod 7 x 24 mod 25
split mod 8 at 313, 317). The 14 spares are unchanged. Verified: kernel coverage
of both holes (0 of 50.8M states remaining), 43/43 set-building orders, integer
sampling of both holes (30,000 each, 0 uncovered; note: use t < 2^48 in the
sampler, int64 overflow produced a false failure at first).

Minimum modulus 45 would additionally require deleting the class of modulus
44 = 4 x 11 (Nielsen's 11-arrow, input '4' on 1 mod 4), a hole of density about
1/99 on the odd branch; see the capacity estimate below.

**Capacity estimate for minimum modulus 45 (20 Oct).** After the M=43/44
patches every fresh prime up to 317 is in use. Covering the modulus-44 hole
K44 = 1 mod 4 x 1 mod 11 x [(1 mod 3) u (3 mod 9)] (density about 0.0100)
with one fresh-prime arrow per progression, using only primes 331..6000 and
the best feasible shape for each prime, reaches a density of only 0.00085
(243 usable primes). The patch scheme therefore stops at minimum modulus 44;
minimum modulus 45 would need a new source of sets (borrowing from Owens'
structure, a reorganised early stage, or a stronger set-building rule), not
more primes.

## Minimum modulus 45: obstruction analysis (20 Oct)

Codex's transcription of Nielsen's prime-11 stage (code/nielsen_prime11.py,
data/nielsen-prime11-transcription.json) is VERIFIED_RELATIVE_ROOT: the ten
inputs are a1 = 4, a2 = 8^, a3 = 3·2+27·1+27^·2, a4 = 3·4+27·4+27^·8^,
a5 = 3·8^+9·8^, a6 = 5^(1,2,3·1,4), a7, a8, a9 = 7^(...), a10, on the region
W x (1 mod 4) with W = (1 mod 3) u (3 mod 9). The only modulus below 45 is
a1 at level 1: 11·4 = 44. Removing it opens K44 = W x 1 mod 4 x 1 mod 11,
density about 0.0100.

Supply for K44: fresh primes > 317 give 0.00085 (one arrow per progression,
best shapes); the 11-power spares 121, 242, 363, 1331, 2662, 3993, ... (Nielsen
never uses the sets 1, 2, 3 alone in the 11-arrow) cover about 27 % of K44
because K44 lies in 1 mod 11; the eight unused {2,3,5,7}-smooth spares cover
about 1 %. Total about 36 %: the patch method cannot reach 45.

Why no redesign of the 11/13 stage helps: on the odd branch the cheap 2-adic
sets 2 and 4 give moduli 22 and 44, so an 11-arrow can only live on a branch
of 2-adic depth >= 3 (sets 8, 16^), i.e. on a quarter (a mod 8) of W x odd,
and a second 11-arrow on another quarter would need 2-adically shifted sets
(16, 32^: no longer covering the quarter) or shifted-down sets (1, 2: below
45). The 13-arrow may use 4 (52) and covers one half (3 mod 4, or 1 mod 4 with
two doubling copies), but a second 13-arrow collides in moduli. Hence W x odd
= three quarters needs three top-level primes; every other prime of Owens is
top-level elsewhere with overlapping set moduli, and a fresh prime would need
>= 96 sets on a branch that offers about 12. Freeing another Owens prime (17
for hole 12, 23 for hole 24, 31 for hole 36) only moves the problem.

Conclusion: minimum modulus 45 requires a new design of the odd branch from
the prime-3 stage onward (Owens' own suggestion: change the sets used with the
prime 3 or the placement of the extra 81^), not a patch. The 42 -> 43 -> 44
patches work because the deleted classes lie on the even branch, where the
7-arrow structure left many spare moduli (49, 98, 147, ...).
