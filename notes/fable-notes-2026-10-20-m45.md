# Minimum modulus 45: determination (Fable, 20 October 2026)

Question (user): the record is believed to be 43; push it to 44 and 45.

## Status of 43 and 44

Both already exist in this repository as kernel-verified candidates, conditional on
Owens' thesis and on Codex's transcription of its Sections 3.1-3.4:

* `data/fable-m43-construction.json` (19 Oct): Owens minus 10 (mod 42) + 14 spare
  {2,3,5,7}-smooth classes + 44 fresh-prime arrows (primes 89, 97..331).
* `data/fable-m44-construction.json` (20 Oct): additionally deletes the single modulus-43
  class (Section 3.14, branch B = 2 mod 4 x 4 mod 5 x 2 mod 3) and covers both holes with
  the same 14 spares + 43 fresh-prime arrows (primes 89, 97..317).

So the candidate record is **44**, not 43. Nothing new was needed for 44 in this session.

## 45: what has to change

Every modulus of the 44-construction is >= 45 except one: the class of modulus 44 = 4 x 11,
the first input ("4") of Nielsen's top-level 11-arrow on the region W x (1 mod 4), where
W = (1 mod 3) u (3 mod 9) (Owens Section 3.5, Nielsen Section 4.5). 42 and 43 occur once
each and are already handled; 44 occurs once. Raising the minimum from 44 to 45 therefore
means: delete that one class and cover what it covered.

What the class 44 covered: x = 1 (mod 4), x = 1 (mod 11), density 1/44. Three parts matter:

1. **K44 = W x 1 mod 4 x 1 mod 11** (density 1/99 = 0.0101 before overlaps): the part the
   11-arrow itself needs.
2. **Hole 12 x 1 mod 11** (hole 12 = 2 mod 3 x 1 mod 4, density 1/132) and **hole 36 x 1 mod
   11** (hole 36 = 6 mod 9 x 1 mod 4, density 1/396): Owens' 17-stage (3.7, via Nielsen's
   19-stage sets) and 31-stage (3.11) both rely on the classes 44 and 88^ ("the 11^ copies need
   only eight entries"). With 44 gone, every 11^ copy inside those stages has an uncovered
   first input: about 2/17 of hole 12 x 1 mod 11 plus the copies inside the composite set and
   inside the 13^ (about 7 cells of size 1/(12*17*13*11)), and 3/31 of hole 36 x 1 mod 11.
   Secondary holes total about 0.0012-0.0014 of density, all lying in 1 mod 11.

## Exact resource accounting (scratchpad m45/frac_account.py)

Resource = exponent tuple over (2,3,5,7) attached to a prime level; each resource is one
modulus and may be used once. Fraction = density of the resource's class inside the arrow's
region, relative to the region. An input of a p-arrow needs fraction >= 1 (minus what the
early stages already cover there).

* 11-arrow on W x 1 mod 4 (moduli 11 s): 10 inputs, 417 resources, fraction sum 9.646
  (inputs 3,4,6,7,8,9,10 are below 1 and rely on Owens' seed classes 45, 125, 135, 375, ...
  which Nielsen marks with x). **Unused 11 s resources with positive fraction: none.**
* 13-arrow on W x 3 mod 4 (moduli 13 s): 10 direct inputs, same structure, fraction 9.646.
  Its odd-modulus sets also cover the 11-side for free: fraction on 1 mod 4 is 0.958 (input 3),
  0.688 (6), 0.604 (8), 0.542 (9), 0.406 (10), 0.250 (7), 0 for inputs 1,2,4,5.
* 143-level (the two 11^ copies filling inputs 11 and 12 of the 13-arrow, Nielsen's
  "doubling"): 420 resources; the collapsed copy (sets 1, 2) covers both halves, so
  K44 n (12 mod 13) is covered by 143 x 1. **Unused 143 s resources with positive fraction:
  none.**

Consequences:

* No 10th set exists for the 11-arrow at the 11-level: the only lost resource is 11 x 4 and
  nothing at 11 s is free.
* A 13^ sub-arrow inside the lost input would need 143 s resources: all used by the doubling
  copies. In fraction units the 143 level supplies about 18.6 units (counting both halves)
  against a demand of about 22 (about 9 for a 13^ in the lost input after borrowing, about 13
  for the two copies), so no re-planning of the 143 level can serve both arrows.
* A 17^ sub-arrow (187 s) would need 16 full sets; the region offers about 13 even if the
  187-level were entirely free, and Owens' 17-stage already uses 187 s moduli on hole 12.
  19^ (18 sets), 23^ (22) are further out of reach.
* Freed level->=2 resources of the 11-arrow (121 x {1,2,3}, 1331 x {1,2,3}, ...) cover
  (1 + 1 + 3/4)/10 = 27.5 % of anything that lies in 1 mod 11 and odd (1 mod 3 for the 363
  family), i.e. 0.00278 of density at most.
* 169-level spares (169 x {1,2,3}, ...) are inefficient because K44 is spread over all classes
  mod 13: under 2 % of K44.
* Fresh primes: the set-building rule with base primes 2,3,5,7,11,13 and primes 331..6000
  (all primes <= 317 are consumed by the 43/44 patches) has capacity sum 1/m = 0.00085
  (scratchpad m45/fresh_capacity.py; adding 13 as a base prime does not raise it: the best
  shapes for q = 331..600 are 2^3 3^2 5 7 11 .. 2^4 3^3 5 7 11).

Kernel figures (scratchpad m45/odd_branch.py; primes 2,3,5,7,11,13, depths 6,4,3,2,2,1, odd
2-adic states only, 137 M states; seed = Owens 3.1-3.4 at minimum 45, Nielsen's 11-arrow with
the 44 class removed, Nielsen's 13-arrow with both doubling copies):

* Slice W x 1 mod 4 x (1 mod 11) at minimum 45 (9,676,800 states, 1,634 rows, modulus-44 row
  removed): **uncovered 4,584,552 states = density 0.005264 = 52.1 % of 1/99**. By class mod 13
  (1..12): classes 1, 2, 4, 5, 11 entirely uncovered (628,992 states each: the 13-arrow's inputs
  there use only 2-adic-specific sets), class 3 and class 12 entirely covered (odd sets, and the
  collapsed copy's 143 x 1), classes 6..10 partly (173,320 / 450,240 / 218,904 / 240,184 /
  356,944). 1 mod 3 part 3,834,768, 3 mod 9 part 749,784.
* Sanity slices: the same model at minimum 44 (modulus-44 row kept) leaves 0 states
  uncovered in the class-1 slice, and at minimum 45 the slices for 11-classes 2 and 6 are fully
  covered (0 states). So the transcription is complete relative to the seed and the hole is
  exactly the class-1 slice above.
* Level->=2 spares of the 11-arrow (121, 242, 363 and their 11-power continuations) cover
  28.4 % of this hole; remaining density 0.00377.
* Free {2,3,5,7}-smooth resources (unused by the seed and by the 14 spares of the 43/44
  patches) with positive coverage: 4 (756, 2268, 5292, 15876), together about 3 % of the hole.
* Fresh-prime capacity (331..6000): 0.00085. **Deficit about 0.0028 of density on K44 alone**,
  before the secondary holes (about 0.001) in the 17- and 31-stages.

The secondary figure (hole 12 x 1 mod 11 = 1/132, hole 36 x 1 mod 11 = 1/396 were covered by
the class 44 outright; Owens re-covers most of them in Sections 3.7 and 3.11 except inside the
11^ copies) is from the thesis text, not from the kernel, since those stages have no transcription.

## Verdict

Minimum modulus 45 is **not reachable by any patch** of the Owens construction: after the
level->=2 spares and all free smooth resources, the remaining part of K44 alone exceeds the
fresh-prime capacity several times over, and the secondary holes in the 17- and 31-stages
(about 0.001) exceed it by themselves. The obstruction is structural: on the odd branch the
only resource lost at 45 is 11 x 4, but Nielsen's 11/13 pair uses every 11 s, 13 s and 143 s
resource, the 17/23/31 stages are exact, and the even branch's primes (19, 29, 37, ...) are
exact as well. Variants examined and rejected by set counting: 11 on W x 1 mod 8 (keeps 10
sets, but W x 5 mod 8 then needs a third top-level prime with ~16 sets where ~11-13 exist);
13 on three classes mod 8 (loses sets); 11 on (1 mod 3) x 1 mod 4 (fits: ~10-11 sets, but
3 mod 9 x 1 mod 4 (density 1/36) has no cover: a 31-type arrow on (3 u 6 mod 9) x 1 mod 4
counts only ~22-25 sets, and 17 on (2 mod 3 u 3 mod 9) x 1 mod 4 counts ~13-14 of 16).

A construction with minimum modulus 45 needs a new early design (Owens' own suggestion: the
prime-3 sets or the placement of the extra 81^), i.e. a different hole geometry on the odd
branch, not more primes. The patch line of 42 -> 43 -> 44 ends at 44.

## 21 Oct: redesign attempt by set-count calculus (user: "keep going on 45")

Instead of patching, re-partition the odd branch. Facts used:

* Hole geometry after the prime-3 stage is invariant up to relabelling: the {2,3}-smooth
  resources >= 45 on the odd branch are 3*2^k (k>=4, one chain), 9*2^k (k>=3, one chain),
  27*2, 27*4^, 81*2, 81*4^, ... (covers exactly one class mod 9 x odd) and (3^k, 2^0) k>=4
  (Owens' 81^, 1/108). Whatever classes they are placed in, the open set is one class mod 3 x odd
  (1/6, "R6"), the other class mod 3 x three classes mod 8 (1/8: a class mod 4 "R1" plus a class
  mod 8 "R2", or three R2), one class mod 9 x odd minus 81^ (5/108, "R5"), one class mod 9 x a
  class mod 4 (1/36, "R3"). Total 0.3657. Unions across these (Nielsen's W = 1 mod 3 u 3 mod 9)
  are the only way to merge regions.
* At 45 the only lost resources are 11*{1,2,3,4}, 13*{1,2,3}, 17*{1,2}, 23*1, 31*1; only 11*4 is
  a loss relative to 42.

Set counts at 45 (base 2-adic + 3-adic sets, then 5^, 7^, 11^, 13^, 17^ sub-arrows, Nielsen-style;
crude counts are 2-3 below what Nielsen achieves with composites and x's, so "fits" below means
fits with margin or is Nielsen's own count):

| region (density) | 11 (10) | 13 (12) | 17 (16) | 23 (22) | 31 (30) |
|---|---|---|---|---|---|
| R1 class mod 3 x class mod 4 (1/12) | 10 tight | 12 | 16-17 (Owens) | ~16-19 | ~17 |
| R2 class mod 3 x class mod 8 (1/24) | ~13 | ~15 | ~17 | 22 (Nielsen) | ~20 |
| R3 class mod 9 x class mod 4 (1/36) | fits | fits | fits | fits | 30-31 (Owens) |
| R5 class mod 9 x odd (1/18) | ~10 tight | ~12 tight | ~14 | ~16 | ~16 |
| R6 class mod 3 x odd (1/6) | ~5 | ~6 | ~7 | ~8 | ~8 |
| W x class mod 4 (1/9) | 9 (42: 10) | 12 only with 11^ copies at 143 borrowing from an 11-arrow on the other half | ~14 | - | - |
| W x class mod 8 (1/18) | 10 (fits) | 12 | ~15 | - | - |
| R1 u R2 = class mod 3 x 3 classes mod 8 | ~6 | ~7 | ~8 | ~8 | ~9 |
| R3 u R5 (two classes mod 9) | - | - | - | - | ~18 |

Consequences:

* R6 needs two primes (no prime covers a class mod 3 x odd), R1 u R2 needs two primes, R3 and R5
  one each: six regions, five primes, unless W-type unions merge R5 into the two halves of R6.
* The W merge needs an 11-arrow on a W half (9 sets at 45) or a 13-arrow whose two 11^ copies
  at 143 borrow from such an 11-arrow. Fraction supply at the 143 level, counted over both
  halves, is ~18.6 units; the two copies need ~13 and a 13^ replacing the lost 11-input needs
  ~9: 22 > 18.6. A 13-arrow on a W half without an 11 partner gets 11 sets (one full copy).
* 11 on W x 1 mod 8 keeps 10 sets, but then W x 5 mod 8 needs 16 sets from 17 (~15 with
  borrowing: the nearest miss) and hole 12 (R1) is left with 23 (needs 22, has ~19).
* 17 on a whole W half needs 16, has ~14 (11^ and 13^ copies at 187/221 compete for the same
  resources from both sides). Even 19 (if it were free) counts ~14-15 on a W half.
* R1-type regions are coverable only by 11, 13, 17 (19); at 45 there are three of them plus R5,
  for three primes.

So no Nielsen-style odd branch exists at 45 with Owens' even branch, which consumes 19, 29, 37,
..., 89 exactly. The remaining avenue is a new even-branch seed (prime-5/7 stages) that frees one
top-level prime, preferably 19, for the odd branch, or a borrowing structure beyond Nielsen's
doubling. Both are new constructions, not modifications.

## 21 Oct (later): computational tools and the density-budget argument

Tools (scratchpad m45/): `engine.py` (exact CP-SAT: max resource-disjoint full covers of a region on the
(2,3,5,7) kernel; too slow beyond ~10k rows), `greedy.py` (constructive min-waste per-cover ILP with
restarts: lower bounds only, finds 9 of Nielsen's 10 at depths (4,4,2,1)), `ubound.py` (rigorous
fraction upper bound: every full cover needs total resource fraction >= 1 - h, each exponent tuple is
one modulus, q^ sub-arrows at level P*q need q-1 full covers there; NO borrowing from other arrows).
Tables: `notes/fable-m45-fraction-bounds-M45.txt`, `...-M42.txt` (format P:bound/needed).

Calibration at 42: the bound reproduces Nielsen's 11 on W x 1 mod 4 (bound 11, actual 10) and 23 on
hole 24 (22/22), and undercounts by 1-3 exactly where Nielsen/Owens borrow heavily (13 on a W half with
the two 11^ copies: 11/12; 17 on hole 12: 15/16; 19 on hole 4 x class mod 5: 17/18). So "x" in the
table means "impossible without substantial borrowing from another top-level arrow on an overlapping
region".

Maximal region density per prime at 45 (regions on which the prime fits, with Nielsen-level borrowing
allowed): 11 -> 1/12 (class mod 3 x class mod 4; a W half of density 1/9 only via a 13^ copy at the
143 level, which collides with the 13-arrow's own copies: 143-level supply 18.6 units over both halves
vs demand >= 15.5 + a 12th input for 13 that does not exist); 13 -> 1/12 (1/9 needs an 11 partner on
the other W half, impossible at 45); 17 -> 1/12 (hole-12 type with borrowing) or 1/18 (W x class mod 8,
bound exactly 16); 23 -> 1/24; 31 -> 1/36. Sum 1/12+1/12+1/12+1/24+1/36 = 0.319 (or 0.292) against the
odd branch's open density 0.366 after the forced {2,3} frame. At 42 the same sum is
1/9+1/9+1/12+1/24+1/36 = 0.366: Owens is exactly at the budget. The deficit at 45 is 0.046 = the
class-mod-9 region R5, i.e. the W trick, and no geometry of the prime-3 stage changes it (the open set
after the {2,3} resources >= 45 is forced up to relabelling: a class mod 3 x odd, three classes mod 8
in the other class mod 3, a class mod 9 x odd minus 81^, a class mod 9 x a class mod 4).

Other sources of coverage checked and rejected: fresh primes as arrows on progression pieces
(capacity 0.0165 for 89..6000, consumed by the 42/43 holes); fresh primes 11*p*s as plain classes
(each covers ~12/p of the K44 hole but the leftover never closes: a sieve, and closing it needs
arrows on astronomically many pieces); partially used levels 11*p for Owens' p=17..89 (plain classes
only, same problem; an arrow there needs p-1 full sets where ~4-7 units are free); level>=2 spares
of the 11-arrow (27.5 % of anything in 1 mod 11, exact); 169-level spares (<2 %); moving 5/7 frame
resources (7*8^, 5*16^) to the odd branch to pre-cover a 7- or 5-class of W x 1 mod 4 (gains at most
one set for the 11-arrow and costs the even branch 1/28 or 1/20, which its primes cannot absorb:
19, 29, 37..47 are at their own maximal densities 1/20, 1/30, 1/60).

Conclusion of this round: within the Nielsen/Owens framework (2^ frame, prime-3 stage on the odd
class, 5/7 stages on the even holes, top-level primes on product regions with doubling/borrowing)
minimum modulus 45 is impossible by the density budget, and the patch of the 44-system is the best
possible (short by about 0.005 of density: K44 after spares 0.0037 + secondary holes 0.001 vs fresh
capacity 0.00085). A 45 would need a construction outside this framework.

## 21 Oct: constructive counts (greedy engine, depths (4,4,2,1), lower bounds) and realized patch capacity

Pipeline results so far (sets found / needed): 11 on W x 1 mod 4: 9/10 at 42 (Nielsen 10), 7/10 at 45;
17 on W x 5 mod 8: 12/16; 17 on hole 12: 12/16 (Owens 16 with borrowing); 23 on hole 12: 13/22;
11 on R5 = (3 mod 9) x odd: **10/10 (fits, constructive)**; 13 on R5: **12/12**; 17 on R5: 14/16;
19 on R5: 13/18; N(R) (all tuples) W4 11, W8 14, R1 14. The greedy undercounts by 1-3 relative to
Nielsen-level designs, so these are consistent with the fraction bounds above.

Realized fresh-prime capacity in the 43/44 patches (from the certificates): the 42-hole H has density
1/42 = 0.0238; the seed covers 29.5 %, the 14 smooth spares bring it to 75.6 %, so the 31 fresh
primes 89..257 delivered only 0.0058 of density, and the 12 primes 263..317 delivered 0.00023 on K.
That is 0.006 for 43 primes, well below the shape-optimal model (0.0157): the hole shapes force deeper
pieces. For 45 the remaining need is 0.0047 (K44 after spares 0.0037 + secondary 0.001) with primes
>= 331 only (model 0.00085, realized far less): short by an order of magnitude.

Best redesign configuration found (odd branch): 11 on U1 = (1 mod 3) x 1 mod 4 u (3 mod 9) x 1 mod 8
(a depth-3 region with 10 sets plus slack, if the union keeps >= 10: estimated 8-10), 13 on W x 3 mod 4
with the two 11^ copies, 17/23/31 on holes 12/24/36 as Owens. Orphan: (3 mod 9) x 5 mod 8, density
1/72 = 0.0139, a product region supporting any prime up to 37 but all of 11..31 are used; as fresh-prime
pieces it would need about 0.014-0.018 of capacity on top of H'. Not reachable.

Final pipeline results (label: found/asked): W1_42 9/10, W1_45 7/10, N_W4 11/24, W8_17_45 12/16, N_W8 14/24, H12_17_45 12/16, H12_23_45 13/22, N_R1 14/28, R5_11_45 10/10, R5_13_45 12/12, R5_17_45 14/16, R5_19_45 13/18, N_R5 16/24, H24_23_45 16/22, N_R2 17/28, H36_31_45 24/30, N_R3 24/34. Raw solutions: data/fable-m45-setcount-pipeline.json.
