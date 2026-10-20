# Minimum modulus of distinct covering systems

Computer-assisted upper bounds for the minimum modulus of a finite distinct covering
system (BBMST distortion method + convex comparison), and explicit constructions on the
lower side. Split off from `covering-weight-mu` on 20 October 2026; that repository keeps
the investigation of BBMST's multiplicative weight.

## Current state (20 October 2026)

* **Upper bound 14,484** (`data/latest-certified-bound.json`, proof
  `src/conditional-stoploss-improvement.tex` Section "standalone certified cutoff 14485",
  witness `data/sparse-certificate-14485.json`, verification `data/sparse-verified-14485.json`).
  The sections below record the earlier certificates 14,496, 14,599, 16,999, 19,999,
  39,999 and 615,299 in the order they were obtained. The target 10,000 remains unproved;
  `notes/fable-notes-2026-10-20-upper.md` explains why the comparison method is
  saturated and what lemma would be needed.
* **Lower side.** Owens' record 42 is the published state of the art. Candidate systems
  with minimum modulus 43 and 44, built by patching Owens' construction with fresh-prime
  arrows, are recorded in `data/fable-m43-construction.json` and
  `data/fable-m44-construction.json` (verified on the finite kernel and by integer
  sampling; they inherit the unverified correctness of Owens' thesis Sections 3.5-3.20).
  Minimum modulus 45 is infeasible for that patch scheme. Codex's searches at 43 (direct),
  50, 60, 100 and 1000 are in `src/construction-search.tex` and the sections below.

Python environment for the `code/fable_*.py` tools: `tmp/fable-venv` (numpy, ortools,
python-sat). Everything else uses the standard library unless a section says otherwise.

## Latest certified bound: 14,496 (15 October 2026)

**Every finite distinct covering system has a modulus at most 14,496.**
The target 10,000 remains unproved.

- `src/conditional-stoploss-improvement.tex`: proof, updated in place.
- `data/latest-certified-bound.json`: strongest completed result.
- `data/prefix-refinement-certificate-14497.json`: manifest of the three component calculations.
- `data/prefix-refinement-verified-14497.json`: exact composition of their completed bounds.
- `data/sparse-verified-14498.json`: completed 250,000-prime base verification.
- `data/prefix-base-14498-verified.json` and `data/prefix-fine-14497-verified.json`: completed prefix checks.

The new bound conditions the restricted-loss complement directly and uses
different moment orders on different terminal product bins. Every omitted
configuration has an explicit bound. A finer verified prefix is combined with
the unchanged base suffix and tail. At 250,000 primes it certifies:

    finite survivor lower bound > 0.0001526084818670576
    all remaining losses        < 0.0001478123106801942
    infinite survivor margin    > 0.0000047961711868634

Reproduce with ordinary Python 3:

    python code/check_conditioning_methods.py
    python code/verify_prefix_refinement.py data/prefix-refinement-certificate-14497.json --output data/prefix-refinement-verified-14497.json

The search uses NumPy; verification uses only the standard library.
The default reproduction command reruns all three components. The delivered
composition report was formed after their separate successful runs; its
`components_replayed_in_this_run: false` refers only to the final composition
step. A composition-only run does not independently establish its input reports.
Experiments and failed searches are recorded separately in
`data/conditioning-experiment-summary.json`. A failed local search is not
an impossibility result. The proof is not independently peer reviewed.

The preceding completed bounds of 14,500 and 14,579 remain in
`data/sparse-verified-14501.json` and `data/sparse-verified-14580.json`.
The unfinished standalone 14,499 check was stopped after the stronger
14,496 certificate completed; no standalone result file is claimed for it.

The proof now implements finite recursive joint atom bounds, including an
additional quadratic density envelope. Local moment enclosures and the new
experiments are recorded in `data/recursive-atom-experiment-summary.json`.
These preserve dependencies between primes; marginal bounds are never
multiplied. The extended optimizer also allows nonzero distortions at primes
2 and 3. The canonical result index changes only after full verification.

## Earlier certified bound: 14,599 (15 October 2026)

**Every finite distinct covering system has a modulus at most 14,599.**
The target 10,000 remains unproved.

- `src/conditional-stoploss-improvement.tex`: mathematical proof, updated in place.
- `data/latest-certified-bound.json`: index of the strongest completed certificate.
- `data/sparse-certificate-14600.json`: exact rational proof inputs.
- `data/sparse-verified-14600.json`: completed integer verification and witness hash.
- `code/verify_sparse_conditioning.py`: standard-library verifier.

The new certificate uses exact early conditioning, jointly optimized early
and later distortions, full-product conditioning, and an infinite-tail moment
of order 11/4. It explicitly retains 85,410 configurations of eight primes
and integrates the remaining product jointly with forbidden divisors.
At 125,000 primes, ending at 1,655,131, outward integer arithmetic certifies:

    finite survivor lower bound > 0.0010344360279985
    all remaining losses        < 0.0008859841856684
    infinite survivor margin    > 0.0001484518423301

The positive `infinite_survivor_lower` certifies noncoverage. The negative
`bbmst_gap` is an unused stopping comparison.

Reproduce with ordinary Python 3:

    python code/check_conditioning_methods.py
    python code/verify_sparse_conditioning.py data/sparse-certificate-14600.json --output data/sparse-verified-14600.json

Search requires NumPy; proof verification uses only the standard library.
The full-product correction is tested against independently enumerated
rational distributions, including higher prime powers. The proof explains
why a smaller explicitly truncated prime set can retain more correction.
Experimental results are in `data/conditioning-experiment-summary.json`.
The intermediate bound 14,799 is also certified in
`data/sparse-verified-14800.json`. All earlier witnesses are preserved.
No global optimum, proof-assistant formalization, independent peer review,
or publication priority is claimed.

## Earlier certified bound: 16,999 (15 October 2026)

**Every finite distinct covering system has a modulus at most 16,999.**
The target 10,000 remains unproved.

The principal result files are:

- `src/conditional-stoploss-improvement.tex`: the complete mathematical argument, updated in place.
- `data/sparse-verified-17000.json`: the completed integer verification, endpoints, and witness hash.
- `data/sparse-certificate-17000.json`: the exact rational parameters and proof settings.
- `code/verify_sparse_conditioning.py`: the reproducible standard-library verifier.

The calculation conditions on 891,573 configurations of primes through 53.
It uses sharper product-bin supports and parameters optimized over the whole
finite horizon. An explicit third-moment tail bound closes the infinite
calculation, using Dusart's prime-counting inequality. At stage 50,000,
prime 611,953, integer arithmetic certifies:

```
finite survivor lower bound > 0.0037894992424377
all remaining losses        < 0.0027157104637365
infinite survivor margin    > 0.0010737887787013
```

The BBMST second-moment stopping comparison does not certify this endpoint;
the negative `bbmst_gap` in the output records that unused comparison.
The positive `infinite_survivor_lower` is the successful certificate.

Reproduce it with ordinary Python 3 (about seven minutes locally):

```powershell
python code/test_reachable_support.py
python code/test_third_moment_tail.py
python code/test_normalized_conditioning.py
python code/verify_sparse_conditioning.py data/sparse-certificate-17000.json --output data/sparse-verified-17000.json
```

The report includes the full-horizon adjoint equations, NumPy search commands,
and unsuccessful experiments at 10,000 and 16,900. Optimization is exploratory;
the rational witness is independently checked with outward integer arithmetic.
No global optimum or proof-assistant formalization is claimed. An intermediate
complete certificate establishes 17,999 in `data/sparse-verified-18000.json`.

## Earlier certified bound: 19,999 (14 October 2026)

**Every finite distinct covering system has a modulus at most 19,999.**
Distinct moduli all at least 20,000 cannot cover the integers. The target
10,000 remains unproved; unsuccessful and incomplete experiments at smaller
cutoffs are kept separate from the integer certificates.

The existing report `src/conditional-stoploss-improvement.tex` now contains
the full mathematical reduction. The improvement conditions on 64,315
configurations of primes through 53, keeps an unrestricted excess bound on
the other configurations, and tracks true first moments inside the product
bins. A finite policy search improves the distortion parameters. BBMST's
distortion construction, generic second-moment product, and analytic stopping
theorem remain in place.

The exact witness is `data/sparse-certificate-20000-50000.json`; its verified
endpoints and SHA-256 hash are in `data/sparse-verified-20000.json`.
At stage 50,000, prime 611,953, integer arithmetic certifies:

```
survivor mass > 0.00982698529
f_k < 4,530,040.397
analytic stopping threshold > 5,203,177.214
certified gap > 673,136.818
```

Reproduce the certificate with ordinary Python 3, using only the standard
library (about seven minutes locally):

```powershell
python code/test_minimum_modulus.py
python code/test_stoploss.py
python code/test_strong_conditioning.py
python code/test_moment_bins.py
python code/test_sparse_conditioning.py
python code/verify_sparse_conditioning.py data/sparse-certificate-20000-50000.json --output data/sparse-verified-20000.json
```

The new independent checks include 48 generalized conditioning cases,
903 exact first-moment bin comparisons, finite path support checks, and
72 exact bounded-correction comparisons. The earlier arithmetic, moment,
and convex-comparison checks also pass. This is a computer-assisted proof
with a written reduction, not a proof-assistant formalization or a claim
of independent peer review or publication priority.

To generate a fresh experimental candidate, use Python with NumPy:

```powershell
python code/sparse_conditioning_experiment.py --K 20000 --start 16 --small 16 --state-cap 3000000 --policy 0.9 --min-stop 50000 --output data/sparse-regenerated-20000.json
python code/verify_sparse_conditioning.py data/sparse-regenerated-20000.json --output data/sparse-regenerated-verified.json
```

Search results are not proofs. `data/sparse-policy-sweep.json` records the
finite policy comparison; the selected policy is 0.9. Trials at 10,000
still exhausted the survivor estimate. Some 18,000 and 19,000 runs retained
positive estimates but did not satisfy the stopping condition at the tested
endpoint. None is claimed as a certificate or an impossibility result.
Intermediate complete certificates also establish 29,999 and 21,999.

## Earlier certified bound: 39,999 (14 October 2026)

**Every finite distinct covering system has a modulus at most 39,999.**
The new certificate excludes distinct moduli all at least 40,000, achieving
the requested scale of improvement. The 615,299 certificate below is retained
as an earlier result.

The proof is in `src/conditional-stoploss-improvement.tex`. It retains BBMST's
distortion measures and analytic stopping theorem, but replaces the quadratic
loss estimate with a direct convex-comparison bound. Conditioning on the
comparison valuations at 2 and 3 preserves part of the minimum-modulus
restriction. A finite product distribution is enclosed with downward bins,
exact Euler-product means, and explicit probability and weighted-mean tails.
The report proves the comparison and the directions of all approximations.

The rational witness is `data/stoploss-certificate-40000.json`; verified
endpoints, tail bounds and its SHA-256 hash are in
`data/stoploss-verified-40000.json`. At k = 50,000 (prime 611,953), the checker
establishes survivor mass greater than 0.0093593 and

```
f_k < 5,053,132.141
analytic stopping threshold > 5,203,177.214
certified gap > 150,045
```

This is a computer-assisted proof with a written mathematical reduction,
not a claim of independent peer review, publication priority, or global
optimality. The arithmetic verifier uses only Python's standard library:

```powershell
python code/test_minimum_modulus.py
python code/test_stoploss.py
python code/verify_stoploss.py
```

The full verification took approximately 254 seconds locally. All numerical
enclosures use integers at scale 10^40. Tests include 361 independent rational
stop-loss checks, 384 dependent-coordinate convex comparisons, the earlier
2,520 interval checks, and 24 restricted-moment checks.

`code/conditional_stoploss_experiment.py` is a separate, nonrigorous NumPy
search. It uses systematic local parameter optimization, not a globally
optimal control claim. With NumPy available, regenerate and then verify with:

```powershell
python code/conditional_stoploss_experiment.py --K 40000 --start 20 --cap 100000 --resolution 256 --min-stop 50000 --output data/stoploss-certificate-40000.json
python code/verify_stoploss.py
```

The supplied witness fixes the exact rational parameters; regenerated
floating-point candidates must pass the verifier before being used as proofs.
This earlier certificate is superseded by the 20,000 cutoff above.

## Certified minimum-modulus improvement (14 October 2026)

**Every finite distinct covering system has a modulus at most 615,299.**
Equivalently, distinct moduli all at least 615,300 cannot cover the integers.
This is a computer-assisted deduction using BBMST's existing distortion method
and Theorem 6.1, with an integer-arithmetic certificate. It is not a claim of
independent peer review, publication priority, or an optimal numerical bound.

The complete mathematical reduction is in `src/minimum-modulus-improvement.tex`
(successfully compiled in the Codex LaTeX editor). The rational distortion
parameters are in `data/minimum-modulus-certificate.json.gz`; the certified
endpoints and witness hash are in `data/minimum-modulus-certificate-verified.json`.
At stage 1,050,000 (prime 16,313,527), the certified upper bound for f is below
191,101,044, whereas BBMST's analytic stopping threshold is above 191,181,612.
That theorem handles all subsequent primes.

Run from the repository root, using Python 3's standard library only:

```powershell
python code/test_minimum_modulus.py
python code/verify_minimum_modulus.py
```

The checker uses outward-rounded integer intervals with 40 decimal places;
neither floating-point moment estimates nor rounded published stopping tables
are trusted. Exact-rational tests independently compare the moment reduction
with BBMST's recursive Theta formula (24 cases), alongside 2,520 interval checks.
Local certificate verification took about 32 seconds.

The experiment script, `code/minimum_modulus_experiment.py`, is explicitly
nonrigorous. `code/optimize_minimum_modulus.py` performs a reproducible six-block
parameter search. To regenerate a candidate and certify it:

```powershell
python code/optimize_minimum_modulus.py
python code/minimum_modulus_experiment.py --K 615300 --blocks data/minimum-modulus-optimization.json --limit 18000000 --min-stop 1050000 --output data/minimum-modulus-certificate.json.gz
python code/verify_minimum_modulus.py
```

The improvement comes from parameter selection and explicit numerical
continuation. A gcd convolution evaluates the restricted second-moment sum
efficiently and rigorously; it does not establish a stronger universal moment
inequality. Search failures at smaller cutoffs are not impossibility results.
The earlier weight investigations below are separate and remain unchanged.

## Construction search toward minimum modulus 43 (15 October 2026)

**No construction improving Owens's 42 record was found in this investigation.**
This direction is separate from the upper-bound work above. The working note
`src/construction-search.tex` gives a finite prime-power-tail search model and
proves that a successful kernel certificate expands to a finite covering with
distinct moduli. It compiled successfully in the desktop LaTeX editor.

`code/construction_kernel.py` is an exact standard-library verifier. Positive
controls of minimum moduli 3 and 6 are in
`data/construction-kernel-control-3.json` and
`data/construction-kernel-control-6.json`. They are validation examples, not
record improvements. All 240 comparisons against direct enumeration and 51
full-period finite-tail tests pass:

```powershell
python code/test_construction_kernel.py
python code/construction_kernel.py data/construction-kernel-control-6.json
```

The searches at 43 with base primes through 13 and 19 did not produce a
verified cover. Sample-feasible candidates retain explicit uncovered states;
time limits do not prove impossibility. Results are recorded separately in
`data/construction-kernel-43-13.json` and
`data/construction-kernel-43-19.json`. The supplementary minimum-10 search also
did not certify a cover. These are initial bounded experiments, not an
exhaustive search of the construction family.

The search needs OR-Tools (installed locally under the ignored
`tmp/construction-deps` directory during this session):

```powershell
python -m pip install --target tmp/construction-deps ortools
python code/search_construction_kernel.py --seconds 180 --initial 128 --output data/construction-kernel-43-13.json
python code/search_construction_kernel.py --primes 2,3,5,7,11,13,17,19 --depths 6,4,3,2,2,2,2,2 --minimum 43 --cap 1000 --seconds 240 --initial 1024 --batch 64 --seed 43 --output data/construction-kernel-43-19.json
```

The stored initial runs used OR-Tools 9.15.6755, four workers, and the same
seed on successive solves. The current script varies the seed on successive
solves; solver versions, thread scheduling, and seeds can change search paths.
Successful certificates must always pass the independent verifier.

Owens's full thesis was subsequently obtained through the institutional page's
browser download and saved as `literature/owens-2014.pdf`. Sections 3.1-3.4
were audited for modulus-family reservations, and the later dependencies were
read. The full record-42 construction is not yet independently verified.

The expanded working note locates the class `10 mod 42` in Section 3.4. It
proves that the conservatively overestimated pool of spare 7-smooth moduli has
conditional capacity only `2759/3600` for covering this whole class. This
excludes a restricted whole-branch replacement, not an overlap-assisted repair
or a rearrangement of Owens's cover. Even releasing the `36*7^d` family gives
only `1153/1200`.

The fresh-prime search gives six finite partial-cover schemas, covering
`10,52,94,136,178,220 mod 840`, with fresh root primes 97,101,103,107,109,113.
Their modulus families are distinct and avoid Owens's symbolic families.
They cover 3/10 of the whole `10 mod 42` branch, not necessarily 3/10 of its
actual uncovered subset. They do **not** form a full repair or a new record.
The note proves the construction rule and explains finite completion.

```powershell
python code/owens_branch_search.py --iterations 2500
python code/test_owens_branch_search.py
```

`data/owens-42-branch-search.json` contains the six checked schemas, exact
density sums, and 72 uniform-split experiments. A rational relaxation rules
out all 72 within the stated frozen-base counting rule; these are restricted
exclusions, not general impossibility results. A joint encoding of the true
residual gap, including later-stage overlaps and early prime-3 rearrangements,
remains to be built.

## Construction target 1000 (15 October 2026)

The new target is a finite distinct covering system with every modulus at
least 1000. **No such construction, or improvement over 42, was obtained.**
The working note now includes a separate target-1000 section.

Exact capacity screening shows that a target-1000 cover needs some prime
divisor at least 59: the reciprocal sum of all 53-smooth moduli at least
1000 is only 0.9843755046..., while the corresponding sum through prime 59
is 1.0537687824.... A separate elementary harmonic-sum bound requires at
least 1718 classes. These necessary conditions do not prove existence.
The Owens audit also identifies 106 allocated early classes below 1000
that would need changing; the whole original cover remains unaudited.

The resource generator now prunes by the base cap, and the exact verifier
uses a sparse digit tree to partition equivalent coordinate states. The
search can handle base-prime supports through 311 without scanning the
full Cartesian product. Tests include 160 direct coordinate-partition
comparisons, 240 complete small-kernel comparisons, 30 resource-generation
comparisons, and 51 finite-tail checks.

One 180-second CP-SAT run and three sample-guided candidate runs at minimum
1000 failed exact verification. The larger sample-guided runs use 262144
training points and primes through 149 and 311. Each stored candidate has
an explicit uncovered progression independently checked by modular
arithmetic, including incompatibility with its finite tail completion.
Training-sample coverage is not a covering certificate or evidence of a
new record. A small minimum-3 control is verified; the additional
minimum-10 diagnostic is rejected.

```powershell
python code/search_large_minimum.py --profile --output data/construction-1000-profile.json
python code/search_large_minimum.py --largest-prime 311 --samples 262144 --sweeps 6 --seconds 180 --seed 1002 --uniform-weights --output data/construction-greedy-1000-311.json
python code/test_construction_kernel.py
python code/verify_construction_searches.py data/construction-kernel-1000-149.json data/construction-greedy-1000-149.json data/construction-greedy-1000-149-large.json data/construction-greedy-1000-311.json data/construction-greedy-control-10.json --output data/construction-1000-witness-checks.json
```

The profile and witness checks use the standard library. The sample search
needs NumPy, available in `tmp/construction-deps` in this session; CP-SAT
also needs OR-Tools. The earlier two weighted sample runs predate the
explicit `weight_policy` metadata field; both multiply uncovered-point
weights by 1.3 after each sweep, as the current default does.

## Structural construction target 100 (15 October 2026)

The latest continuation uses **no sampled states**. An exact decision
diagram represents the uncovered kernel, and each placement computes its
new coverage and all overlaps with integer arithmetic. A verified
minimum-10 positive control covers 91,968,307,200 kernel states using 215
resources. Its finite completion is specified in
`data/construction-structural-control-10-verified.json`; this is a machinery
test, not a record improvement.

For minimum 100, the strongest base run here places 257 rows on primes
through 31 before reaching its diagram limit. Its exact uncovered kernel
fraction is 0.391747707953315.... This is a partial construction, not an
impossibility result. Alternative exact placement runs also fail, with
independently checked uncovered progressions.

The new structural repair assigns a separate fresh prime `q` to each block.
Choose `q-1` distinct divisors `n_d` with `q*n_d >= 100`, and put the
congruence `a mod n_d` on root tail digit `d`. Distinct roots prevent modulus
collisions. All root states are covered exactly when the base state lies
in `a mod lcm(n_d)`. This identity permits exact overlap accounting and an
explicit CRT witness whenever a projected gap remains.

Nineteen such blocks, rooted at primes 37 through 113, add 1414 rows. The
fraction of base states admitting an uncovered extension decreases to
0.389603370030882.... This is a projected kernel fraction, not the natural
density of uncovered integers. The full 1671-row candidate is rigorously
rejected: its saved CRT progression misses every row and every finite tail
completion. The block schemas are checked separately.

There is also a rigorous obstruction to completing this particular repair.
With the 257-row base frozen, even allowing all original base primes in
prefix divisors and every admissible fresh root, **one independent divisor
block per root** can remove at most
`9001713966181/3032474610765600 = 0.002968438361932...` of the base-state
space. The audit exhausts 40320 prefix moduli and 4220 roots. It excludes
this restricted repair family, not minimum 100 generally. Changing the
base, using multiple interacting blocks at a root, or making later blocks
depend on earlier root primes remains outside this obstruction.

The proof and limitations are in the new target-100 section of
`src/construction-search.tex`. **No minimum-100 cover or improvement over
42 has been obtained.** Reproduce the main experiment and audit with:

```powershell
python code/structural_cover_search.py --minimum 100 --largest-prime 31 --cap 1000 --seconds 120 --nodes 200000 --order prime-stage --placement digits --output data/construction-structural-100-31-digits.json
python code/structural_prime_patches.py data/construction-structural-100-31-digits.json --seconds 120 --maximum-root 113 --output data/construction-structural-100-prime-patches.json
python code/structural_prime_patches.py data/construction-structural-100-31-digits.json --capacity-audit --output data/construction-structural-100-patch-capacity.json
python code/verify_construction_searches.py data/construction-structural-100-prime-patches.json --output data/construction-structural-patch-witness-check.json
python code/test_structural_cover_search.py
python code/test_structural_prime_patches.py
```

All these commands use the standard library. The new tests include 900
exhaustive marginal/residual comparisons, 2880 projected-state checks,
14 small capacity audits, and row round-trips in three coordinate orders.
The earlier finite-kernel and finite-tail tests also pass. Time and node
limits are recorded in the artifacts; a limit is never reported as proof
of coverage or impossibility.

## Construction target 60 (15 October 2026)

Five exact, unsampled searches were run at minimum 60. **No complete cover
was found.** The density screen first permits the prime range through 13
(capacity 1.0356386372...), while primes through 11 are insufficient
(capacity 0.7948533807...). These are necessary conditions only. At least
103 distinct classes are required by the elementary reciprocal-sum bound.

The best partial candidate in this batch places 260 rows on primes through
31 and leaves an exact kernel fraction of 0.271856093967964... uncovered.
It stopped at the diagram limit. Its saved uncovered progression has been
checked against every row and its finite completion. All four other
candidates also have independently verified failure witnesses in
`data/construction-60-witness-checks.json`.

The new `disjoint` placement rule universally quantifies all unselected
digits and accepts only whole cylinders inside the residual. This controls
overlap by forbidding it within the kernel. It examined all 2098 allowed
resources through prime 43 with base cap 10000, placed 180 rows, and left
74.5026866...% of its kernel uncovered. All 16110 pairs of its cylinders
were independently checked for disjointness, and their total mass checked
against the reported residual. Finite closing classes may still overlap.
This greedy failure is not an impossibility result for disjoint constructions.

For the frozen 260-row base, an exhaustive audit also excludes repair by
one independent fixed-divisor block per fresh root: the sum of maximum
block masses is only 0.003066838588637... of the base-state space. The audit
includes 27648 prefix moduli and 3005 roots. Interacting blocks and changes
to the base are outside this restricted obstruction.

The working note now has a separate target-60 section. To reproduce the
main searches and the restricted capacity audit:

```powershell
python code/structural_cover_search.py --minimum 60 --largest-prime 31 --cap 1000 --seconds 120 --nodes 200000 --order prime-stage --placement digits --output data/construction-structural-60-31-digits.json
python code/structural_cover_search.py --minimum 60 --largest-prime 43 --cap 10000 --seconds 90 --nodes 100000 --order prime-stage --placement disjoint --layout large-prime-first --output data/construction-structural-60-43-disjoint.json
python code/structural_prime_patches.py data/construction-structural-60-31-digits.json --capacity-audit --output data/construction-structural-60-patch-capacity.json
python code/test_structural_cover_search.py
```

The updated tests compare the new whole-cylinder calculation with literal
enumeration in 900 cases. The evidence so far establishes neither existence
nor impossibility at minimum 60.

## Squarefree construction target 30 (15 October 2026)

**No minimum-30 squarefree cover was obtained.** This restriction changes
the problem substantially. Dalton and Trifonov's *Extreme Covering Systems*
(2022), Problem 4, asks whether every distinct squarefree cover must have
least modulus 2. A minimum-30 construction would answer that problem
negatively. Cummings, Filaseta, and Trifonov prove an upper bound of 118;
that bound does not decide existence at 30. Sources and links are recorded
in the working note's squarefree section and `literature/SOURCES.md`.

The previous prime-power tail model does not preserve squarefreeness, even
if its depth is set to 1. The new search uses ordinary residue coordinates
with p states per prime and literal squarefree moduli. No tail completion
is used. The independent verifier checks prefix rows only.

Exact density screening shows that some modulus must have a prime divisor
at least 37. All squarefree moduli at least 30 supported on primes through
31 have reciprocal sum only `31160796824/33426748355 < 1`.

The bounded search through prime 37 places 121 classes before its diagram
limit. Its missed density is 0.359330739535289..., and the progression
`5045898743977 mod 7420738134810` is directly checked to miss every class.
This rejects the saved candidate, not the target in general. A minimum-2
squarefree control with 14 classes was verified by examining every residue
modulo 210. Tests compare 200 small systems against direct enumeration.

```powershell
python code/search_squarefree_construction.py --minimum 30 --largest-prime 37 --seconds 60 --nodes 100000 --output data/construction-squarefree-30-37.json
python code/search_squarefree_construction.py --minimum 2 --largest-prime 7 --output data/construction-squarefree-control-2.json
python code/test_squarefree_construction.py
```

## Squarefree minimum 3: joint parity construction (15 October 2026)

**No complete minimum-3 cover was found.** The new approach normalizes
the class modulo each odd prime p to residue 0, and the class modulo 2p
to residue 1 modulo p with a free parity. Adding unused resources and
permuting CRT coordinates makes this normalization lossless for existence.
All remaining composite classes are then optimized jointly on the two
parity layers, with every unit-box point constrained, without sampling.

This yields a rigorous counting obstruction: any distinct squarefree
system with all moduli at least 3 and prime factors at most 11 misses at
least 22 residues modulo 2310. Thus a putative cover must use a prime
divisor at least 13. The proof enumerates eight parity patterns and bounds
the maximum residual coverage of each composite resource. It does not
depend on a solver's infeasibility or optimality report.

The joint solver returned partial candidates through primes 7, 11, and 13.
A second structural step removes classes d and 2d together and computes
their best joint replacement exactly. The prime-13 gap decreased from
1921 to 829 residues modulo 30030. Extending to prime 17 and allowing all
old pairs to move produced 126 distinct squarefree congruences, missing
8040 of 510510 residues (about 1.575%). All integers in
`404 mod 510510` remain uncovered. The full period and the separate
coordinate-cylinder verifier both check this failure. Solver time limits
and local stagnation do not establish global optimality.

For a frozen residual R modulo L, adding one fresh prime q can supply at
most `sum_{d|L} max_a |R intersect (a mod d)|` incidences on its q copies.
The final candidate would require `q <= 5093/670 < 19`, so no fresh prime
can finish it while its old classes stay fixed. This excludes only that
extension route, not a joint reconstruction or a general minimum-3 cover.

```powershell
python code/squarefree_parity_search.py --largest-prime 13 --seconds 60 --output data/squarefree-parity-3-13.json
python code/squarefree_parity_repair.py data/squarefree-parity-3-13.json --output data/squarefree-parity-3-13-repaired.json
python code/squarefree_parity_repair.py data/squarefree-parity-3-13-repaired.json --extend-prime 17 --sweeps 10 --output data/squarefree-parity-3-17-repaired.json
python code/squarefree_parity_repair.py --anchor-audit --output data/squarefree-parity-anchor-audit.json
python code/test_squarefree_parity.py
```

The solver requires OR-Tools; the repair and proof audits use the standard
library. Tests include 500 exhaustive pair-optimization comparisons,
100 normalization comparisons, and 31 direct anchor-capacity checks.
The full proof is in the new squarefree-minimum-3 section of the research
note. No novelty claim is made for the derived necessary conditions.

## Active unrestricted target: minimum 50 (15 October 2026)

The current objective is an explicit finite distinct covering with every
modulus at least 50. No such cover has been found or verified.

Two structural methods were added. `nested_resource_exchange.py` computes
the exact best joint placement of two nested cylinder resources. It made
three strict improvements to a 73-resource candidate, which remains
independently rejected. `graded_arrow_blocks.py` constructs disjoint input
covers of a progression, records their minimum moduli, and relocates an
entire group of progression branches together using exact residue counts.
Neither method uses sampled coverage as a certificate.

The graded construction through prime 89 produces 37 guaranteed progressions
on period 151200, leaving 41308 points outside their union. Expanded rows
can cover additional points: this projected gap is not a counterexample
to all expanded congruences. An exact capacity audit proves that the
specific shape and prime-order menu cannot cover the entire projected
period: even its stagewise density maxima sum to only `4975/6048`.

The continuation therefore uses coverage inside the individual inputs of
earlier prime arrows. `borrowed_arrow_probe.py` found three sufficient
relative input counts for a fresh root 97 among 342 tested progressions.
For example, `33 mod 720` has 96 relative input sets, meeting the 96-input
requirement. This witness has now been exported and independently audited
as a relative replacement using root 97. It is a partial replacement,
not a cover of all integers.

```powershell
python code/test_nested_resource_exchange.py
python code/test_graded_arrow_blocks.py
python code/graded_arrow_blocks.py --minimum 50 --largest-prime 89 --core-exponents 5,3,2,1 --sweeps 3 --output data/construction-50-graded-blocks.json
python code/borrowed_arrow_probe.py data/construction-50-graded-blocks.json --root 97 --maximum-modulus 720 --per-shape 4 --output data/construction-50-borrowed-probe-97.json
```

The new tests passed 400 exhaustive pair-placement comparisons and 48
finite-kernel supply audits covering 228 resource-disjoint input sets.
The graded method also passed a complete minimum-3 positive control.
See the unrestricted-minimum-50 section of `src/construction-search.tex`
for the proofs, scope of the obstructions, and remaining verification work.

### Verified relative-arrow continuation

`relative_arrow_search.py` now checks every borrowed input using an
independent CRT calculation and reconstructs the disjoint resource supply.
The frozen-base prime-97 replacement has 29 checked dependencies and reduces
the guaranteed progression gap from 41308 to 41098 modulo 151200.

Rebuilding stages with branch-dependent input counts and preferring broader
admissible inputs gives 36 root stages with 2119 checked dependencies. Its
guaranteed union still leaves 33689 residues modulo 151200. The run allowed
roots through 229; the last successful stage uses 173. Gaps in this projected
union are not claimed as counterexamples to the expanded congruences.

The dependency tests exhaust a 384-state example, check 405 containment
comparisons and 7008 borrowed-coverage implications, and verify 12 disjoint
exponent families. The small example covers all its states.

An exact early-core construction in `exact_core_array.py` places 283
resources using primes 2,3,5,7 and examines all 7257600 kernel states. It
leaves 2998392 uncovered and passes a separate negative-witness check.
It preserves partial tail coverage omitted by the progression-only method.
This core is now integrated using target regions that can fix a tail digit,
with an independent audit of the resulting dependencies (see below).

```powershell
python code/test_relative_arrows.py
python code/relative_arrow_search.py data/construction-50-graded-blocks.json --largest-prime 229 --maximum-modulus 2520 --per-shape 4 --output data/construction-50-relative-229.json
python code/exact_core_array.py --minimum 50 --primes 2,3,5,7 --depths 6,4,3,2 --cap 50000 --order mass --output data/construction-50-exact-core.json
```

Minimum 50 remains the active objective. No complete minimum-50 cover is
claimed by any of these experiments.

### Exact regions and joint replacement inputs

`core_region_arrows.py` integrates the exhaustive core with relative prime
arrows, including chronological residual conditioning and branch-specific
containment checks. A minimum-20 control closes all 241920 kernel states
using 16 stages and 370 checked dependencies. Its saved finite-completion
schema is `data/construction-region-control-20.json`; the full finite list
has not been expanded.

The strongest minimum-50 run so far uses a core concentrated on the even
class, then 26 stages through prime 113. Its guaranteed residual is
812946 of 7257600 kernel states, with 1529 checked dependencies. This is
only a partial schema, not an integer-density claim or a complete cover.
See `data/construction-50-focus-2-regions-113.json`.

The new structural search permits inputs formed from unions of prime-power
classes. `prime_forest_supply.py` packs one-prime antichain covers without
reusing an exponent/input resource. The region audit checks these forests
independently. Initial bounded forest and union-target runs remain partial
and have not improved the best minimum-50 result.

`joint_root_repair.py` goes further by choosing all inputs of one root
jointly, using a shared fixed-exponent resource pool and exhaustive
projection of the target residual. Its first root-11 run has 784 projected
points, 60 resource types, and 4315 candidate classes. The solver returned
UNKNOWN, so it establishes neither a successful replacement nor an
impossibility result. Small positive controls pass direct integer-period
checks, and corrupted certificates are rejected.

A subsequent independently checked integer-weight certificate excludes
that specific 60-resource menu: its capacity is at most
6351876232/1000000179, below the ten required inputs. This explains the
failed attempt without excluding richer constructions. Tail families,
deeper exponents, and other targets remain available. See
`data/construction-50-joint-root-11-capacity.json`.

An independent audit of the explicit Owens-inspired early seed verifies
three local branch claims at cutoff 42. Removing the resources with bases
42, 45, and 48 reopens those branches. This is an audit of the seed, not of
Owens's full thesis. See `data/owens-seed-branch-audit.json`.

```powershell
python code/test_core_region_arrows.py
python code/test_core_region_arrows.py --forest --unions
python code/prime_forest_supply.py
python code/test_joint_root_repair.py
python code/test_owens_seed_core.py
```

The active objective remains an explicit, rigorously verified distinct
covering system with minimum modulus at least 50. It has not been reached.

### Verified local replacement and published-stage audit

The tail-aware joint model finds a root-11 replacement for the whole
`10 mod 42` branch relative to the saved minimum-50 early core. It uses
40 distinct exponent resources in ten inputs; every new modulus is at
least 55. The exact audit checks all 40 projected states in every input,
covering 102240 fine core states. The explicit certificate is
`data/construction-50-joint-tail-repair-42.json`.

This is a local repair: it consumes prime-11 resources required elsewhere
in Owens's construction. Its bounded integrated continuation reaches
root 61 and still leaves 1487500 guaranteed missing states. It does not
improve the best complete-cover record or the best partial run here.

`nielsen_prime11.py` now explicitly reconstructs the published ten-input
prime-11 pattern, with Owens's prime-5 permutation. Its cutoff-42 audit
verifies all 417 exponent resources are distinct and checks every one of
628992 target states in each input. This verifies that local stage only.
At cutoff 50, the exact diagnostic finds the illegal modulus-44 family
and holes caused by removing modulus 45 in inputs 6, 7, 9, and 10.
See `data/nielsen-prime11-transcription.json` and
`data/nielsen-prime11-minimum50-diagnostic.json`.

The joint attempt to repair the three early holes together remains
unresolved. Both explicit class assignment and implicit residue-choice
refinement return no verified replacement; all trial assignments undergo
an exhaustive target check. Tail-aware projection and borrowing tests are
in `code/test_joint_tail_repair.py`.

### Modulus-90 relocation and partial borrowing

A structural core change now restores the odd targets used by the published
prime-11 and prime-13 stages: move the class `78 mod 90` to `3 mod 90`.
The new class restores the odd part of the deleted `3 mod 45` class. The
former even class remains a repair obligation. An exhaustive comparison
of 1612800 kernel states verifies that these odd targets have exactly the
same residual as the cutoff-42 core.

This produces nine verified admissible prime-11 inputs with 116 distinct
resources. Keeping the valid higher powers of its remaining input isolates
the exceptional modulus-44 class instead of discarding the entire family.
These are local certificates, not a complete root or covering system:

- `data/construction-50-relocated90-nine-inputs.json`
- `data/construction-50-relocated90-prime11-truncated.json`
- `data/construction-50-relocated90-odd-agreement.json`

The prime-13 doubling step is now explicitly verified and supported by the
construction checker. It combines partial earlier coverage with new rows
inside each prime-11 digit, without copying borrowed modulus families.
Its compact minimum-50 certificate has 237 own resources and minimum new
modulus 52 (`data/construction-50-relocated90-prime13.json`).

Interleaving useful earlier prime stages with core primes lets a prime-17
stage close the remaining odd target. The bounded minimum-50 continuation
through root 59 is still partial, leaving 1314132 guaranteed missing states.
See `data/construction-50-relocated90-continuation.json`. This does not
improve the best partial gap reported above, and the full objective remains
unresolved.
