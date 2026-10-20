# Source and search record (minimum modulus project)

## Construction target 1000, 2026-10-15

- Read the primary paper Raj Agrawal, Prarthana Bhatia, Kratik Gupta,
  Powers Lamb, Andrew Lott, Alex Rice, Christine Rose Ward,
  *Computations and Observations on Congruence Covering Systems*,
  Integers 24A (2024), A1:
  https://math.colgate.edu/~integers/a1Proc23/a1Proc23.pdf .
  Printed page 2 reports Owens's record 42. This is supporting literature
  context, not an exhaustive current-record verification.
- The exact smooth-number capacity screen, harmonic-sum cardinality bound,
  106-class count from the early Owens reservation audit, new search
  implementation, and candidate rejections are local deductions and
  experiments in `src/construction-search.tex` and the target-1000 data
  files. No source is claimed to establish the desired construction.
- Search results at primes through 149 and 311 are incomplete negative
  experiments, not impossibility theorems. Each retained candidate is
  rejected by an independently checked uncovered arithmetic progression.

## Lower-bound construction search, 2026-10-15

- Tyler Owens, *A Covering System with Minimum Modulus 42*, BYU MS thesis,
  2014: https://scholarsarchive.byu.edu/etd/4329/ . The institutional abstract
  was accessible and confirms 42. After direct requests returned HTTP 403,
  the institutional page's browser Download link succeeded with the user's
  explicit download permission. The file is `literature/owens-2014.pdf`
  (26 PDF pages; SHA256
  `250FA5A3ACBA9D3C08E8064BDED193DE3ECAD1508F3362E38CB652A4E8B007A9`).
  Read the notation, all construction sections, and conclusion. Visually
  inspected PDF pages 14-16 for the 7-stage formulas. Printed pages 8-10
  locate the modulus-42 class, the 7-smooth allocations, and the extra
  `36*7^d` family. The full construction has not been independently verified.
- Pace P. Nielsen, *A covering system whose smallest modulus is 40*,
  J. Number Theory 129 (2009), 640-666,
  https://doi.org/10.1016/j.jnt.2008.09.016 . Read the local copy
  `C:/Users/mekak/OneDrive/Mathmatics/Papers/covering system/Nielsen_2009_minimum-modulus40_web.pdf`,
  particularly Sections 3-5; visually inspected the finite-tail explanation
  on PDF page 5. The symbolic arrow machinery and the warning about simply
  deleting the smallest class come from this source.
- The finite-kernel formulation and fresh-prime completion proof are written
  in `src/construction-search.tex`, using Nielsen's finite-tail idea.
  No novelty claim is made. The method is validated with small positive
  controls; no covering system of minimum modulus at least 43 was obtained.
- The follow-up reservation audit, `2759/3600` conditional-density bound,
  fresh-prime counting search, and six partial branch covers are deductions
  recorded in the same working note, not results attributed to Owens.
  The ambiguity in the compressed placement of the final 125 factors on
  printed pages 8-9 is avoided by omitting those reservations from the
  upper bound on spare capacity. Printed pages 18-19 explicitly suggest
  changing the early prime-3 construction or relocating the extra 81 tail.
- Current-record searches found no verified improvement over 42. This was
  a bounded literature search, not an exhaustive priority check. Later
  expositions and the Erdős Problems page still cite Owens's record.

The previous upper-bound work is separate and was not re-audited in this turn.

## Recursive atoms and first-prime parameters, 2026-10-15

The finite recursive atom calculation and quadratic density envelope are
implemented and proved in the existing proof document. They use an analytic
infinite LCM sum and subtract only explicitly enclosed finite improvements.
Their actual-measure first and second moments are kept separate from the
independent-comparison moments in the stop-loss identity. New exact-rational
and actual-measure checks are recorded by the method test suite.

The further search includes nonzero distortion at primes 2 and 3, a finer
prefix option, explicitly controlled exponent truncation, and a separately
recheckable common future distortion parameter. These are local deductions
from the same BBMST distortion construction and convex comparison; no new
external mathematical input is introduced. The primary BBMST text was
rechecked at https://arxiv.org/html/1811.03547v1 . Numerical discoveries and
verification status are in `data/recursive-atom-experiment-summary.json`.
The 10,000 target remains unproved.

## Conditioned-complement certificate and tail refinements, 2026-10-15

The completed manifest `data/prefix-refinement-certificate-14497.json` and
result `data/prefix-refinement-verified-14497.json` prove a minimum-modulus
upper bound of **14,496**. The survivor margin exceeds 0.0000047961711868634.
This combines the independently verified 250,000-prime base at cutoff
14,498 with an exactly matching old prefix and a finer prefix at 14,497.
The default replay command recomputes all three component checks.
The same proof document has been updated and compiled in place.

The new local argument conditions the restricted-loss complement directly
and bounds omitted configurations using divisor pairing. Further local
work proves dyadic moment bounds through order 5 (subject to explicit
analytic stopping hypotheses) and a tail estimate that
chooses moment orders separately across retained product bins, with an
explicit bound for every discarded path. The search now differentiates
that stratified terminal objective and also tests coordinate and block
updates. These mathematical deductions use the same BBMST comparison and
Dusart prime-counting input already cited below. BBMST's primary text was
rechecked at https://arxiv.org/html/1811.03547v1 .

A further local deduction bounds the fibre density by
1/(1-min(alpha,delta)) <= min(nu,1+nu*alpha), giving a recursive
least-common-multiple atom bound and an integer-moment refinement.
In that earlier run, only its marginal diagnostic was evaluated; it was
not inserted into the independent-product certificate. Multiplying the marginal caps would
not be justified. This is a local deduction, not a new external premise.

Other numerical candidates remain experimental until their complete
integer verification succeeds. No 10,000 result or global
optimization theorem is claimed.

## Full-product conditioning certificate, 2026-10-15

An earlier completed local certificate gives minimum modulus at most
**14,599**, excluding distinct moduli all at least 14,600. Its exact inputs
and completed result are `data/sparse-certificate-14600.json` and
`data/sparse-verified-14600.json`. The latest-result index now points to the
stronger completed certificate. The same proof document was updated and
compiled successfully.

The new full-product correction tracks the remaining divisor-count product
jointly with forbidden-divisor contributions. Its complementary expectation
has explicit finite support. This permits retaining 85,410 configurations
of only eight primes, while integrating the remaining factors jointly.
The proof also establishes monotonicity of the exact retained correction
when the explicitly truncated prime set is reduced. Together with jointly
optimized early/later parameters and the 11/4-moment tail, this yields
at prime 1,655,131 a finite survivor greater than 0.0010344360279985 and
an infinite loss below 0.0008859841856684. The certified positive margin
exceeds 0.0001484518423301.

The full-product implementation passes 108 independent exact-rational
checks, including higher prime powers and zero divisor coefficients.
Four-prime conditioning and the 14,500 target remain separate experiments.
No 10,000 result or global optimization claim is made. The mathematical
external inputs are unchanged; all further inequalities are proved locally.

## Exact early conditioning and fractional moments, 2026-10-15

A completed local certificate improves the minimum-modulus upper bound to
**14,799**. Its witness and outward-integer result are
`data/sparse-certificate-14800.json` and `data/sparse-verified-14800.json`.
The same proof document, `src/conditional-stoploss-improvement.tex`, now
contains the finite-support proof for exact early complements, a statewise
maximum of conditional corrections, the shifted-log prime tail, and the
fractional-moment extension. The certificate uses order 11/4 and has a
certified infinite-survivor margin greater than 0.0003138744607032.

The early and later parameters are varied jointly, including derivatives
of the conditioning probabilities. Additional full-product conditioning
and a smaller explicitly enumerated prime set are explored separately.
Independent checks and source hashes are recorded in
`data/conditioning-method-checks.json`; numerical searches are not certificates.

The external inputs remain BBMST's distortion construction and Dusart's
Theorem 6.9, inequality (6.6), checked again in the primary papers:
https://arxiv.org/html/1811.03547v1 and https://arxiv.org/pdf/1002.0442 .
All new inequalities are local deductions proved in the report. The target
10,000 remains unproved; no priority or independent-review claim is made.

## Dynamic conditioning and joint prime presence, 2026-10-15

An earlier completed local certificate excludes distinct moduli all at least
14,900, giving minimum modulus at most **14,899**. The proof remains in
`src/conditional-stoploss-improvement.tex`. Exact inputs and completed integer
verification are `data/sparse-certificate-14900.json` and
`data/sparse-verified-14900.json`. Intermediate complete certificates also
establish 15,499 and 14,999. The latest index points to the newer result above.

The additional local deductions are dynamic comparison weights for nonzero
distortions on conditioned primes, a joint Bernoulli presence correction that
preserves dependence between forbidden divisors and product size, and convex
interpolation of the third-moment tail growth exponent. The finite calculation
uses 891,573 conditioned states and 150,000 primes. At prime 2,015,177, its
certified finite survivor exceeds 0.0012690976240386, while the infinite tail
is below 0.0012352021012060. The certified positive margin exceeds
0.0000338955228326. All proof arithmetic uses outward integer bounds.

The search optimizes the terminal criterion jointly after a fixed prefix of
100 primes; it does not claim an optimum over all early parameters or a global
minimum. Separate early-prefix perturbations, earlier distortion starts, finer
bins, and the 10,000 target were tested without obtaining further certificates.
The experimental records are indexed in `data/conditioning-experiment-summary.json`.
The external premises are unchanged from the BBMST and Dusart inputs below.
No independent peer review, proof-assistant formalization, or publication
priority is asserted.

## Full-horizon optimization and third-moment tail, 2026-10-15

An earlier completed certificate excludes distinct moduli all at least 17,000,
giving minimum modulus at most 16,999. Its exact inputs and verification are
`data/sparse-certificate-17000.json` and `data/sparse-verified-17000.json`.
The proof is in the existing `src/conditional-stoploss-improvement.tex`.

The new local deductions are rounded-product graph support bounds, adjoint
optimization of the terminal criterion, a decreasing-event normalization
inequality, and an explicit third-moment bound for the infinite prime tail.
The normalization refinement is proved and tested but is not needed by the
17,000 certificate. The third-moment stopping argument is needed by this
particular certificate; BBMST's older stopping comparison does not pass.

Primary sources checked directly:

- BBMST, https://arxiv.org/html/1811.03547v1 : distortion construction,
  prefix preservation, and the union bound in equation (13), in the proof
  of Theorem 3.1. The earlier 18,000 certificate still uses Theorem 6.1.
- Pierre Dusart, https://arxiv.org/pdf/1002.0442 , Theorem 6.9, equation (6.6):
  pi(x) <= x/(log(x)-1.1) for x >= 60184. This unconditional explicit bound
  is the additional external premise in the new tail proof.

The tail calculation and finite certificate use outward integer intervals.
The 10,000 and 16,900 experiments remain unsuccessful numerical searches,
not impossibility results. No publication-priority or independent-review
claim is made.

## Stronger conditioning certificate for cutoff 20,000, 2026-10-14

The existing conditional-loss report has been extended in place. The new
local deductions are: general conditioning on an initial prime set; support
and chord bounds using true first moments within downward product bins; and
a bounded conditioning correction with the unrestricted excess bound retained
outside the finite conditioning set. Full proofs and all rounding directions
are in `src/conditional-stoploss-improvement.tex`.

BBMST's original arXiv version was reread, especially Theorem 6.1 and Sections
6–8: https://arxiv.org/html/1811.03547v1 . Its distortion construction and
analytic termination theorem are unchanged external inputs. The new finite
loss bounds are proved locally and are not attributed to BBMST.

`data/sparse-certificate-20000-50000.json` passes the standard-library integer
verifier `code/verify_sparse_conditioning.py`, with result and witness hash in
`data/sparse-verified-20000.json`. The exclusion cutoff is 20,000, hence the
integer upper bound is 19,999. At k=50,000 the certified stopping margin
exceeds 673,136.818. Intermediate full runs certify cutoffs 30,000 and 22,000.

Parameter screening, larger conditioning sets, first-moment bin experiments,
dynamic conditioning, and an experimental association correction are recorded
in the new scripts and data files. There is no certificate for 10,000, and
the exploratory 18,000/19,000 endpoints are not presented as proofs.
No claim of publication priority, global optimization, independent peer
review, or proof-assistant formalization is made.

## Conditional convex-loss certificate for cutoff 40,000, 2026-10-14

The new report `src/conditional-stoploss-improvement.tex` proves a product
convex-comparison lemma from the conditional atom caps of BBMST's measures.
The prior weight report in this repository suggested the comparison approach;
the new report includes the complete proof rather than assuming that local
research deduction as an external theorem.

The external premises remain BBMST's distortion construction, marginal
preservation, and Theorem 6.1 in https://arxiv.org/html/1811.03547v1 .
The conditional excess-loss bound, retained 2-and-3 divisors, downward product
bins, and explicit omitted-probability treatment are derived in the report.
The exact rational parameter witness passes the integer verifier at k=50,000,
establishing the exclusion cutoff 40,000 and integer upper bound 39,999.

Additional searches combined covering systems, convex moments, distortion,
minimum modulus, and 616000. They located BBMST's author-hosted Oxford lecture
slides at https://www.memphis.edu/msci/people/pbalistr/oxford_2018.pdf and later
number-field/function-field generalizations, but no external theorem is used
to assert the new cutoff. This is a bounded search, not a novelty guarantee.
The result has not been independently peer reviewed or formally verified.

## Minimum-modulus numerical certificate, 2026-10-14

Primary mathematical source read for this investigation:
https://arxiv.org/html/1811.03547v1 (also https://arxiv.org/pdf/1811.03547).
The relevant premises are Lemmas 2.1, 2.2, 3.3, 3.6 and 3.7; Theorem 6.1;
the recursion and Table 1 in Section 6; and the restricted smooth sums in
Section 8. Section 7 illustrates the same stopping criterion but contributes
no additional hypothesis to the new certificate. Numbering is the 2018 v1
preprint's. The final published paper was not separately audited.

The new proof uses the analytic Theorem 6.1 directly after explicitly processing
1,050,000 primes. It does not trust Table 1's rounded numerical thresholds.
The reported cutoff 615,300 and integer upper bound 615,299 are local deductions
verified with the supplied rational witness and integer interval code.
The first-moment sum and restricted second-moment majorant are BBMST's;
the finite gcd-convolution evaluation is derived in the accompanying report.

Searches of the title, author names, and the 616,000 bound also located the
authorial exposition arXiv:2211.01417, later number-field and function-field
generalizations, and https://www.erdosproblems.com/2 . No priority or exhaustive
literature claim is made. Only BBMST's primary mathematical results are used
as external premises for the improved-bound proof.

## Squarefree construction target 30, 2026-10-15

- Jack Dalton and Ognian Trifonov, *Extreme Covering Systems*, Journal of
  Integer Sequences 25 (2022), Article 22.9.1. Problem 4 asks whether every
  distinct squarefree covering has least modulus 2; the introduction also
  describes the minimum-3 existence question as open. Read from the
  [journal PDF](https://cs.uwaterloo.ca/journals/JIS/VOL25/Trifonov/trif3.pdf).
- Maria Cummings, Michael Filaseta, and Ognian Trifonov, *An upper bound for
  the minimum modulus in a covering system with squarefree moduli*,
  [arXiv:2211.08548](https://arxiv.org/abs/2211.08548), 2022 preprint,
  Theorem 1.1: every distinct squarefree covering has a modulus at most 118.
  This bound does not decide the proposed minimum 30.

These are literature statements, not consequences of the local search.
The minimum-30 reciprocal-capacity screen, finite experiment, and explicit
failure witness are local calculations recorded in the construction note.

