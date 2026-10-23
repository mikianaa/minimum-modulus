# Progress

Status of the determination of M, the largest possible minimum modulus of a finite
distinct covering system. Published state of the art: 42 <= M <= 616,000.
Last updated 23 October 2026. Work in progress: independent verification of the results
below is a separate, unfinished task (see the AI disclosure in `README.md`).

## Current state

| Side | Result | Status | Where |
| --- | --- | --- | --- |
| Upper | M <= 14,484 | Certified (integer-arithmetic verifier), audited independently 17 Oct | `data/latest-certified-bound.json`, `src/conditional-stoploss-improvement.tex` |
| Upper target | M <= 9,999 | Not reached. Every single-stage refinement of the comparison lemma is closed; a new idea is needed at the early primes | `notes/fable-notes-2026-10-21-upper.md` |
| Lower | M >= 44 (candidate) | Covering verified on the finite kernel and by integer sampling; distinctness argued structurally; conditional on Owens' thesis and the transcription of its Sections 3.1-3.4 | `data/fable-m44-construction.json` (43: `data/fable-m43-construction.json`) |
| Lower target | 45 | Proved out of reach of the Nielsen/Owens framework (patches, spares, fresh primes, re-partitions, density budget) | `notes/fable-notes-2026-10-20-m45.md`, `notes/fable-m45-fraction-bounds-*.txt` |

Uncommitted work as of 23 Oct: the 21 Oct adversary-toy scripts and results, the M=45
set-count pipeline results, the 21 Oct upper-side note and the 21 Oct addition to the M=45
note.

## Timeline

| Date (2026) | Side | Result | Files |
| --- | --- | --- | --- |
| 14 Oct | upper | First certificate 615,299 (BBMST method, rational witness, integer verifier). Then 39,999 (convex comparison, conditioning at 2 and 3) and 19,999 (conditioning on 64,315 configurations through 53). | `src/minimum-modulus-improvement.tex`, `data/stoploss-verified-40000.json`, `data/sparse-verified-20000.json` |
| 15 Oct | upper | 16,999, 14,599, 14,496, then 14,484 (standalone cutoff 14,485, 300,000 stages). | `data/sparse-verified-17000.json`, `-14600`, `data/prefix-refinement-verified-14497.json`, `data/sparse-verified-14485.json` |
| 15 Oct | lower | Codex's searches: kernel/CP-SAT at 43, structural targets 1000, 100, 60, 50 (relative arrows, joint root repair, modulus-90 relocation), squarefree 30 and 3. No improvement over 42; several restricted obstructions. | `src/construction-search.tex`, `notes/codex-readme-log-2026-10-14-15.md` |
| 17 Oct | upper | Independent audit of 14,484: mathematics sound, replay agrees, method saturated at these parameters (loss budget 0.40 below p=100, 0.45 for 100<p<1000). | `data/conditioning-method-checks.json` |
| 18 Oct | upper | Per-stage loss profile; cap-inflation ceiling 0.126 found, then modulated caps refuted by Monte Carlo (realised gain 1-2 %). | `data/fable-loss-profile-14485.json`, `data/fable-nu-cap-diagnostic-14485.json`, `code/fable_modulated_cap_mc.py` |
| 18 Oct | side | Minimal-LCM frontier collected (m=3..7: 120, 360, 1440, 5040, 10080). SAT and CP-SAT cover tools built; they reproduce m<=4 and fail to find the known m=7 system in 1 h. | `code/fable_sat_cover.py`, `code/fable_cpsat_cover.py`, `code/fable_lcm_scan.py`, `notes/fable-notes-2026-10-18.md` |
| 19 Oct | lower | Candidate with minimum modulus 43: Owens minus 10 mod 42, 14 smooth spares, 44 fresh-prime arrows (89, 97..331). Kernel coverage 0 remaining; 40,000-sample integer check; 44/44 orders certified. | `data/fable-m43-construction.json` |
| 19 Oct | side | HKLT Problem 1 ({2,3,5}-smooth): the 3-stage distortion bound reproduces m<=9 exactly, does not settle m=9. | `code/fable_smooth235_bound.py` |
| 20 Oct | lower | Candidate with minimum modulus 44: additionally deletes the modulus-43 class; 43 fresh arrows (89, 97..317). Kernel-verified. | `data/fable-m44-construction.json` |
| 20 Oct | lower | 45 by patching is infeasible: hole after spares 0.0037 + secondary 0.001 of density vs fresh-prime capacity 0.00085 (primes 331..6000). | `notes/fable-notes-2026-10-20-m45.md` |
| 20 Oct | upper | Truncation alpha<=1 measured at 21-25 % of the stage loss for the tower adversary; class collisions 2-14 %; branching adversaries recover nothing robustly. | `code/fable_class_collision_*.py`, `code/fable_self_saturated.py`, `notes/fable-notes-2026-10-20-upper.md` |
| 21 Oct | lower | 45 impossible inside the Nielsen/Owens framework: rigorous fraction upper bounds per prime and region; odd-branch density sum 0.319 < 0.366 needed (Owens at 42 is exactly at budget). Set-count pipeline confirms (11 on W x 1 mod 4: 7 of 10 sets at 45). | `notes/fable-m45-fraction-bounds-M45.txt`, `-M42.txt`, `data/fable-m45-setcount-pipeline.json` |
| 21 Oct | upper | Truncation lead closed: exact small-instance adversaries reach 93-96 % of the untruncated comparison, so the comparison is within 5-10 % of the true single-stage adversary; 10,000 needs 11 %. Block tail would give only 14,482. | `code/fable_adversary_toy*.py`, `data/fable-adversary-toy-results.json`, `notes/fable-notes-2026-10-21-upper.md` |

## Closed routes (do not retry)

Upper side (all measured or refuted; details in the notes cited above):

* Parameter and cap tuning of the 14,484 certificate; the optimiser has equalised marginal losses.
* Fibre-dependent (modulated) caps: ceiling 9-28 %, realised 1-2 %.
* Lineage double counting / nested-model accounting and per-stage truncation inside the lemma: not valid (Owens' system has sieve mass exactly 1; two-event toy).
* Truncation alpha<=1, class collisions, measure shaping under the caps: exact adversaries realise 93-96 % of the comparison.
* Block tail refinement: 14,484 -> about 14,482.

Lower side:

* Patching the 44-system toward 45 with level>=2 spares, 169-level spares, fresh primes (as arrows or plain classes), 143/187-level sub-arrows, moving 5/7 frame resources.
* Re-partitions: 11 on W x 1 mod 8; 13 on three classes mod 8; 11 on (1 mod 3) x 1 mod 4; 11 on a W half via a 13^ copy.
* Codex's direct kernel/CP-SAT and structural searches at 43, 50, 60, 100, 1000 (restricted obstructions recorded in `src/construction-search.tex`).
* Brute SAT/CP-SAT for minimal LCM beyond N of order 1e5.

## Open directions

* Upper: a two-stage toy (consecutive stages q < p with nu_q > 1) to price the coupling that single-stage comparisons miss; a different distortion rule; arguments that use coverage itself (global residue consistency) rather than moments. The irreducible first-moment part below p = 30 is 0.12.
* Lower: a new odd-branch early design (Owens' hints: prime-3 sets, placement of 81^); a computer search at the set-count level over region assignments; independent verification of Owens' Sections 3.5-3.20, which all candidates depend on.
* Side: minimal LCM for m = 8 (HKLT give 172,800); HKLT's open m = 9 for {2,3,5}-smooth LCM.

## Where the details are

* Upper-side proofs: `src/conditional-stoploss-improvement.tex`, `src/minimum-modulus-improvement.tex`.
* Lower-side models and obstructions: `src/construction-search.tex`.
* Codex's dated log (14-15 Oct, formerly the README body): `notes/codex-readme-log-2026-10-14-15.md`.
* Fable's dated notes: `notes/fable-notes-2026-10-18.md`, `-2026-10-20-upper.md`, `-2026-10-20-m45.md`, `-2026-10-21-upper.md`.
