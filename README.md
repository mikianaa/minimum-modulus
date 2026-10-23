# Minimum modulus of distinct covering systems

How large can the smallest modulus of a finite covering system with distinct moduli be?
Write M for the supremum. Known before this project: M >= 42 (Owens 2014, explicit
construction) and M <= 616,000 (Balister, Bollobas, Morris, Sahasrabudhe, Tiba, "BBMST").
This repository attacks the question from both sides with computer-assisted proofs and
explicit constructions. Progress and status live in [PROGRESS.md](PROGRESS.md); this file
only says what the project is and how it is laid out.

## What is being done

* **Upper side.** BBMST's distortion method is reimplemented with a convex-comparison
  stage bound, exact rational witnesses and an integer-arithmetic verifier. Each certified
  cutoff K proves "every finite distinct covering system has a modulus < K". The current
  certificate is K = 14,485, i.e. M <= 14,484. The working target is 10,000.
* **Lower side.** Owens' construction is transcribed as a finite kernel (prime-power tail
  model) and patched: classes with modulus below the target are deleted and the holes are
  refilled with spare smooth moduli and arrows on fresh primes. Candidate systems with
  minimum modulus 43 and 44 exist; 45 is shown to be unreachable inside this framework.
* **Side results.** Minimal-LCM questions for small minimum modulus (SAT and CP-SAT tools),
  a density bound for the {2,3,5}-smooth problem of HKLT, squarefree variants.

## Layout

| Path | Contents |
| --- | --- |
| `PROGRESS.md` | Current state, timeline, closed routes, open directions. |
| `src/*.tex` | Proofs: `conditional-stoploss-improvement.tex` (upper-bound certificates), `minimum-modulus-improvement.tex` (first 615,299 certificate), `construction-search.tex` (lower-side search models and obstructions). |
| `code/` | Searches, verifiers and tests. Files `fable_*.py` are Fable's tools; the rest are Codex's. Verifiers use only the standard library; searches need NumPy or OR-Tools. |
| `data/` | Exact witnesses, verification reports and experiment results (JSON). Fable's results are `fable-*.json`. |
| `notes/` | Dated working notes: `fable-notes-*.md` (Fable), `codex-readme-log-2026-10-14-15.md` (Codex's detailed log, formerly the body of this README), fraction-bound tables. |
| `literature/` | `SOURCES.md` with references; PDFs are kept locally and ignored by git. |
| `tmp/` | Ignored. Python environments `tmp/fable-venv` (numpy, scipy, ortools, python-sat) and `tmp/construction-deps`. |

## Key artefacts

* Upper bound: `data/latest-certified-bound.json` indexes the strongest certificate
  (witness `data/sparse-certificate-14485.json`, report `data/sparse-verified-14485.json`).
  Replay (about 70 minutes, standard library only):

      python code/check_conditioning_methods.py
      python code/verify_sparse_conditioning.py data/sparse-certificate-14485.json --output data/sparse-verified-14485.json

* Lower bound: `data/fable-m43-construction.json` and `data/fable-m44-construction.json`
  describe the candidate systems (Owens' seed, deleted classes, spare classes, fresh-prime
  pieces) and record how coverage, distinctness and set-building orders were checked.
  Both are conditional on the correctness of Owens' thesis and on the transcription of its
  Sections 3.1-3.4 in `code/owens_seed_core.py`; Sections 3.5-3.20 are not machine-verified.
* Kernel verifier for constructions: `code/construction_kernel.py` (tests in
  `code/test_construction_kernel.py`).

## AI disclosure

This repository is produced with two AI coding agents working in the same tree under
the author's direction: OpenAI Codex and Claude Code (Anthropic; model Claude Fable 5.1,
referred to as "Fable" in file names and notes). Codex owns `src/`, the non-`fable_` files
in `code/` and `data/`, and `literature/`. Fable owns `code/fable_*.py`, `data/fable-*.json`
and `notes/fable-*.md`. Proofs, searches, certificates, verifiers and notes were written
by these agents; the author sets the goals and reviews the results.

**Work in progress.** The results here are intermediate. Independent verification of the
claims (replay of the certificates, audit of the constructions and of the transcription of
Owens' thesis) is being carried out as a separate piece of work and is not finished. Until
it is, every bound in this repository should be read as a computer-assisted claim with the
stated dependencies, not as an established result. Nothing is peer reviewed.

## Conventions

Record results in `PROGRESS.md` (one line per milestone) and in a dated note, not here.
