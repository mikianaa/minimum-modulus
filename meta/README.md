# Meta-analysis log of the agents' research progress

Records how the two agents (Codex, Claude Code "Fable") advance the determination of the
minimum modulus M, for later analysis of agent research behaviour. Separate from
`PROGRESS.md` (which states the mathematics) and from the dated notes (which give the
details).

## Files

| File | Role |
| --- | --- |
| `agent-progress.jsonl` | Source of truth. One JSON record per milestone or attempt, hand-written at the time. Lines starting with `#` are comments. |
| `meta_stats.py` | Reads the local Claude Code and Codex transcripts, aggregates token usage, and renders the log. Standard library only. |
| `token-usage.json` | Output: tokens, API-list-price cost and wall-clock per session and per day, per agent. |
| `agent-progress.md` | Output: the log as a table with the cumulative tokens spent at each milestone. |

Regenerate the two outputs with `python meta/meta_stats.py` from the repository root.
Transcripts are read from `~/.claude/projects/<project>/*.jsonl` (usage per API call) and
`~/.codex/sessions/**/rollout-*.jsonl` (per-turn `token_count` events, filtered to sessions
whose working directory is this project or its parent `covering-weight-mu`), from 14 Oct 2026.

## Metrics and how each is recorded

| Metric (user's list) | Field(s) in `agent-progress.jsonl` | How it is measured |
| --- | --- | --- |
| Final progress: rigorously proved bounds | `bounds_after.lower_rigorous`, `.lower_conditional`, `.upper_certified` | State of knowledge after the entry. `lower_rigorous` stays at Owens' published 42 until a candidate is verified end to end; `lower_conditional` is the best kernel-verified candidate (depends on Owens' thesis and its transcription); `upper_certified` is the best cutoff passed by the integer verifier. |
| Compute resources: tokens, cost, time | `compute.tokens` (computed), `compute.cost_usd` (computed), `compute.wall_clock_hours` (hand-entered, null if unknown) | Tokens and cost come from the transcripts via `meta_stats.py`; cost is API list price (Fable 5.1: 10 / 50 USD per MTok in / out, cache read 0.25; Codex rate unset), so it is an estimate, not a bill. |
| Path to discovery: at which token the improvement appeared | computed column "tokens at discovery" in `agent-progress.md` | Cumulative tokens of that agent (and of both agents) on this project up to the entry's `timestamp`. |
| Search efficiency: candidates and success rate | `search.candidates`, `search.successes` | Hand-entered: number of distinct candidates / designs / runs tried in the attempt, and how many led to a verified result. Null when the attempt had no enumerable candidates. |
| Mathematical method | `method` (list of tags) | `parameter-optimization`, `inequality-improvement`, `new-lemma`, `conditioning`, `construction-patch`, `exact-search`, `density-bound`, `audit`, `numerical-experiment`, `literature`, `tooling`. |
| Cause of failure | `failure_cause` | `numerical` (a numerical estimate was wrong or a numerical refutation), `logical-gap` (no proof route for an idea), `search-dead-end` (solver / search limits), `capacity-shortfall` (a counting or density obstruction), `none`. |

Other fields: `id`, `timestamp` (ISO 8601 with offset; `approx: true` when reconstructed from
file dates or memory rather than logged), `agent` (`codex` or `claude-code`), `model`,
`side` (`upper`, `lower`, `side`), `kind` (`certificate`, `construction`, `obstruction`,
`lead`, `refutation`, `audit`, `tool`), `result.summary`, `result.rigor`
(`certified`, `conditional`, `candidate`, `none`), `notes`, `files`.

## Adding an entry

Append one line to `agent-progress.jsonl` when a milestone is reached or an attempt is
closed, then run `python meta/meta_stats.py`. Keep entries factual; the computed columns
should never be written by hand.
