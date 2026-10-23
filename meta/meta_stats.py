"""Aggregate agent token usage and annotate the progress log for meta-analysis.

Reads
  * Claude Code transcripts  (~/.claude/projects/<project>/*.jsonl)
  * Codex rollouts           (~/.codex/sessions/YYYY/MM/DD/rollout-*.jsonl)
  * meta/agent-progress.jsonl (one record per milestone / attempt, hand-written)

Writes
  * meta/token-usage.json   per-session and per-day token usage, totals, cost
  * meta/agent-progress.md  the progress log rendered as a table, with the
                            cumulative tokens spent (per agent and overall) at
                            each milestone's timestamp

Only the standard library is used.  Run from the repository root:

    python meta/meta_stats.py

Pricing is API-list pricing and is only an estimate of cost: the Claude Code
sessions ran on a subscription, and Codex rates are not recorded locally
(set CODEX_RATES below if known).
"""
from __future__ import annotations

import datetime as dt
import glob
import json
import os
import sys
from collections import defaultdict

HOME = os.path.expanduser("~")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
META = os.path.join(ROOT, "meta")

# Working directories whose sessions belong to this project (the project was
# split out of covering-weight-mu on 20 Oct 2026; sessions in that directory
# from START_DATE on are minimum-modulus work).
CWD_KEYS = ("covering-weight-mu", "minimum-modulus")
CLAUDE_PROJECT_DIRS = [
    os.path.join(HOME, ".claude", "projects", "c--Users-mekak-code-covering-weight-mu"),
    os.path.join(HOME, ".claude", "projects", "c--Users-mekak-code-minimum-modulus"),
]
CODEX_SESSIONS = os.path.join(HOME, ".codex", "sessions")
START_DATE = dt.datetime(2026, 10, 14, tzinfo=dt.timezone.utc)

# USD per million tokens.  Claude Fable 5.1 list price (claude-api skill, cached
# 2026-10-06): input 10, output 50, cache read 0.25; cache writes 1.25x (5 min)
# and 2x (1 h) of input.
CLAUDE_RATES = {
    "claude-fable-5-1": {"input": 10.0, "output": 50.0, "cache_read": 0.25,
                         "cache_write_5m": 12.5, "cache_write_1h": 20.0},
    "claude-sonnet-5-5": {"input": 2.0, "output": 10.0, "cache_read": 0.20,
                          "cache_write_5m": 2.5, "cache_write_1h": 4.0},
}
CLAUDE_DEFAULT_MODEL = "claude-fable-5-1"
CODEX_RATES = None  # e.g. {"input": ..., "cached": ..., "output": ...}


def parse_ts(s: str) -> dt.datetime:
    return dt.datetime.fromisoformat(s.replace("Z", "+00:00"))


def claude_events():
    """One event per API call: (ts, session, model, counts)."""
    seen = set()
    for d in CLAUDE_PROJECT_DIRS:
        # top-level session files plus subagent transcripts under <session>/subagents/
        for f in sorted(glob.glob(os.path.join(d, "**", "*.jsonl"), recursive=True)):
            session = os.path.relpath(f, d).split(os.sep)[0][:8]
            with open(f, encoding="utf-8") as fh:
                for line in fh:
                    if '"assistant"' not in line:
                        continue
                    try:
                        x = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    if x.get("type") != "assistant":
                        continue
                    msg = x.get("message") or {}
                    u = msg.get("usage")
                    mid = msg.get("id")
                    if not u or not mid or mid in seen:
                        continue
                    seen.add(mid)
                    cc = u.get("cache_creation") or {}
                    yield {
                        "ts": parse_ts(x["timestamp"]),
                        "agent": "claude-code",
                        "session": session,
                        "model": msg.get("model"),
                        "input": u.get("input_tokens", 0),
                        "cache_write_5m": cc.get("ephemeral_5m_input_tokens", 0),
                        "cache_write_1h": cc.get("ephemeral_1h_input_tokens", 0),
                        "cache_read": u.get("cache_read_input_tokens", 0),
                        "output": u.get("output_tokens", 0),
                    }


def codex_events():
    """One event per token_count message (the per-turn delta)."""
    for f in sorted(glob.glob(os.path.join(CODEX_SESSIONS, "*", "*", "*", "rollout-*.jsonl"))):
        session = os.path.basename(f)[8:27]
        cwd = None
        model = None
        events = []
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                try:
                    x = json.loads(line)
                except json.JSONDecodeError:
                    continue
                p = x.get("payload")
                if not isinstance(p, dict):
                    continue
                cwd = cwd or p.get("cwd")
                model = model or p.get("model")
                if p.get("type") == "token_count" and p.get("info"):
                    last = p["info"]["last_token_usage"]
                    events.append({
                        "ts": parse_ts(x["timestamp"]),
                        "agent": "codex",
                        "session": session,
                        "model": None,
                        "input": last.get("input_tokens", 0) - last.get("cached_input_tokens", 0),
                        "cache_write_5m": 0,
                        "cache_write_1h": 0,
                        "cache_read": last.get("cached_input_tokens", 0),
                        "output": last.get("output_tokens", 0),
                    })
        if not cwd or not any(k in cwd for k in CWD_KEYS):
            continue
        for e in events:
            e["model"] = model
            yield e


def total(e):
    return e["input"] + e["cache_write_5m"] + e["cache_write_1h"] + e["cache_read"] + e["output"]


def cost(e):
    if e["agent"] == "claude-code":
        r = CLAUDE_RATES.get(e["model"] or "", CLAUDE_RATES[CLAUDE_DEFAULT_MODEL])
        return (e["input"] * r["input"] + e["output"] * r["output"] + e["cache_read"] * r["cache_read"]
                + e["cache_write_5m"] * r["cache_write_5m"] + e["cache_write_1h"] * r["cache_write_1h"]) / 1e6
    if CODEX_RATES:
        r = CODEX_RATES
        return (e["input"] * r["input"] + e["cache_read"] * r["cached"] + e["output"] * r["output"]) / 1e6
    return None


def summarise(events):
    per_session = {}
    per_day = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    per_model = defaultdict(lambda: {"total": 0, "output": 0, "calls": 0})
    for e in events:
        k = (e["agent"], e["session"])
        s = per_session.setdefault(k, {"agent": e["agent"], "session": e["session"], "models": defaultdict(int),
                                       "start": e["ts"], "end": e["ts"], "calls": 0,
                                       "input": 0, "cache_write": 0, "cache_read": 0, "output": 0,
                                       "total": 0, "cost_usd": 0.0 if cost(e) is not None else None})
        s["models"][e["model"] or "?"] += total(e)
        pm = per_model[f"{e['agent']}/{e['model'] or '?'}"]
        pm["total"] += total(e)
        pm["output"] += e["output"]
        pm["calls"] += 1
        s["start"] = min(s["start"], e["ts"])
        s["end"] = max(s["end"], e["ts"])
        s["calls"] += 1
        s["input"] += e["input"]
        s["cache_write"] += e["cache_write_5m"] + e["cache_write_1h"]
        s["cache_read"] += e["cache_read"]
        s["output"] += e["output"]
        s["total"] += total(e)
        c = cost(e)
        if c is not None:
            s["cost_usd"] += c
        day = e["ts"].astimezone(dt.timezone(dt.timedelta(hours=9))).date().isoformat()
        per_day[day][e["agent"]]["total"] += total(e)
        per_day[day][e["agent"]]["output"] += e["output"]
    for s in per_session.values():
        s["hours"] = round((s["end"] - s["start"]).total_seconds() / 3600, 2)
        s["models"] = dict(sorted(s["models"].items(), key=lambda kv: -kv[1]))
        s["start"] = s["start"].isoformat()
        s["end"] = s["end"].isoformat()
        if s["cost_usd"] is not None:
            s["cost_usd"] = round(s["cost_usd"], 2)
    totals = defaultdict(lambda: {"total": 0, "output": 0, "cost_usd": 0.0, "sessions": 0})
    for s in per_session.values():
        t = totals[s["agent"]]
        t["total"] += s["total"]
        t["output"] += s["output"]
        t["sessions"] += 1
        if s["cost_usd"] is not None:
            t["cost_usd"] += s["cost_usd"]
        else:
            t["cost_usd"] = None
    return {
        "generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "start_date": START_DATE.date().isoformat(),
        "rates_usd_per_mtok": {"claude-code": CLAUDE_RATES, "codex": CODEX_RATES},
        "totals": {k: dict(v) for k, v in totals.items()},
        "per_model": dict(per_model),
        "per_day": {d: {a: dict(v) for a, v in agents.items()} for d, agents in sorted(per_day.items())},
        "per_session": sorted(per_session.values(), key=lambda s: s["start"]),
    }


def cumulative_at(events, ts):
    """Tokens spent by each agent (and in total) before timestamp ts."""
    out = defaultdict(int)
    for e in events:
        if e["ts"] <= ts:
            out[e["agent"]] += total(e)
            out["all"] += total(e)
    return out


def load_progress():
    path = os.path.join(META, "agent-progress.jsonl")
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                rows.append(json.loads(line))
    return rows


def fmt_m(n):
    return f"{n/1e6:.1f}M" if n else "-"


def render_progress(rows, events):
    lines = ["# Agent progress (rendered by meta/meta_stats.py; source: agent-progress.jsonl)", "",
             "Cumulative tokens are all tokens (input, cache, output) spent by the agent on this",
             "project from 14 Oct 2026 up to the milestone's timestamp. `rigor`: certified = integer",
             "verifier passed; conditional = depends on an unverified source; candidate = not yet",
             "checked; none = no new bound.", "",
             "| # | time (JST) | agent | side | kind | result | rigor | method | candidates / successes | failure cause | tokens at discovery (agent / all) | wall-clock |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(rows, 1):
        ts = parse_ts(r["timestamp"])
        cum = cumulative_at(events, ts)
        agent = r["agent"]
        res = r.get("result", {})
        search = r.get("search", {})
        cand = search.get("candidates")
        succ = search.get("successes")
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} / {} | {} |".format(
            i, ts.astimezone(dt.timezone(dt.timedelta(hours=9))).strftime("%m-%d %H:%M"),
            agent, r.get("side", ""), r.get("kind", ""),
            res.get("summary", ""), res.get("rigor", ""),
            ", ".join(r.get("method", [])),
            f"{cand} / {succ}" if cand is not None else "-",
            r.get("failure_cause") or "-",
            fmt_m(cum[agent]), fmt_m(cum["all"]),
            f"{r['compute']['wall_clock_hours']} h" if r.get("compute", {}).get("wall_clock_hours") is not None else "-",
        ))
    lines.append("")
    lines.append("Notes per entry:")
    lines.append("")
    for i, r in enumerate(rows, 1):
        if r.get("notes"):
            lines.append(f"{i}. {r['notes']}" + (f" Files: {', '.join(r['files'])}." if r.get("files") else ""))
    return "\n".join(lines) + "\n"


def main():
    events = [e for e in list(claude_events()) + list(codex_events()) if e["ts"] >= START_DATE]
    events.sort(key=lambda e: e["ts"])
    usage = summarise(events)
    with open(os.path.join(META, "token-usage.json"), "w", encoding="utf-8") as fh:
        json.dump(usage, fh, indent=1)
    rows = load_progress()
    with open(os.path.join(META, "agent-progress.md"), "w", encoding="utf-8") as fh:
        fh.write(render_progress(rows, events))
    for agent, t in usage["totals"].items():
        c = f"${t['cost_usd']:.0f}" if t["cost_usd"] is not None else "n/a"
        print(f"{agent:12s} sessions {t['sessions']:3d}  tokens {t['total']/1e6:8.1f}M  output {t['output']/1e6:6.2f}M  cost {c}")
    print(f"{len(rows)} progress entries rendered")


if __name__ == "__main__":
    sys.exit(main())
