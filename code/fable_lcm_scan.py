"""Scan for the minimal LCM of a distinct covering system with minimum modulus M.

For each candidate N (smooth over --primes, in [--lo, --hi]) in increasing order:
  1. reciprocal filter: sum_{d|N, d>=M} 1/d >= 1, else UNSAT;
  2. coarse filters: for g = N/p (p prime | N) solve the exact relaxation at
     modulus g: each modulus d chooses one residue s mod gcd(d,g) and covers
     weight gcd(d,g)/d of every class c = s (mod gcd(d,g)) of Z_g; every class
     needs weight >= 1.  Infeasible => no covering system with LCM N (proof);
  3. full CP-SAT model (code/fable_cpsat_cover.py) with a time limit.
Results are appended as JSON lines to --log.  Stops at the first SAT.
Usage: python code/fable_lcm_scan.py --M 8 --lo 10080 --hi 172800 --primes 2,3,5,7 --time 600 --log out.jsonl
Fable-owned file; see memory fable-loop-plan.
"""
import argparse, json, math, sys, time
from fractions import Fraction
from ortools.sat.python import cp_model
sys.path.insert(0, __import__("os").path.dirname(__file__))
from fable_cpsat_cover import divisors, solve as full_solve


def smooth(primes, lo, hi):
    out = [1]
    for p in primes:
        new = []
        for n in out:
            m = n * p
            while m <= hi:
                new.append(m); m *= p
        out += new
    return sorted(n for n in out if lo <= n <= hi)


def coarse(N, M, g, time_limit=60.0, workers=8):
    """Exact relaxation at modulus g. Returns 'INFEASIBLE', 'FEASIBLE' or 'UNKNOWN'."""
    D = [d for d in divisors(N) if d >= M]
    m = cp_model.CpModel()
    y = {}
    for d in D:
        h = math.gcd(d, g)
        for s in range(h):
            y[(d, s)] = m.NewBoolVar(f"y{d}_{s}")
        m.AddAtMostOne(y[(d, s)] for s in range(h))
    scale = N  # weight gcd(d,g)/d * N is an integer
    for c in range(g):
        terms = []
        for d in D:
            h = math.gcd(d, g)
            terms.append((scale * h // d, y[(d, c % h)]))
        m.Add(sum(w * v for w, v in terms) >= scale)
    if N % g == 0 and N in D:
        for s in range(1, math.gcd(N, g)): m.Add(y[(N, s)] == 0)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = time_limit; s.parameters.num_workers = workers
    st = s.Solve(m)
    return {cp_model.INFEASIBLE: "INFEASIBLE", cp_model.OPTIMAL: "FEASIBLE", cp_model.FEASIBLE: "FEASIBLE"}.get(st, "UNKNOWN")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--M", type=int, required=True); ap.add_argument("--lo", type=int, default=1); ap.add_argument("--hi", type=int, required=True)
    ap.add_argument("--primes", default="2,3,5,7"); ap.add_argument("--time", type=float, default=600.0)
    ap.add_argument("--coarse-time", type=float, default=60.0); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--log", required=True); ap.add_argument("--skip-full", action="store_true")
    a = ap.parse_args()
    primes = [int(x) for x in a.primes.split(",")]
    for N in smooth(primes, a.lo, a.hi):
        t = time.time(); rec = dict(N=N, M=a.M)
        D = [d for d in divisors(N) if d >= a.M]
        if sum(Fraction(1, d) for d in D) < 1:
            rec.update(status="UNSAT", how="reciprocal"); 
        else:
            rec["coarse"] = {}
            proven = False
            for p in sorted(set(q for q in primes if N % q == 0), reverse=True):
                g = N // p
                r = coarse(N, a.M, g, a.coarse_time, a.workers); rec["coarse"][str(g)] = r
                if r == "INFEASIBLE":
                    rec.update(status="UNSAT", how=f"coarse g={g}"); proven = True; break
            if not proven:
                if a.skip_full:
                    rec.update(status="UNKNOWN", how="coarse only")
                else:
                    f = full_solve(N, a.M, a.time, a.workers)
                    rec.update(status=f["status"], how="full", full_seconds=f["seconds"])
                    if "system" in f: rec["system"] = f["system"]
        rec["seconds"] = round(time.time() - t, 1)
        print(json.dumps({k: v for k, v in rec.items() if k != "system"}), flush=True)
        with open(a.log, "a") as fh: fh.write(json.dumps(rec) + "\n")
        if rec["status"] == "SAT":
            print("FIRST SAT at N =", N); break


if __name__ == "__main__":
    main()
