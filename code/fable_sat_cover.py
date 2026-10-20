"""Exact SAT search for distinct covering systems with prescribed LCM.

Question answered: given N and M, is there a set of congruence classes
r_d (mod d), one for each d in some subset of the divisors d >= M of N,
covering every integer?  Because all moduli divide N it suffices to cover
the residues 0..N-1.  The encoding is exact, so UNSAT is a proof that no
distinct covering system with min modulus >= M has all moduli dividing N.

Usage:
  python code/fable_sat_cover.py --N 120 --M 3
  python code/fable_sat_cover.py --scan --M 4 --limit 20000 [--primes 2,3,5,7]
Requires python-sat (tmp/fable-venv).  Fable-owned file; see memory fable-loop-plan.
"""
import argparse, json, math, sys, time
from fractions import Fraction
from pysat.solvers import Cadical153, Glucose4
from pysat.card import CardEnc, EncType


def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i * i != n:
                ds.append(n // i)
        i += 1
    return sorted(ds)


def solve(N, M, max_literals=200_000_000, verbose=False, time_limit=None):
    """time_limit: seconds; returns status UNKNOWN when exceeded."""
    D = [d for d in divisors(N) if d >= M]
    if sum(Fraction(1, d) for d in D) < 1:
        return dict(N=N, M=M, status="UNSAT", reason="sum of 1/d below 1", seconds=0.0)
    if N * len(D) > max_literals:
        return dict(N=N, M=M, status="SKIPPED", reason="encoding too large", literals=N * len(D))
    t = time.time()
    var = {}
    nxt = 1
    for d in D:
        for r in range(d):
            var[(d, r)] = nxt
            nxt += 1
    # CaDiCaL 1.5.3 in python-sat does not honour interrupt(); use Glucose4 for timed runs
    s = Glucose4() if time_limit else Cadical153()
    top = nxt - 1
    for d in D:
        lits = [var[(d, r)] for r in range(d)]
        cnf = CardEnc.atmost(lits=lits, bound=1, top_id=top, encoding=EncType.seqcounter)
        top = max(top, cnf.nv)
        for c in cnf.clauses:
            s.add_clause(c)
    for n in range(N):
        s.add_clause([var[(d, n % d)] for d in D])
    # symmetry breaking: translation x -> x - r lets the class with modulus N, if used, have residue 0
    for r in range(1, N):
        s.add_clause([-var[(N, r)]])
    if time_limit:
        import threading
        timer = threading.Timer(time_limit, lambda: s.interrupt())
        timer.start()
        ok = s.solve_limited(expect_interrupt=True)
        timer.cancel()
    else:
        ok = s.solve()
    status = "SAT" if ok else ("UNSAT" if ok is False else "UNKNOWN")
    res = dict(N=N, M=M, status=status, seconds=round(time.time() - t, 2), moduli_available=len(D))
    if ok:
        model = set(l for l in s.get_model() if l > 0)
        system = sorted((d, r) for (d, r), v in var.items() if v in model)
        # independent check
        assert len({d for d, _ in system}) == len(system), "moduli not distinct"
        assert all(any(n % d == r for d, r in system) for n in range(N)), "coverage check failed"
        res["system"] = [dict(modulus=d, residue=r) for d, r in system]
        res["moduli_used"] = len(system)
        res["min_modulus"] = min(d for d, _ in system)
    s.delete()
    return res


def candidates(M, limit, primes):
    """N <= limit, built from the given primes, ordered by size."""
    out = [1]
    for p in primes:
        new = []
        for n in out:
            m = n * p
            while m <= limit:
                new.append(m)
                m *= p
        out += new
    return sorted(n for n in out if n >= M)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int)
    ap.add_argument("--M", type=int, required=True)
    ap.add_argument("--scan", action="store_true", help="smallest feasible N up to --limit")
    ap.add_argument("--limit", type=int, default=10000)
    ap.add_argument("--primes", default="2,3,5,7,11,13")
    ap.add_argument("--output")
    ap.add_argument("--time-limit", type=float, default=None, help="seconds per N; UNKNOWN when exceeded")
    a = ap.parse_args()
    results = []
    if a.scan:
        primes = [int(x) for x in a.primes.split(",")]
        for N in candidates(a.M, a.limit, primes):
            r = solve(N, a.M, time_limit=a.time_limit)
            if r["status"] != "UNSAT" or r.get("reason") != "sum of 1/d below 1":
                print(N, r["status"], r.get("seconds"), r.get("reason", ""), flush=True)
            results.append(r)
            if r["status"] == "SAT":
                print("smallest feasible N found:", N, "system:", [(x["residue"], x["modulus"]) for x in r["system"]])
                break
    else:
        r = solve(a.N, a.M, verbose=True, time_limit=a.time_limit)
        results.append(r)
        print(json.dumps(r, indent=1))
    if a.output:
        json.dump(dict(question="distinct covering system, all moduli | N, min modulus >= M", results=results), open(a.output, "w"), indent=1)


if __name__ == "__main__":
    main()
