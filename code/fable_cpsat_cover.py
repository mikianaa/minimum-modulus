"""CP-SAT (OR-Tools) search for distinct covering systems with moduli dividing N
and minimum modulus >= M.  Same exact model as code/fable_sat_cover.py:
one Boolean per (modulus d | N, d >= M, residue r), at most one residue per
modulus, every residue class mod N covered.  Symmetry: modulus N, if used,
has residue 0 (translation).  UNSAT (status INFEASIBLE) is a proof of nonexistence.

Usage: python code/fable_cpsat_cover.py --M 4 --N 360 [--time 600] [--workers 8]
Fable-owned file; see memory fable-loop-plan.
"""
import argparse, json, math, time
from fractions import Fraction
from ortools.sat.python import cp_model


def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i * i != n: ds.append(n // i)
        i += 1
    return sorted(ds)


def solve(N, M, time_limit=600.0, workers=8, require=(), forbid=(), redundant=False, linearization=1):
    D = [d for d in divisors(N) if d >= M and d not in forbid]
    if sum(Fraction(1, d) for d in D) < 1:
        return dict(N=N, M=M, status="INFEASIBLE", reason="sum of 1/d below 1", seconds=0.0)
    t = time.time()
    m = cp_model.CpModel()
    x = {(d, r): m.NewBoolVar(f"x{d}_{r}") for d in D for r in range(d)}
    for d in D:
        m.AddAtMostOne(x[(d, r)] for r in range(d))
    for n in range(N):
        m.AddBoolOr([x[(d, n % d)] for d in D])
    for r in range(1, N):
        m.Add(x[(N, r)] == 0)
    for d in require:
        m.AddBoolOr([x[(d, r)] for r in range(d)])
    if redundant:
        # implied counting constraints: for g=N/p every class c mod g needs total weight >= 1,
        # where modulus d contributes gcd(d,g)/d if its residue is = c mod gcd(d,g).
        primes = [q for q in range(2, 200) if N % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
        for q in primes:
            g = N // q
            agg = {}
            for d in D:
                h = math.gcd(d, g)
                for s_ in range(h):
                    agg[(d, s_)] = sum(x[(d, r)] for r in range(s_, d, h))  # residue r = s_ mod h
            for c in range(g):
                m.Add(sum((N * math.gcd(d, g) // d) * agg[(d, c % math.gcd(d, g))] for d in D) >= N)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_workers = workers
    s.parameters.linearization_level = linearization
    st = s.Solve(m)
    name = s.StatusName(st)
    res = dict(N=N, M=M, status=name, seconds=round(time.time() - t, 1), moduli_available=len(D))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        system = sorted((d, r) for (d, r), v in x.items() if s.Value(v))
        assert len({d for d, _ in system}) == len(system)
        assert all(any(n % d == r for d, r in system) for n in range(N)), "coverage check failed"
        res["status"] = "SAT"
        res["system"] = [dict(modulus=d, residue=r) for d, r in system]
        res["moduli_used"] = len(system); res["min_modulus"] = min(d for d, _ in system)
    elif st == cp_model.INFEASIBLE:
        res["status"] = "UNSAT"
    else:
        res["status"] = "UNKNOWN"
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--M", type=int, required=True); ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--time", type=float, default=600.0); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--output"); ap.add_argument("--redundant", action="store_true"); ap.add_argument("--lin", type=int, default=1)
    a = ap.parse_args()
    r = solve(a.N, a.M, a.time, a.workers, redundant=a.redundant, linearization=a.lin)
    print(json.dumps({k: v for k, v in r.items() if k != "system"}))
    if "system" in r: print("system:", [(e["residue"], e["modulus"]) for e in r["system"]])
    if a.output: json.dump(r, open(a.output, "w"), indent=1)
