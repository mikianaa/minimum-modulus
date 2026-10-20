"""Independent modular-arithmetic checks of rejected construction candidates.

For each row, exhibit a coordinate excluding its whole infinite tail family
and all finite closing classes. This does not use the trie's state partition.
"""
from pathlib import Path
import json
from math import prod
import argparse

ROOT = Path(__file__).resolve().parents[1]


def verify_uncovered_progression(data):
    ps, ds, rows = data['primes'], data['depths'], data['rows']
    witness = data['verification']['uncovered_progression']
    R, L = witness['residue'], witness['modulus']
    assert L == prod(p**E for p,E in zip(ps,ds)) and 0 <= R < L
    assert len({tuple(r['exponents']) for r in rows}) == len(rows)
    checks = 0
    for row in rows:
        assert prod(p**e for p,e in zip(ps,row['exponents'])) >= data['minimum']
        excluded = False
        for p,E,e,r in zip(ps,ds,row['exponents'],row['residues']):
            if not e:
                continue
            if e < E:
                excluded = R % p**e != r
            else:
                assert e == E
                a,d0 = divmod(r,p-1)
                earliest = a+(d0+1)*p**(E-1)
                # Later tail classes and all closing classes require the
                # prefix a modulo p^E. The first class requires earliest.
                excluded = R % p**E not in (a,earliest)
            if excluded:
                break
        assert excluded, 'Witness did not exclude a completed row'
        checks += 1
    return {'status':'UNCOVERED_PROGRESSION_VERIFIED', 'rows_checked':checks,
            'minimum_target':data['minimum'], 'uncovered_progression':witness}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files',nargs='+')
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    results = {}
    for filename in args.files:
        data = json.loads((ROOT/filename).read_text())
        results[filename] = verify_uncovered_progression(data)
    (ROOT/args.output).write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps({name:{k:v for k,v in result.items() if k != 'uncovered_progression'}
                      for name,result in results.items()},indent=2))
