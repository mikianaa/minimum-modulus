"""Full enumeration of placements checks the exact coupled exchange rule."""
from itertools import product
import json
from pathlib import Path
import random

from nested_resource_exchange import best_nested_pair
from structural_cover_search import Diagram,ONE


def run():
    rng = random.Random(50)
    cases = 0
    for domains in [[2,3,2],[1,2,3],[3,3],[2,2,2,2]]:
        points = list(product(*(range(n) for n in domains)))
        for _ in range(100):
            diagram,root = Diagram(domains),ONE
            residual = set(points)
            for _ in range(rng.randrange(9)):
                selected = rng.sample(range(len(domains)),rng.randrange(len(domains)+1))
                a = {v:rng.randrange(domains[v]) for v in selected}
                root = diagram.erase(root,a)
                residual = {x for x in residual if not all(x[v]==r for v,r in a.items())}
            fine = sorted(rng.sample(range(len(domains)),rng.randrange(len(domains)+1)))
            coarse = sorted(rng.sample(fine,rng.randrange(len(fine)+1)))
            score,a,b = best_nested_pair(diagram,root,coarse,fine)
            assert set(a)==set(coarse) and set(b)==set(fine)
            def hits(x,assignment):
                return all(x[v]==r for v,r in assignment.items())
            actual = sum(hits(x,a) or hits(x,b) for x in residual)
            brute = 0
            for av in product(*(range(domains[v]) for v in coarse)):
                aa = dict(zip(coarse,av))
                for bv in product(*(range(domains[v]) for v in fine)):
                    bb = dict(zip(fine,bv))
                    brute = max(brute,sum(hits(x,aa) or hits(x,bb) for x in residual))
            assert score == actual == brute
            diagram.count.cache_clear()
            cases += 1
    result = dict(status='PASS',exhaustive_pair_placement_comparisons=cases,
                  includes_empty_residual=True,includes_unit_domains=True,
                  includes_equal_supports=True)
    Path('data/nested-resource-exchange-tests.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)


if __name__ == '__main__':
    run()
