"""Independent finite enumeration for exact residual and marginal operations."""
from itertools import product
from math import prod
import random
from structural_cover_search import Diagram,ONE,coordinates,condition_from_row,row_from_condition


def main():
    rng = random.Random(100)
    cases = 0
    for domains in [[2,3,2],[1,2,4],[2,2,3,2]]:
        for _ in range(25):
            diagram = Diagram(domains)
            root = ONE
            points = set(product(*(range(d) for d in domains)))
            for _ in range(12):
                selected = sorted(rng.sample(range(len(domains)),rng.randrange(len(domains)+1)))
                actual,assignment,_ = diagram.maximum_cofactor(root,selected)
                scores = []
                for values in product(*(range(domains[v]) for v in selected)):
                    scores.append(sum(all(x[v] == a for v,a in zip(selected,values)) for x in points))
                assert actual == max(scores)
                wg,wa,_ = diagram.whole_cofactor(root,selected)
                volume = prod(domains)//prod(domains[v] for v in selected)
                assert wg == (volume if max(scores) == volume else 0)
                if wg:
                    assert wg == sum(all(x[v] == a for v,a in wa.items()) for x in points)
                dg,da,_ = diagram.digitwise_cofactor(root,selected)
                assert dg == sum(all(x[v] == a for v,a in da.items()) for x in points)
                assert dg <= actual
                hit = {x for x in points if all(x[v] == a for v,a in assignment.items())}
                assert len(hit) == actual
                root = diagram.erase(root,assignment)
                points -= hit
                assert diagram.count(root) == len(points)
                compact,new_root = diagram.compact(root)
                assert compact.count(new_root) == len(points)
                cases += 1
    for ps,ds in [([2,3,5],[3,2,1]),([3,7],[2,2])]:
        for layout in ['prime-first','large-prime-first','digit-level']:
            domains,blocks = coordinates(ps,ds,layout)
            for es in product(*(range(E+1) for E in ds)):
                chosen = {v:rng.randrange(domains[v]) for block,e in zip(blocks,es) for v in block[:e]}
                row = row_from_condition(es,chosen,ps,ds,blocks)
                assert condition_from_row(row,ps,ds,blocks) == chosen
    print({'status':'PASS','exact_marginal_and_residual_comparisons':cases})


if __name__ == '__main__':
    main()
