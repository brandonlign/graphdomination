#!/usr/bin/env python3
"""Count the guide's full-view multiset independently of delivered generators.

Standard library only. Every stub is a distinct (neighbor, instance) pair.
First-appearance block order removes block-order permutations, but equal
incidence sets in distinct blocks are NOT merged. Different partitions of
distinct stub instances can produce identical recorded incidence lists;
these remain separate occurrences, as required for this multiset count.
This count corroborates the grammar; it is not a proof of graph coverage
or a replay of features, polynomial signs, or the large certificate.
"""
from itertools import combinations,product
from pathlib import Path
import hashlib,json

HERE=Path(__file__).resolve().parent

def require(ok,message):
    if not ok:raise RuntimeError(message)

def partitions(stubs,blocks=()):
    """Restricted partitions, canonically ordered by first appearance."""
    if not stubs:
        yield blocks
        return
    stub=stubs[0];neighbor=stub[0]
    for i,block in enumerate(blocks):
        if len(block)<3 and all(old[0]!=neighbor for old in block):
            yield from partitions(stubs[1:],blocks[:i]+(block+(stub,),)+blocks[i+1:])
    yield from partitions(stubs[1:],blocks+((stub,),))

def boundary_choices(size):
    # Flags are (0,0),(1,0),(1,1); if boundary degree equals
    # block size, only (0,0),(1,1) are allowed.
    return sum(2 if degree==size else 3 for degree in (2,3) if degree>=size)

def main():
    require([boundary_choices(k) for k in (1,2,3)]==[6,5,2],
            'Unexpected boundary-degree/flag choices')
    rows=[];totals={}
    for d in (2,3):
        total=0;possible_edges=list(combinations(range(1,d+1),2))
        for bits in product((0,1),repeat=len(possible_edges)):
            edges=[edge for edge,bit in zip(possible_edges,bits) if bit]
            internal={u:sum(u in edge for edge in edges) for u in range(1,d+1)}
            for degree_tuple in product((2,3),repeat=d):
                remaining=[degree_tuple[u-1]-1-internal[u] for u in range(1,d+1)]
                if min(remaining)<0:
                    continue
                stubs=tuple((u,i) for u in range(1,d+1) for i in range(1,remaining[u-1]+1))
                partition_count=0;full_count=0
                for blocks in partitions(stubs):
                    partition_count+=1;ways=1
                    require(sum(map(len,blocks))==len(stubs),'Lost stub')
                    for block in blocks:
                        require(1<=len(block)<=3,'Invalid block size')
                        require(len({stub[0] for stub in block})==len(block),
                                'Repeated neighbor label within block')
                        ways*=boundary_choices(len(block))
                    full_count+=ways
                rows.append(dict(root_degree=d,neighbor_edges=edges,
                    neighbor_degrees=list(degree_tuple),remaining_stub_counts=remaining,
                    restricted_partitions=partition_count,
                    full_views_after_boundary_degree_and_flag_choices=full_count))
                total+=full_count
        totals[str(d)]=total
    require(totals=={'2':2713,'3':229103},'Guide count does not match independent grammar')
    report=dict(status='PASS',method='Distinct-stub first-appearance restricted partitions; independent standard-library implementation',
        imports_from_delivered_generators=False,
        boundary_degree_and_flag_choices_by_block_size={'1':6,'2':5,'3':2},
        full_views_by_root_degree=totals,full_views_total=sum(totals.values()),
        structural_cases=rows,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Combinatorial multiset count only. No feature/sign replay or claim that counts prove universal graph coverage.')
    (HERE/'VIEW_GRAMMAR_COUNTS.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('PASS: independent grammar count gives 2,713 degree-two-root and 229,103 degree-three-root full views; total 231,816.')

if __name__=='__main__':main()
