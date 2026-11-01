#!/usr/bin/env python3
"""Independent structural count: no import or execution of primary generator.

Count boundary-decorated partial matchings for degree-two roots and
boundary-decorated set partitions on unmatched neighbor sets for high roots.
Exact integer recursion; no solver, graph fixture expansion or activities.
"""
if not __debug__:
    raise RuntimeError('Exact audit refuses optimized Python')
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
import hashlib
import gzip
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def decorations(k,incident_high):
    result=0
    for h in (2,3,4):
        if h<k or incident_high and h!=2:
            continue
        result+=2 if h==k else 3
    return result

# Labelled left stubs either have their own endpoint or join exactly one
# unconsumed right stub. This never constructs a full set partition.
def bipartite_stub_count(p,q,neighbor_edge):
    a=p-1-neighbor_edge
    b=q-1-neighbor_edge
    if min(a,b)<0 or neighbor_edge and min(p,q)>2:
        return 0
    L=decorations(1,p>2)
    R=decorations(1,q>2)
    P=decorations(2,p>2 or q>2)
    @lru_cache(None)
    def count(left,right_mask):
        if left==0:
            return R**right_mask.bit_count()
        value=L*count(left-1,right_mask)
        for j in range(b):
            if right_mask>>j&1:
                value+=P*count(left-1,right_mask^(1<<j))
        return value
    return count(a,(1<<b)-1)

# Pick the least remaining boundary stub; independently sum all possible
# blocks containing it, weighted by the possible degrees and flags.
@lru_cache(None)
def partition_weight(mask):
    if mask==0:return 1
    low=mask&-mask
    rest=mask^low
    value=0
    sub=rest
    while True:
        block=sub|low
        value+=decorations(block.bit_count(),False)*partition_weight(rest^sub)
        if sub==0:break
        sub=(sub-1)&rest
    return value

def neighbor_matchings(mask):
    if mask==0:
        yield (),0
        return
    low=mask&-mask
    i=low.bit_length()-1
    rest=mask^low
    # i has no inside edge, and contributes its one boundary stub.
    for edges,unmatched in neighbor_matchings(rest):
        yield edges,unmatched|low
    for j in range(rest.bit_length()):
        if rest>>j&1:
            for edges,unmatched in neighbor_matchings(rest^(1<<j)):
                yield ((i,j),)+edges,unmatched

by_tuple={}
by_edges=Counter()
for p,q in product((2,3,4),repeat=2):
    vals=[bipartite_stub_count(p,q,e) for e in (0,1)]
    by_tuple[f'2:{p}{q}']=sum(vals)
    for e,v in enumerate(vals):by_edges[f'2:{e}']+=v
for d in (3,4):
    total=0
    for edges,unmatched in neighbor_matchings((1<<d)-1):
        c=partition_weight(unmatched)
        total+=c
        by_edges[f'{d}:{len(edges)}']+=c
    by_tuple[f'{d}:'+str(2)*d]=total
by_root={str(d):sum(v for key,v in by_tuple.items() if key.startswith(f'{d}:'))
         for d in (2,3,4)}
assert by_root=={'2':5216,'3':977,'4':11360}
assert sum(by_root.values())==17553
primary=json.loads((ROOT/'enumeration.json').read_text())
assert primary['by_root']==by_root
assert primary['by_degree_tuple']==by_tuple
assert primary['by_neighbor_edge_count']==dict(by_edges)

# Prospective raw rows. Degree and spin are both explicit; impossible
# configurations may yield zero rows. Do not confuse raw size with rank.
states=list(product((2,3,4),(0,1)))
def allowed(a,b):return a==2 or b==2
stars=[]
for d in (2,3,4):
    for spin in (0,1):
        for neighbors in combinations_with_replacement(states,d):
            if all(allowed(d,h) for h,_ in neighbors):
                stars.append((d,spin,neighbors))
reversal=[]
for a,c in combinations(range(6),2):
    for b in range(6):
        for triangle in (0,1):
            if not allowed(states[a][0],states[b][0]):continue
            if not allowed(states[b][0],states[c][0]):continue
            if triangle and not allowed(states[a][0],states[c][0]):continue
            reversal.append((a,b,c,triangle))
classes=[(d,h,k) for d,h in product((2,3,4),repeat=2)
         for k in range(1,min(d,h)+1)]
patterns=[]
for d in (2,3,4):
    for neighbors in combinations_with_replacement((2,3,4),d):
        if all(allowed(d,h) for h in neighbors):patterns.append((d,neighbors))
assert (len(stars),len(reversal),len(classes),len(patterns))==(60,56,23,8)
rows=dict(rooted_stars=len(stars),two_step_reversal=len(reversal),
          degree_m_q=6,mixed_edge=2,equalities=60+56+6+2,
          boundary_classes=len(classes),exterior_rows=2*len(classes),
          association_patterns=len(patterns),
          exterior_rows_identically_zero_from_h_eq_k=sum(h==k for d,h,k in classes))
# Independently inspect the stored descriptions and evaluate each forbidden
# closed-neighborhood clause with ordinary sets rather than primary masks.
stored=0
stored_histogram=Counter()
seen=set()
with gzip.open(ROOT/'views.jsonl.gz','rt') as source:
    for line in source:
        v=json.loads(line)
        d=v['d'];deg=v['degrees'];edges=[set(e) for e in v['edges']]
        assert len(deg)==d and d in (2,3,4)
        assert all(allowed(d,h) for h in deg)
        assert all(len(e)==2 and e <= set(range(1,d+1)) for e in edges)
        assert len({tuple(sorted(e)) for e in edges})==len(edges)
        assert all(allowed(deg[a-1],deg[b-1]) for a,b in v['edges'])
        rem=[h-1-sum(i in e for e in edges) for i,h in enumerate(deg,1)]
        assert min(rem)>=0
        stubs=[i for i,r in enumerate(rem,1) for _ in range(r)]
        assert stubs==v['stubs']
        parts=v['parts']
        assert sorted(j for part in parts for j in part)==list(range(len(stubs)))
        assert all(part==sorted(part) and part for part in parts)
        assert [part[0] for part in parts]==sorted(part[0] for part in parts)
        blocks=[[stubs[j] for j in part] for part in parts]
        assert blocks==v['blocks']
        assert len(v['choices'])==len(blocks)
        for T,(h,x,b) in zip(blocks,v['choices']):
            assert len(set(T))==len(T)
            assert h in (2,3,4) and h>=len(T)
            assert all(allowed(h,deg[i-1]) for i in T)
            assert (x,b) in ((0,0),(1,0),(1,1))
            assert h!=len(T) or x==b
        forbidden=[set(range(d+1))]
        for i in range(1,d+1):
            exterior_selected=all(x for T,(h,x,b) in zip(blocks,v['choices']) if i in T)
            if exterior_selected:
                C={0,i}
                for e in edges:
                    if i in e:C|=e
                forbidden.append(C)
        forbidden += [set(T) for T,(h,x,b) in zip(blocks,v['choices']) if b]
        expected=[]
        for k in range(1<<(d+1)):
            S={i for i in range(d+1) if k>>i&1}
            if all(not C<=S for C in forbidden):expected.append(k)
        assert expected==v['masks']
        key=(d,tuple(deg),tuple(map(tuple,v['edges'])),tuple(map(tuple,parts)),tuple(map(tuple,v['choices'])))
        assert key not in seen
        seen.add(key)
        stored+=1
        stored_histogram[str(d)]+=1
assert stored==17553 and dict(stored_histogram)==by_root
result=dict(all_passed=True,total=sum(by_root.values()),by_root=by_root,
            stored_descriptions_independently_checked=stored,
            by_degree_tuple=by_tuple,by_neighbor_edge_count=dict(by_edges),
            independent_partition_weights=[partition_weight((1<<k)-1) for k in range(5)],
            row_dimensions=rows,
            scope='Independent count and prospective raw dimensions only; no certificate or rank computation.',
            enumeration_rules_sha256=hashlib.sha256((ROOT/'PROTOCOL.md').read_bytes()).hexdigest(),
            primary_generator_sha256=hashlib.sha256((ROOT/'enumerate_validate.py').read_bytes()).hexdigest(),
            primary_count_sha256=hashlib.sha256((ROOT/'enumeration.json').read_bytes()).hexdigest(),
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'grammar_count_independent.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,indent=2,sort_keys=True))
