#!/usr/bin/env python3
"""Exact fixed-fixture checks only; no graph census or optimization.

K2,4 and K2,5 are the saved family. The final three records preserve a
specific failed max-degree-four construction inspired by splitting K2,5's
hubs; they are not a search result or proof of monotonicity.
"""
if not __debug__:
    raise RuntimeError('Exact fixture audit refuses optimized Python')

from collections import Counter
from fractions import Fraction
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent

def polynomial(n,edges):
    neighborhoods=[1<<v for v in range(n)]
    for u,v in edges:
        assert u!=v
        neighborhoods[u]|=1<<v
        neighborhoods[v]|=1<<u
    counts=Counter(S.bit_count() for S in range(1<<n)
                   if all(S&N for N in neighborhoods))
    return [counts[k] for k in range(n+1)]

def degrees(n,edges):
    return [sum(v in edge for edge in edges) for v in range(n)]

def mean(p):
    return Fraction(sum(k*c for k,c in enumerate(p)),sum(p))

def check(n,edges,added):
    assert added not in edges and tuple(reversed(added)) not in edges
    p=polynomial(n,edges)
    q=polynomial(n,edges+[added])
    dh=degrees(n,edges);dg=degrees(n,edges+[added])
    assert min(dh)>=2
    assert all(qk>=pk for pk,qk in zip(p,q))
    return dict(n=n,edges=edges,added_edge=added,H_degrees=dh,G_degrees=dg,
                H_edge_minimal=all(dh[u]==2 or dh[v]==2 for u,v in edges),
                D_H=p,D_G=q,mean_H_at_one=str(mean(p)),
                mean_G_at_one=str(mean(q)),mean_difference_at_one=str(mean(q)-mean(p)))

results=[]
for k in (4,5):
    E=[(h,l) for h in (0,1) for l in range(2,k+2)]
    item=check(k+2,E,(2,3))
    item['name']=f'K2,{k} plus one edge in the k-part'
    diff=[q-p for p,q in zip(item['D_H'],item['D_G'])]
    assert diff==[2 if j==k-1 else 0 for j in range(k+3)]
    results.append(item)
assert results[0]['mean_difference_at_one']=='-6/329'
assert results[1]['mean_difference_at_one']=='18/9215'
# Symbolic K2,4 numerator 3D - tD': no numerical activities chosen.
p=results[0]['D_H']
numerator=[(3-k)*c for k,c in enumerate(p)]
assert numerator==[0,0,9,0,-15,-12,-3]
# One explicit degree-four construction, with its three natural leaf-edge
# variants preserved after all failed to increase the activity-one mean.
E=[(0,4),(4,1),(2,5),(5,3)]
E += [(0,i) for i in (6,7)]+[(1,i) for i in (8,9,10)]
E += [(2,i) for i in (6,7)]+[(3,i) for i in (8,9,10)]
for e in [(8,9),(6,7),(6,8)]:
    item=check(11,E,e)
    assert max(item['H_degrees'])==4 and max(item['G_degrees'])==4
    assert item['H_edge_minimal']
    assert Fraction(item['mean_difference_at_one'])<0
    item['name']='split-hub K2,5-inspired negative construction'
    results.append(item)
# Explicit equality/strictness fixtures, chosen from the analytic proof.
for n,E,e,label,equal in [
    (6,[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3)],(0,3),'two triangles plus a cross edge',True),
    (4,[(0,1),(1,2),(2,3),(3,0)],(0,2),'C4 plus a diagonal',False),
    (8,[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4)],(0,4),'two C4s plus a cross edge',False)
]:
    item=check(n,E,e)
    assert (item['D_H']==item['D_G'])==equal
    item['name']=label
    results.append(item)
output=dict(all_passed=True,scope='Eight fixed graph pairs, no census; all subsets counted exactly.',
            K24_3D_minus_tDprime=numerator,fixtures=results,
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'transfer_fixture_checks.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(dict(all_passed=True,fixtures=len(results),
    mean_differences=[r['mean_difference_at_one'] for r in results]),indent=2))
