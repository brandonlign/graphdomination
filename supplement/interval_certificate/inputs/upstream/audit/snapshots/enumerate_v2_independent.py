#!/usr/bin/env python3
"""Independent full-degree V2 universe, ordinary sets, exact integer columns."""
import sys
if not __debug__: raise RuntimeError("Optimized Python is forbidden for V3 exact checks")
from itertools import product,combinations
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
from math import gcd
from functools import reduce
import json,gzip,time
from base_views_v1_copy import subsets
from validate_v2 import view_counts,reverse
HERE=Path(__file__).resolve().parent
COORDS=tuple((t,a,b,c) for t in (0,1) for a in range(4) for b in range(4) for c in range(a+1,4))

def partitions(items):
    if not items:
        yield ();return
    first=items[0];others=items[1:]
    for r in (0,1,2):
        for extra in combinations(others,r):
            B=(first,)+extra
            if len({u for u,j in B})<len(B):continue
            for tail in partitions(tuple(x for x in others if x not in extra)):
                yield (B,)+tail

def family(d,edges,blocks):
    U=set(range(1,d+1));R={0}|U;forbidden=[R]
    for u in U:
        if all(B[2] for B in blocks if u in B[0]):
            forbidden.append({0,u}|{w for w in U if (min(u,w),max(u,w)) in edges})
    forbidden.extend(set(B[0]) for B in blocks if B[3])
    return tuple(S for S in subsets(R) if all(not T<=S for T in forbidden))

def feature_key(d,deg,edges,blocks):
    sets=family(d,edges,blocks);z=len(sets)
    counts=[sum(u in S for S in sets) for u in range(d+1)]
    nums=[6*counts[0]]
    for k in (2,3):nums.append((6//k)*z*(k*int(d==k)-sum(h==k for h in deg)))
    for k in (2,3):nums.append((6//k)*(k*int(d==k)*counts[0]-sum(counts[u] for u in range(1,d+1) if deg[u-1]==k)))
    j=0
    for u in range(1,d+1):
        sign=int(d==2 and deg[u-1]==3)-int(d==3 and deg[u-1]==2)
        j+=sign*sum(0 in S and u in S for S in sets)
    nums.append(6*j)
    V={'d':d,'neighbor_degrees':deg,'induced_edges':edges,
       'blocks':[{'T':B[0],'degree':B[1],'x':B[2],'b':B[3]} for B in blocks]}
    cts=Counter()
    for S in sets:cts.update(view_counts(V,S))
    nums.extend(6*(cts[K]-cts[reverse(K)]) for K in COORDS)
    denominator=6*z;g=reduce(gcd,nums,denominator)
    return tuple(n//g for n in nums)+(denominator//g,)

def actual_key(V):
    return feature_key(V['d'],tuple(V['neighbor_degrees']),frozenset(tuple(e) for e in V['induced_edges']),
               tuple(sorted((tuple(B['T']),B['degree'],B['x'],B['b']) for B in V['blocks'])))

def run():
    start=time.time();cache={};mult=Counter();full_count=0;bydegree=Counter()
    for d in (2,3):
        for E in subsets(combinations(range(1,d+1),2)):
            for deg in product((2,3),repeat=d):
                amounts=[deg[u-1]-1-sum(u in e for e in E) for u in range(1,d+1)]
                if min(amounts)<0:continue
                stubs=tuple((u,j) for u,a in enumerate(amounts,1) for j in range(a))
                for P in partitions(stubs):
                    T=tuple(tuple(u for u,j in B) for B in P)
                    opts=[[(h,x,b) for h in (2,3) if h>=len(B) for x,b in ((0,0),(1,0),(1,1)) if h!=len(B) or x==b] for B in T]
                    for choice in product(*opts):
                        blocks=tuple(sorted((B,*C) for B,C in zip(T,choice)))
                        K=(d,deg,E,blocks)
                        if K not in cache:cache[K]=feature_key(*K)
                        mult[cache[K]]+=1;full_count+=1;bydegree[d]+=1
    assert full_count==231816
    records=json.loads((HERE/'graph_views_v2.json').read_text());coordindex={tuple(K):i for i,K in enumerate(product((0,1),range(4),range(4),range(4)))}
    for R in records:
        K=actual_key(R['view']);assert K in mult
        want=tuple(map(F,R['original_features']))+tuple(F(R['pattern_residuals'][coordindex[C]]) for C in COORDS)
        assert tuple(F(n,K[-1]) for n in K[:-1])==want
    data={'all_passed':True,'independent_full_views':full_count,'full_views_by_root_degree':dict(bydegree),
          'unique_combinatorial_views':len(cache),'unique_columns':len(mult),
          'actual_exposures_checked':len(records),'coordinate_order':[list(C) for C in COORDS],
          'rows':[{'denominator':K[-1],'numerators':list(K[:-1]),'full_views':n} for K,n in sorted(mult.items())]}
    with gzip.open(HERE/'independent_columns_v2.json.gz','wt') as f:json.dump(data,f,sort_keys=True)
    summary={k:v for k,v in data.items() if k not in ('rows','coordinate_order')};summary['elapsed_seconds']=time.time()-start
    (HERE/'independent_enumeration_v2.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':run()
