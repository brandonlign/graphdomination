#!/usr/bin/env python3
"""Independent V3 ordinary-set features, with all constraints explicitly ordered."""
import sys
if not __debug__: raise RuntimeError('Optimized Python is forbidden for V3 exact checks')
from itertools import combinations_with_replacement
from fractions import Fraction as F
from math import gcd,lcm
from functools import reduce
from collections import Counter
from enumerate_v2_independent import feature_key as v2_key, family
from validate_v2 import symbol
SIGS=tuple((d,x,nb) for d in (2,3) for x in (0,1) for nb in combinations_with_replacement(range(4),d))
CLASSES=tuple((d,h,k) for d in (2,3) for h in (2,3) for k in range(1,h+1))
PATTERNS=((2,0,F(1,3)),(2,1,F(9,25)),(2,2,F(5,13)),(3,0,F(79,207)),(3,1,F(69,197)),(3,2,F(15,47)),(3,3,F(13,45)))

def signatures(d,deg,E,blocks,S):
    root=(d,int(0 in S),tuple(sorted(symbol(deg[u-1],int(u in S)) for u in range(1,d+1))))
    neighbors=[]
    for u in range(1,d+1):
        ns=[symbol(d,int(0 in S))]
        ns.extend(symbol(deg[w-1],int(w in S)) for w in range(1,d+1) if (min(u,w),max(u,w)) in E)
        ns.extend(symbol(h,x) for T,h,x,b in blocks if u in T)
        assert len(ns)==deg[u-1]
        neighbors.append((deg[u-1],int(u in S),tuple(sorted(ns))))
    return root,neighbors

def feature_vector(d,deg,E,blocks):
    old=v2_key(d,deg,E,blocks);features=[F(n,old[-1]) for n in old[:-1]]
    sets=family(d,E,blocks);z=len(sets);star={s:F(0) for s in SIGS}
    for S in sets:
        root,neighbors=signatures(d,deg,E,blocks,S);star[root]+=1
        for sig in neighbors:star[sig]-=F(1,sig[0])
    features.extend(star[s]/z for s in SIGS)
    for dr,h,k in CLASSES:
        actual=[(T,hh,x,b) for T,hh,x,b in blocks if d==dr and hh==h and len(T)==k]
        sx=sum(x for T,hh,x,b in actual);sb=sum(b for T,hh,x,b in actual)
        features.extend((F(len(actual),2)-sx,F(sx,2**(h-k))-sb))
    p0=features[0]
    for dr,r,L in PATTERNS:
        actualr=sum(h==(3 if dr==2 else 2) for h in deg)
        features.append((p0-L) if d==dr and actualr==r else F(0))
    assert len(features)==141
    return tuple(features)

def canonical(features):
    D=lcm(*(q.denominator for q in features))
    nums=tuple(q.numerator*(D//q.denominator) for q in features)
    g=reduce(gcd,nums,D)
    return tuple(n//g for n in nums)+(D//g,)

def from_actual(V):
    return feature_vector(V['d'],tuple(V['neighbor_degrees']),frozenset(tuple(e) for e in V['induced_edges']),
          tuple(sorted((tuple(B['T']),B['degree'],B['x'],B['b']) for B in V['blocks'])))

if __name__=='__main__':
    print('Feature module; 1 objective,113 zero-mean coordinates,27 nonnegative coordinates.')
