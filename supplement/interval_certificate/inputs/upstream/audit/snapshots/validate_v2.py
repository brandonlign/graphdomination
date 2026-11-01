#!/usr/bin/env python3
"""Independent ordinary-set graph/exposure validation of all128 V2 equations."""
import sys
if not __debug__: raise RuntimeError("Optimized Python is forbidden for V3 exact checks")
from itertools import product
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import json
from base_views_v1_copy import GRAPHS,subsets,allowed,view_data,view_allowed
HERE=Path(__file__).resolve().parent
PATTERNS=tuple(product((0,1),range(4),range(4),range(4)))

def symbol(degree,bit):
    assert degree in (2,3) and bit in (0,1)
    return 2*(degree-2)+int(bit)

def direct_counts(A,v,I):
    out=Counter()
    for u in A[v]:
        for w in A[u]-{v}:
            out[(int(w in A[v]),symbol(len(A[v]),int(v in I)),
                 symbol(len(A[u]),int(u in I)),symbol(len(A[w]),int(w in I)))]+=1
    return out

def view_counts(V,S):
    """S uses canonical local labels0,...,d."""
    out=Counter();d=V['d'];deg=V['neighbor_degrees'];edges=set(tuple(e) for e in V['induced_edges'])
    s0=symbol(d,int(0 in S))
    for u in range(1,d+1):
        s1=symbol(deg[u-1],int(u in S))
        for w in range(1,d+1):
            if (min(u,w),max(u,w)) in edges:
                out[(1,s0,s1,symbol(deg[w-1],int(w in S)))]+=1
        for B in V['blocks']:
            if u in B['T']:
                out[(0,s0,s1,symbol(B['degree'],B['x']))]+=1
    return out

def reverse(K):
    a,b,c,d=K
    return a,d,c,b

def original_features(V,sets):
    d=V['d'];deg=V['neighbor_degrees'];z=len(sets)
    p=[F(sum(i in S for S in sets),z) for i in range(d+1)]
    out=[p[0]]
    for k in (2,3):out.append(F(d==k)-sum(F(1,k) for h in deg if h==k))
    for k in (2,3):out.append(int(d==k)*p[0]-sum(p[i]/k for i in range(1,d+1) if deg[i-1]==k))
    j=F(0)
    for u in range(1,d+1):
        sign=int(d==2 and deg[u-1]==3)-int(d==3 and deg[u-1]==2)
        j+=sign*F(sum(0 in S and u in S for S in sets),z)
    return out+[j]

def run():
    graphresults=[];records=[]
    for name,A in GRAPHS.items():
        n=len(A);allsets=[I for I in subsets(A) if allowed(A,I)];Z=len(allsets)
        totals={K:F(0) for K in PATTERNS};old=[F(0) for j in range(6)];count=0
        for v in A:
            R={v}|A[v]
            for E in subsets(set(A)-R):
                local=[S for S in subsets(R) if allowed(A,E|S)]
                if not local:continue
                count+=1;V=view_data(A,v,E);z=len(local)
                assert local==[S for S in subsets(R) if view_allowed(A,V,S)]
                canonical=[frozenset(0 if x==v else V['neighbors'].index(x)+1 for x in S) for S in local]
                cnt=Counter()
                for S,C in zip(local,canonical):
                    dcount=direct_counts(A,v,E|S);vcount=view_counts(V,C)
                    assert dcount==vcount,(name,v,E,S,dcount,vcount)
                    cnt.update(vcount)
                residuals={K:F(cnt[K]-cnt[reverse(K)],z) for K in PATTERNS}
                weight=F(z,n*Z)
                for K in PATTERNS:totals[K]+=weight*residuals[K]
                f=original_features(V,canonical)
                old[0]+=weight
                for j in range(1,6):old[j]+=weight*f[j]
                records.append({'graph':name,'view':V,'outside':sorted(E),'weight':str(weight),
                    'local_subsets':[sorted(C) for C in canonical],'original_features':list(map(str,f)),
                    'pattern_residuals':[str(residuals[K]) for K in PATTERNS]})
        assert all(x==0 for x in totals.values()),(name,totals)
        assert old==[1,0,0,0,0,0],(name,old)
        density=F(sum(map(len,allsets)),n*Z)
        if name=='C4':assert density==F(4,11)
        graphresults.append({'graph':name,'n':n,'count':Z,'density':str(density),
            'feasible_root_exposures':count,'all128_reversal_sums_zero':True,
            'original_consistency_sums':list(map(str,old))})
    out={'all_passed':True,'protocol':'PROTOCOL_V2_AUDIT.md','graphs':graphresults,
         'ordered_pattern_equations':len(PATTERNS),'actual_exposures':len(records),
         'symbol_code':'2*(degree-2)+selected_bit','pattern_order':[list(K) for K in PATTERNS]}
    (HERE/'validation_v2.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'graph_views_v2.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='pattern_order'},indent=2))
if __name__=='__main__':run()
