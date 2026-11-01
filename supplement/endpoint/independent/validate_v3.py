#!/usr/bin/env python3
"""Validate full V3 family on actual graph/exposure distributions."""
import sys
if not __debug__: raise RuntimeError('Optimized Python is forbidden for V3 exact checks')
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import json
from base_views_v1_copy import GRAPHS,subsets,allowed,view_data,view_allowed
from validate_v2 import PATTERNS as WALK_PATTERNS, direct_counts,reverse,symbol
from enumerate_v2_independent import COORDS
from features_v3_independent import SIGS,CLASSES,PATTERNS,signatures,from_actual,canonical
HERE=Path(__file__).resolve().parent

def actual_signature(A,u,I):
    return len(A[u]),int(u in I),tuple(sorted(symbol(len(A[w]),int(w in I)) for w in A[u]))

def run():
    records=[];results=[];c4=[]
    for name,A in GRAPHS.items():
        n=len(A);feasible=[I for I in subsets(A) if allowed(A,I)];Z=len(feasible)
        sums=[F(0) for _ in range(141)];normalization=F(0);cases=0
        for v in A:
            R={v}|A[v]
            for E in subsets(set(A)-R):
                local=[S for S in subsets(R) if allowed(A,E|S)]
                if not local:continue
                cases+=1;V=view_data(A,v,E)
                assert local==[S for S in subsets(R) if view_allowed(A,V,S)]
                fs=from_actual(V);weight=F(len(local),n*Z)
                # Verify root and every neighbor's full signature against original graph.
                counts=Counter();star={s:F(0) for s in SIGS}
                blocks=tuple((tuple(B['T']),B['degree'],B['x'],B['b']) for B in V['blocks'])
                for S in local:
                    I=E|S;C=frozenset(0 if x==v else V['neighbors'].index(x)+1 for x in S)
                    root,neighbors=signatures(V['d'],V['neighbor_degrees'],set(tuple(e) for e in V['induced_edges']),blocks,C)
                    assert root==actual_signature(A,v,I)
                    assert neighbors==[actual_signature(A,u,I) for u in V['neighbors']]
                    star[root]+=1
                    for s in neighbors:star[s]-=F(1,s[0])
                    counts.update(direct_counts(A,v,I))
                assert fs[6:54]==tuple(F(counts[K]-counts[reverse(K)],len(local)) for K in COORDS)
                assert fs[54:114]==tuple(star[s]/len(local) for s in SIGS)
                # Independently check structural block counts using actual outside vertices.
                boundary=set().union(*(A[u]-R for u in A[v]));g=[]
                for d,h,k in CLASSES:
                    ws=[w for w in boundary if len(A[v])==d and len(A[w])==h and len(A[w]&R)==k]
                    sx=sum(w in E for w in ws)
                    sb=sum(({w}|(A[w]-R))<=E for w in ws)
                    g.extend((F(len(ws),2)-sx,F(sx,2**(h-k))-sb))
                assert fs[114:134]==tuple(g)
                p=F(sum(v in S for S in local),len(local));g=[]
                for d,r,L in PATTERNS:
                    ar=sum(len(A[u])==(3 if d==2 else 2) for u in A[v])
                    g.append(p-L if len(A[v])==d and ar==r else F(0))
                assert fs[134:]==tuple(g)
                normalization+=weight;sums=[a+weight*b for a,b in zip(sums,fs)]
                item={'graph':name,'view':V,'outside':sorted(E),'weight':str(weight),'features':list(map(str,fs))}
                records.append(item)
                if name=='C4':c4.append(item)
        assert normalization==1 and all(x==0 for x in sums[1:114])
        assert all(x>=0 for x in sums[114:]),(name,sums[114:])
        density=F(sum(map(len,feasible)),n*Z);assert sums[0]==density
        if name=='C4':assert density==F(4,11)
        results.append({'graph':name,'n':n,'admissible_count':Z,'feasible_root_exposures':cases,
            'density':str(density),'normalization':str(normalization),'all113_equalities_zero':True,
            'inequality_sums':list(map(str,sums[114:]))})
    out={'all_passed':True,'protocol':'PROTOCOL_V3_AUDIT.md','graphs':results,'actual_exposures':len(records),
         'equalities_excluding_normalization':113,'inequalities':27,'c4_optimum_candidate':'4/11','c4_feasible':True,
         'star_signature_order':[[d,x,list(ns)] for d,x,ns in SIGS],
         'exterior_class_order':[list(x) for x in CLASSES],
         'association_order':[[d,r,str(L)] for d,r,L in PATTERNS]}
    (HERE/'validation_v3.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'graph_views_v3.json').write_text(json.dumps(records,indent=2)+'\n')
    (HERE/'c4_feasible_mixture_v3.json').write_text(json.dumps(c4,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if not k.endswith('_order')},indent=2))
if __name__=='__main__':run()
