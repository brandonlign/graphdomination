#!/usr/bin/env python3
"""Independent ordinary-set validation of all graph/exposure families."""
import sys
if not __debug__: raise RuntimeError("Optimized Python is forbidden for V3 exact checks")
from itertools import combinations
from fractions import Fraction as F
from pathlib import Path
import json

OUT=Path(__file__).resolve().parent

def subsets(items):
    a=tuple(sorted(items))
    for r in range(len(a)+1):
        for t in combinations(a,r):
            yield frozenset(t)

def graph(n,edges):
    A={i:set() for i in range(n)}
    for u,v in edges:
        assert u!=v and v not in A[u]
        A[u].add(v); A[v].add(u)
    assert all(len(A[v]) in (2,3) for v in A)
    return A

GRAPHS={
'C4':graph(4,[(0,1),(1,2),(2,3),(3,0)]),
'K2_3':graph(5,[(u,v) for u in (0,1) for v in (2,3,4)]),
'K4':graph(4,list(combinations(range(4),2))),
'Q3':graph(8,[(u,v) for u,v in combinations(range(8),2) if (u^v) in (1,2,4)]),
'K4_one_subdivision':graph(5,[e for e in combinations(range(4),2) if e!=(0,1)]+[(0,4),(1,4)]),
}

def allowed(A,I):
    return all(not ({z}|A[z])<=I for z in A)

def view_data(A,v,E):
    U=sorted(A[v]); R={v}|set(U)
    B=sorted(set().union(*(A[u]-R for u in U)))
    index={u:i+1 for i,u in enumerate(U)}
    induced=[(index[u],index[w]) for u,w in combinations(U,2) if w in A[u]]
    blocks=[]
    for w in B:
        T=A[w]&R
        assert v not in T and T and len(T)<=3
        x=int(w in E)
        b=int(({w}|(A[w]-R))<=E)
        assert b<=x
        if len(A[w])==len(T): assert b==x
        blocks.append({'vertex':w,'T':sorted(index[u] for u in T),'degree':len(A[w]),'x':x,'b':b})
    assert sum(len(q['T']) for q in blocks)==sum(len(A[u])-1-len(A[u]&set(U)) for u in U)
    return {'root':v,'d':len(U),'neighbors':U,'neighbor_degrees':[len(A[u]) for u in U],
            'induced_edges':induced,'blocks':blocks}

def view_allowed(A,V,L,omit_b=False,omit_neighbor=False):
    v=V['root']; U=V['neighbors']; R={v}|set(U)
    if R<=L: return False
    if not omit_neighbor:
        for u in U:
            needed=[q for q in V['blocks'] if V['neighbors'].index(u)+1 in q['T']]
            if all(q['x'] for q in needed) and ({v,u}|(A[u]&set(U)))<=L:
                return False
    if not omit_b:
        for q in V['blocks']:
            if q['b'] and {U[i-1] for i in q['T']}<=L:
                return False
    return True

def fmt(x): return str(x)

def run():
    results=[]; all_views=[]; failed_controls={'omit_b':[],'omit_neighbor':[]}
    for name,A in GRAPHS.items():
        n=len(A); feasible=[I for I in subsets(A) if allowed(A,I)]; Z=len(feasible)
        coeff=[sum(len(I)==r for I in feasible) for r in range(n+1)]
        sums=[F(0) for _ in range(6)] # norm, degree2,degree3,occ2,occ3,joint
        exposure_count=0; comparisons=0
        for v in A:
            R={v}|A[v]; outside=set(A)-R
            for E in subsets(outside):
                actual=[L for L in subsets(R) if allowed(A,E|L)]
                if not actual: continue
                exposure_count+=1
                V=view_data(A,v,E)
                predicted=[L for L in subsets(R) if view_allowed(A,V,L)]
                assert actual==predicted,(name,v,E,actual,predicted)
                comparisons+=len(list(subsets(R)))
                for control in failed_controls:
                    bad=[L for L in subsets(R) if view_allowed(A,V,L,**{control:True})]
                    if bad!=actual and not any(w['graph']==name for w in failed_controls[control]):
                        failed_controls[control].append({'graph':name,'root':v,'outside':sorted(E),
                            'extra_subsets':[sorted(L) for L in bad if L not in actual]})
                weight=F(len(actual),n*Z); d=len(A[v]); rootocc=F(sum(v in L for L in actual),len(actual))
                fs=[F(1)]
                for k in (2,3): fs.append(F(int(d==k))-sum(F(1,k) for u in A[v] if len(A[u])==k))
                for k in (2,3):
                    fs.append(int(d==k)*rootocc-sum(F(sum(u in L for L in actual),len(actual)*k) for u in A[v] if len(A[u])==k))
                joint=F(0)
                for u in A[v]:
                    sign=int(d==2 and len(A[u])==3)-int(d==3 and len(A[u])==2)
                    joint+=sign*F(sum(v in L and u in L for L in actual),len(actual))
                fs.append(joint)
                sums=[a+weight*b for a,b in zip(sums,fs)]
                V.update({'graph':name,'outside':sorted(E),'weight':fmt(weight),'root_occupancy':fmt(rootocc),
                          'features':[fmt(x) for x in fs],
                          'feasible_subsets':[[0 if x==v else V['neighbors'].index(x)+1 for x in sorted(L)] for L in actual]})
                all_views.append(V)
        assert sums==[1,0,0,0,0,0],(name,sums)
        density=F(sum(map(len,feasible)),n*Z)
        results.append({'graph':name,'n':n,'edges':[[u,w] for u in A for w in sorted(A[u]) if u<w],
            'complement_coefficients':coeff,'complement_count':Z,'complement_size_sum':sum(map(len,feasible)),
            'density':fmt(density),'exposures':exposure_count,'subset_comparisons':comparisons,
            'consistency_sums':[fmt(x) for x in sums]})
        if name=='C4': assert density==F(4,11)
    assert all(failed_controls.values()),failed_controls
    data={'all_passed':True,'protocol':'PROTOCOL.md','graphs':results,'negative_controls':failed_controls,
          'total_feasible_root_exposures':len(all_views),'views_file':'graph_views.json'}
    (OUT/'validation_results.json').write_text(json.dumps(data,indent=2)+'\n')
    (OUT/'graph_views.json').write_text(json.dumps(all_views,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__':run()
