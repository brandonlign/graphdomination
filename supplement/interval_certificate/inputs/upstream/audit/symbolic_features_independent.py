#!/usr/bin/env python3
"""Independent exact cardinality-indexed full weighted feature polynomials.

No primary interval evaluator imported. Ordinary-set feasible families and
walks use the previously independently audited structural implementation.
"""
from __future__ import annotations
if not __debug__: raise RuntimeError('Independent continuum checks refuse optimized Python')
from pathlib import Path
from itertools import combinations, combinations_with_replacement, product
from fractions import Fraction as F
from collections import Counter
from math import comb
import sys, json, gzip, hashlib, time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'snapshots'))
from base_views_v1_copy import subsets
from enumerate_v2_independent import family, partitions
from validate_v2 import symbol, view_counts, reverse
SIGS=tuple((d,x,nb)for d in(2,3)for x in(0,1)for nb in combinations_with_replacement(range(4),d))
SI={s:i for i,s in enumerate(SIGS)}
CLASSES=tuple((d,h,k)for d in(2,3)for h in(2,3)for k in range(1,h+1))
PATTERNS=tuple([(2,r)for r in range(3)]+[(3,r)for r in range(4)])
COORDS=tuple((t,a,b,c)for t in(0,1)for a in range(4)for b in range(4)for c in range(a+1,4))
def require(p,msg):
    if not p: raise RuntimeError(msg)
def trim(v):
    v=list(v)
    while len(v)>1 and v[-1]==0:v.pop()
    return tuple(v) or (0,)
def peval(v,z):
    ans=F(0)
    for c in reversed(v):ans=ans*z+c
    return ans

def signatures(d,deg,edges,blocks,S):
    rt=(d,int(0 in S),tuple(sorted(symbol(deg[u-1],u in S)for u in range(1,d+1))))
    ns=[]
    for u in range(1,d+1):
        sym=[symbol(d,0 in S)]
        sym.extend(symbol(deg[v-1],v in S)for v in range(1,d+1)if(min(u,v),max(u,v))in edges)
        sym.extend(symbol(h,x)for T,h,x,b in blocks if u in T)
        require(len(sym)==deg[u-1],'Star size mismatch')
        ns.append((deg[u-1],int(u in S),tuple(sorted(sym))))
    return rt,ns

def symbolic(d,deg,edges,blocks):
    sets=family(d,edges,blocks)
    require(bool(sets)and frozenset() in sets,'Missing empty local set')
    size=d+2
    Z=[0]*size;O=[0]*size
    polys=[[0]*size for _ in range(113)]
    V=dict(d=d,neighbor_degrees=deg,induced_edges=edges,blocks=[dict(T=T,degree=h,x=x,b=b)for T,h,x,b in blocks])
    for S in sets:
        n=len(S);Z[n]+=1;O[n]+=0 in S
        fs=[]
        for k in(2,3):fs.append(6*int(d==k)-(6//k)*sum(h==k for h in deg))
        for k in(2,3):fs.append(6*int(d==k)*int(0 in S)-(6//k)*sum(u in S for u in range(1,d+1)if deg[u-1]==k))
        joint=sum((int(d==2 and deg[u-1]==3)-int(d==3 and deg[u-1]==2))*int(0 in S and u in S)for u in range(1,d+1))
        fs.append(6*joint)
        counts=view_counts(V,S)
        fs.extend(6*(counts[K]-counts[reverse(K)])for K in COORDS)
        star=[0]*60
        root,neighbors=signatures(d,deg,edges,blocks,S)
        star[SI[root]]+=6
        for sig in neighbors:star[SI[sig]]-=6//sig[0]
        fs.extend(star)
        require(len(fs)==113,'Wrong transport count')
        for row,x in zip(polys,fs):row[n]+=x
    exterior=[]
    for rd,h,k in CLASSES:
        selected=[B for B in blocks if d==rd and B[1]==h and len(B[0])==k]
        exterior.append((len(selected),sum(B[2]for B in selected),sum(B[3]for B in selected)))
    r=sum(h==(3 if d==2 else 2)for h in deg)
    return (d,trim(Z),trim(O),tuple(trim(p)for p in polys),tuple(exterior),PATTERNS.index((d,r)))

def enumerate_inputs():
    for d in(2,3):
        for E in subsets(combinations(range(1,d+1),2)):
            for deg in product((2,3),repeat=d):
                amount=[deg[u-1]-1-sum(u in e for e in E)for u in range(1,d+1)]
                if min(amount)<0:continue
                stubs=tuple((u,j)for u,a in enumerate(amount,1)for j in range(a))
                for partition in partitions(stubs):
                    T=tuple(tuple(u for u,j in B)for B in partition)
                    opts=[[(h,x,b)for h in(2,3)if h>=len(B)for x,b in((0,0),(1,0),(1,1))if h!=len(B)or x==b]for B in T]
                    for choice in product(*opts):
                        blocks=tuple(sorted((B,*C)for B,C in zip(T,choice)))
                        yield d,deg,E,blocks

def from_primary_rep(rep):
    blocks=tuple(sorted((tuple(T),h,x,b)for T,h,(x,b)in zip(rep['blocks'],rep['boundary_degrees'],rep['flags'])))
    return symbolic(rep['d'],tuple(rep['neighbor_degrees']),frozenset(tuple(e)for e in rep['neighbor_edges']),blocks)

def independent_association():
    import sympy as sp
    z=sp.Symbol('z');q=z/(1+z);out=[]
    for d,r in PATTERNS:
        deg=([3]*r+[2]*(d-r))if d==2 else([2]*r+[3]*(d-r))
        T=0
        # Direct disjoint probabilities of proper selected-neighbor subsets.
        for S in subsets(range(d)):
            if len(S)==d:continue
            term=q**len(S)*(1-q)**(d-len(S))
            for u in S:term*=1-q**(deg[u]-1)
            T+=term
        n,den=sp.fraction(sp.cancel(z*T/(1+z*T)))
        if den.subs(z,0)<0:n=-n;den=-den
        cv=lambda p:trim([int(x)for x in reversed(sp.Poly(p,z).all_coeffs())])
        N,D=cv(n),cv(den)
        require(D[0]>0 and all(x>=0 for x in D),'Association denominator positivity coefficients')
        require(N[0]==0 and all(x>=0 for x in N),'Association numerator nonnegative coefficients')
        out.append(dict(pattern=[d,r],N=N,D=D))
    return out

def eval_complete(K,z,association=None):
    d,Z,O,T,ext,p=K
    Z=peval(Z,z);O=peval(O,z)
    out=[O/Z]+[peval(t,z)/(6*Z)for t in T]
    q=z/(1+z)
    for (rd,h,k),(N,X,B)in zip(CLASSES,ext):out.extend([q*N-X,q**(h-k)*X-B])
    aa=association or independent_association()
    for j,L in enumerate(aa):out.append(O/Z-peval(L['N'],z)/peval(L['D'],z)if j==p else F(0))
    require(len(out)==141,'Wrong full feature count')
    return tuple(out)

def run():
    start=time.time();cache={};mult=Counter();degrees=Counter()
    for count,K in enumerate(enumerate_inputs(),1):
        if K not in cache:cache[K]=symbolic(*K)
        mult[cache[K]]+=1;degrees[K[0]]+=1
        if count%50000==0:print('Independent symbolic full views',count,flush=True)
    require(count==231816,'Wrong independent structural count')
    rows=[dict(key=K,full_views=mult[K])for K in sorted(mult)]
    with gzip.open(HERE/'symbolic_columns_independent.json.gz','wt')as f:json.dump(dict(rows=rows),f,separators=(',',':'))
    association=independent_association()
    (HERE/'association_independent.json').write_text(json.dumps(association,indent=2)+'\n')
    report=dict(all_passed=True,full_views=count,full_views_by_degree=dict(degrees),combinatorial_views=len(cache),symbolic_columns=len(rows),seconds=time.time()-start,
                proof='For each canonical key, local partition, root numerator, all 113 scaled transport polynomials, 10 exterior triples and association pattern are preserved exactly; complete rational functions determined.')
    (HERE/'symbolic_enumeration_independent.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':run()
