#!/usr/bin/env python3
"""V1 local-view enumeration and rational LP certificate replay.

Only the `solve` subcommand uses scipy. Exact modes use Python's Fraction.
The full degree expansion is represented by exact positive multiplicities.
"""
from __future__ import annotations
if not __debug__:
    raise RuntimeError("V3 exact checks refuse optimized Python (-O/-OO)")
import argparse
import gzip
import hashlib
import itertools
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
NAMES=['normalization','m2','m3','q2','q3','j']

def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def partitions(stubs):
    """Restricted-growth partition generator, first-appearance block order."""
    blocks=[]
    def rec(t):
        if t==len(stubs):
            yield tuple(tuple(b) for b in blocks)
            return
        u=stubs[t]
        for b in blocks:
            if len(b)<3 and all(stubs[s]!=u for s in b):
                b.append(t)
                yield from rec(t+1)
                b.pop()
        blocks.append([t])
        yield from rec(t+1)
        blocks.pop()
    yield from rec(0)

def flag_options(block):
    """(x,b), multiplicity of admissible boundary degrees, representative degree."""
    s=len(block)
    if s==1:
        return [((0,0),2,2),((1,0),2,2),((1,1),2,2)]
    if s==2:
        return [((0,0),2,2),((1,0),1,3),((1,1),2,2)]
    assert s==3
    return [((0,0),1,3),((1,1),1,3)]

def feasible_masks(d,edges,blocks,flags):
    forbidden=[(1<<(d+1))-1]
    for u in range(1,d+1):
        if all(flags[t][0] for t,T in enumerate(blocks) if u in T):
            mask=1|(1<<u)
            for a,b in edges:
                if a==u:mask |= 1<<b
                if b==u:mask |= 1<<a
            forbidden.append(mask)
    for T,(x,b) in zip(blocks,flags):
        if b:forbidden.append(sum(1<<u for u in T))
    return tuple(S for S in range(1<<(d+1)) if all(S&Q!=Q for Q in forbidden))

def features(d,degrees,masks):
    z=len(masks)
    p=[F(sum(bool(S&(1<<u)) for S in masks),z) for u in range(d+1)]
    joint=[F(sum(bool(S&1) and bool(S&(1<<u)) for S in masks),z) for u in range(1,d+1)]
    m=[F(d==k)-F(sum(g==k for g in degrees),k) for k in (2,3)]
    q=[(p[0] if d==k else 0)-sum((p[u]/k for u in range(1,d+1) if degrees[u-1]==k),F(0)) for k in (2,3)]
    j=sum((joint[u-1] for u in range(1,d+1) if d==2 and degrees[u-1]==3),F(0))-sum((joint[u-1] for u in range(1,d+1) if d==3 and degrees[u-1]==2),F(0))
    return (p[0],*m,*q,j)

def all_views():
    for d in (2,3):
        pairs=list(itertools.combinations(range(1,d+1),2))
        for emask in range(1<<len(pairs)):
            edges=tuple(p for j,p in enumerate(pairs) if emask>>j&1)
            internal=[sum(u in e for e in edges) for u in range(1,d+1)]
            for degrees in itertools.product((2,3),repeat=d):
                ns=[degrees[u]-1-internal[u] for u in range(d)]
                if min(ns)<0:continue
                stubs=tuple(u+1 for u,n in enumerate(ns) for _ in range(n))
                for parts in partitions(stubs):
                    blocks=tuple(tuple(stubs[s] for s in P) for P in parts)
                    opts=[flag_options(T) for T in blocks]
                    for choices in itertools.product(*opts):
                        flags=tuple(c[0] for c in choices)
                        mult=1
                        for c in choices:mult*=c[1]
                        bdegrees=tuple(c[2] for c in choices)
                        masks=feasible_masks(d,edges,blocks,flags)
                        rep=dict(d=d,neighbor_degrees=degrees,neighbor_edges=edges,
                                 stub_partition=parts,blocks=blocks,boundary_degrees=bdegrees,
                                 flags=flags,feasible_masks=masks)
                        yield rep,mult

def enumerate_views():
    start=time.time(); families={};columns={};counts=Counter();sha=hashlib.sha256()
    for rep,mult in all_views():
        d=rep['d'];degrees=tuple(rep['neighbor_degrees']);masks=tuple(rep['feasible_masks'])
        fkey=(d,degrees,masks)
        if fkey not in families:
            feat=features(d,degrees,masks)
            families[fkey]=feat
            if feat not in columns:
                columns[feat]=dict(representative=rep,coarse_views=0,full_views=0,families=0)
            columns[feat]['families']+=1
        feat=families[fkey]
        columns[feat]['coarse_views']+=1
        columns[feat]['full_views']+=mult
        counts['coarse_views']+=1;counts['full_views']+=mult
        counts[f'coarse_views_degree_{d}']+=1;counts[f'full_views_degree_{d}']+=mult
        sha.update((json.dumps([rep,mult],sort_keys=True,separators=(',',':'))+'\n').encode())
    rows=[]
    for feat,item in sorted(columns.items()):
        rows.append(dict(id=len(rows),features=[str(x) for x in feat],**item))
    data=dict(schema='V1: c,m2,m3,q2,q3,j; full degrees compressed by multiplicity',
              counts=dict(counts),families=len(families),columns=len(rows),
              generation_sha256=sha.hexdigest(),rows=rows)
    with gzip.open(HERE/'views_v1.json.gz','wt') as f:json.dump(data,f,sort_keys=True)
    summary={k:v for k,v in data.items() if k!='rows'}
    summary['elapsed_seconds']=time.time()-start
    dump(HERE/'enumeration_v1.json',summary)
    print(json.dumps(summary,indent=2))

def load_views():
    with gzip.open(HERE/'views_v1.json.gz','rt') as f:return json.load(f)

def solve():
    # This gate prevents use before independent all-exposure validation.
    gate=HERE/'audit'/'validation_results.json'
    if not gate.exists():raise RuntimeError('Independent validation gate absent: '+str(gate))
    validation=json.loads(gate.read_text())
    if validation.get('all_passed') is not True:raise RuntimeError('Validation gate did not pass')
    import numpy as np
    from scipy.optimize import linprog
    import sympy as sp
    data=load_views();rows=data['rows']
    C=[tuple(F(x) for x in r['features']) for r in rows]
    A=np.array([[1.0]*len(C)]+[[float(f[k]) for f in C] for k in range(1,6)])
    c=np.array([float(f[0]) for f in C]);b=np.array([1.,0,0,0,0,0])
    result=linprog(c,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
    raw=dict(success=bool(result.success),status=int(result.status),message=result.message,
             objective=float(result.fun) if result.success else None)
    if result.success:
        support=[i for i,x in enumerate(result.x) if x>1e-9]
        raw.update(support=[dict(id=i,weight=float(result.x[i])) for i in support],
                   equality_marginals=result.eqlin.marginals.tolist(),
                   equality_residual=result.eqlin.residual.tolist())
    dump(HERE/'lp_v1_numerical.json',raw)
    print(json.dumps(raw,indent=2))
    if not result.success:return
    matrix=sp.Matrix([[1 for _ in support]]+[[sp.Rational(C[i][k].numerator,C[i][k].denominator) for i in support] for k in range(1,6)])
    target=sp.Matrix([1,0,0,0,0,0]);weights=matrix.gauss_jordan_solve(target)[0]
    if any(x.free_symbols for x in weights):raise RuntimeError('Primal support not a unique basic solution')
    exact_weights=[F(str(w)) for w in weights]
    dual=[F(float(x)).limit_denominator(10**8) for x in result.eqlin.marginals]
    cert=dict(schema='V1 exact primal and dual',generation_sha256=data['generation_sha256'],
              target='4/11',support=[dict(id=i,weight=str(w),features=rows[i]['features'],representative=rows[i]['representative']) for i,w in zip(support,exact_weights)],
              objective=str(sum((w*C[i][0] for i,w in zip(support,exact_weights)),F(0))),
              dual=dict(zip(NAMES,map(str,dual))))
    dump(HERE/'certificate_v1.json',cert)
    verify(regenerate=False)

def verify(regenerate=True):
    data=load_views();cert=json.loads((HERE/'certificate_v1.json').read_text());rows=data['rows']
    assert cert['generation_sha256']==data['generation_sha256']
    sums=[F(0)]*6;obj=F(0)
    for item in cert['support']:
        r=rows[item['id']];assert item['features']==r['features'];assert item['representative']==r['representative']
        w=F(item['weight']);assert w>=0
        f=list(map(F,r['features']));sums[0]+=w
        for k in range(1,6):sums[k]+=w*f[k]
        obj+=w*f[0]
    assert sums==[F(1),F(0),F(0),F(0),F(0),F(0)]
    assert obj==F(cert['objective'])
    dual=[F(cert['dual'][n]) for n in NAMES]
    min_slack=None;tight=0
    for r in rows:
        f=list(map(F,r['features']));slack=f[0]-dual[0]-sum((dual[k]*f[k] for k in range(1,6)),F(0))
        assert slack>=0,(r['id'],slack)
        if slack==0:tight+=1
        if min_slack is None or slack<min_slack:min_slack=slack
    assert obj==dual[0]
    out=dict(all_passed=True,primal_sums=list(map(str,sums)),objective=str(obj),
             target_gap=str(obj-F(4,11)),dual_objective=str(dual[0]),
             unique_columns_checked=len(rows),tight_unique_columns=tight,
             min_dual_slack=str(min_slack),regenerated=False)
    if regenerate:
        sha=hashlib.sha256();count=0;full=0;family_cache={};known={tuple(F(x) for x in r['features']) for r in rows}
        for rep,mult in all_views():
            key=(rep['d'],tuple(rep['neighbor_degrees']),tuple(rep['feasible_masks']))
            if key not in family_cache:family_cache[key]=features(*key)
            f=family_cache[key]
            assert f in known
            assert f[0]-dual[0]-sum((dual[k]*f[k] for k in range(1,6)),F(0))>=0
            sha.update((json.dumps([rep,mult],sort_keys=True,separators=(',',':'))+'\n').encode())
            count+=1;full+=mult
        assert sha.hexdigest()==data['generation_sha256']
        assert count==data['counts']['coarse_views'];assert full==data['counts']['full_views']
        out.update(regenerated=True,coarse_views_checked=count,full_views_covered=full,
                   generation_sha256=sha.hexdigest())
    dump(HERE/('verification_v1.json' if regenerate else 'verification_v1_columns.json'),out)
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['enumerate','solve','verify']);a=p.parse_args()
    {'enumerate':enumerate_views,'solve':solve,'verify':verify}[a.command]()
