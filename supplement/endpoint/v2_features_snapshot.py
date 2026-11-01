#!/usr/bin/env python3
"""Full V2 typed nonbacktracking walk LP, exact replay included."""
from __future__ import annotations
if not __debug__:
    raise RuntimeError("V3 exact checks refuse optimized Python (-O/-OO)")
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import gcd
from functools import reduce
from collections import Counter
import argparse,gzip,hashlib,json,time
import v1_universe_snapshot as v1
HERE=Path(__file__).resolve().parent
REV=[(t,a,b,c) for t in (0,1) for a in range(4) for b in range(4) for c in range(a+1,4)]
RID={p:i for i,p in enumerate(REV)}
NAMES=['normalization','m2','m3','q2','q3','j']+[f'R_{t}_{a}_{b}_{c}' for t,a,b,c in REV]

def dump(path,x):Path(path).write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')

def all_full_views():
    for rep,mult in v1.all_views():
        choices=[]
        for T,(x,b) in zip(rep['blocks'],rep['flags']):
            choices.append([k for k in (2,3) if k>=len(T) and (k!=len(T) or x==b)])
        assert mult==reduce(lambda a,b:a*b,map(len,choices),1)
        for degrees in product(*choices):
            yield dict(rep,boundary_degrees=degrees)

def numerator_features(rep):
    d=rep['d'];ds=rep['neighbor_degrees'];sets=rep['feasible_masks'];z=len(sets)
    occ=[sum(bool(S&(1<<u)) for S in sets) for u in range(d+1)]
    joint=[sum(bool(S&1) and bool(S&(1<<u)) for S in sets) for u in range(1,d+1)]
    nums=[6*occ[0],6*z*(d==2)-3*z*sum(k==2 for k in ds),6*z*(d==3)-2*z*sum(k==3 for k in ds),
          6*occ[0]*(d==2)-3*sum(occ[u] for u in range(1,d+1) if ds[u-1]==2),
          6*occ[0]*(d==3)-2*sum(occ[u] for u in range(1,d+1) if ds[u-1]==3),
          6*(sum(joint[u-1] for u in range(1,d+1) if d==2 and ds[u-1]==3)-sum(joint[u-1] for u in range(1,d+1) if d==3 and ds[u-1]==2))]
    paths=[] # middle, endpoint_inside_or_0, endpoint_degree, exterior_bit, adjacent
    for u,w in rep['neighbor_edges']:
        paths.append((u,w,ds[w-1],0,1));paths.append((w,u,ds[u-1],0,1))
    for T,k,(x,b) in zip(rep['blocks'],rep['boundary_degrees'],rep['flags']):
        for u in T:paths.append((u,0,k,x,0))
    assert len(paths)==sum(k-1 for k in ds)
    rev=[0]*48
    for S in sets:
        a=2*(d-2)+(S&1)
        for u,w,k,x,t in paths:
            b=2*(ds[u-1]-2)+int(bool(S&(1<<u)))
            c=2*(k-2)+(int(bool(S&(1<<w))) if w else x)
            if a<c:rev[RID[t,a,b,c]]+=6
            elif c<a:rev[RID[t,c,b,a]]-=6
    nums+=rev;den=6*z
    g=reduce(gcd,nums,den)
    return den//g,tuple(x//g for x in nums)

def enumerate_views():
    start=time.time();columns={};counts=Counter();sha=hashlib.sha256()
    for rep in all_full_views():
        key=numerator_features(rep)
        if key not in columns:columns[key]=dict(representative=rep,full_views=0)
        columns[key]['full_views']+=1
        counts['full_views']+=1;counts[f'degree_{rep["d"]}']+=1
        sha.update((json.dumps(rep,sort_keys=True,separators=(',',':'))+'\n').encode())
    rows=[]
    for (den,nums),info in sorted(columns.items()):rows.append(dict(id=len(rows),denominator=den,numerators=nums,**info))
    data=dict(schema='V2 common denominator; c,m2,m3,q2,q3,j then48 reversals',
              symbol_order=[[2,0],[2,1],[3,0],[3,1]],reversal_order=REV,
              counts=dict(counts),columns=len(rows),generation_sha256=sha.hexdigest(),rows=rows)
    with gzip.open(HERE/'views_v2.json.gz','wt') as f:json.dump(data,f,sort_keys=True)
    out={k:v for k,v in data.items() if k!='rows'};out['elapsed_seconds']=time.time()-start
    dump(HERE/'enumeration_v2.json',out);print(json.dumps(out,indent=2))

def load():
    with gzip.open(HERE/'views_v2.json.gz','rt') as f:return json.load(f)

def feats(r):return tuple(F(n,r['denominator']) for n in r['numerators'])

def solve():
    gate=HERE/'audit'/'validation_v2.json'
    if not gate.exists() or json.loads(gate.read_text()).get('all_passed') is not True:raise RuntimeError('V2 independent validation gate missing/not passed')
    import numpy as np
    import sympy as sp
    from scipy.optimize import linprog
    data=load();rows=data['rows'];C=[feats(r) for r in rows];m=len(C[0])
    A=np.array([[1.]*len(C)]+[[float(f[k]) for f in C] for k in range(1,m)])
    c=np.array([float(f[0]) for f in C]);b=np.array([1.]+[0.]*(m-1))
    ans=linprog(c,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
    raw=dict(success=bool(ans.success),status=int(ans.status),message=ans.message,objective=float(ans.fun) if ans.success else None)
    if ans.success:
        support=[i for i,x in enumerate(ans.x) if x>1e-9]
        raw.update(support=[dict(id=i,weight=float(ans.x[i])) for i in support],
                   equality_marginals=ans.eqlin.marginals.tolist(),equality_residual=ans.eqlin.residual.tolist())
    dump(HERE/'lp_v2_numerical.json',raw);print(json.dumps(raw,indent=2))
    if not ans.success:return
    matrix=sp.Matrix([[1 for i in support]]+[[sp.Rational(C[i][k].numerator,C[i][k].denominator) for i in support] for k in range(1,m)])
    weights=matrix.gauss_jordan_solve(sp.Matrix([1]+[0]*(m-1)))[0]
    assert not any(x.free_symbols for x in weights)
    weights=[F(str(x)) for x in weights]
    dual=[F(float(x)).limit_denominator(10**9) for x in ans.eqlin.marginals]
    cert=dict(schema='V2 exact rational primal/dual candidate',generation_sha256=data['generation_sha256'],
              target='4/11',support=[dict(id=i,weight=str(w),features=list(map(str,C[i])),representative=rows[i]['representative']) for i,w in zip(support,weights)],
              objective=str(sum((w*C[i][0] for i,w in zip(support,weights)),F(0))),
              dual=dict(zip(NAMES,map(str,dual))))
    dump(HERE/'certificate_v2_attempt1.json',cert)
    # Recover from ALL tight columns, not only the sparse primal support.
    import runpy
    runpy.run_path(str(HERE/'recover_dual_v2.py'),run_name='__main__')

def verify(regenerate=True):
    data=load();rows=data['rows'];cert=json.loads((HERE/'certificate_v2.json').read_text())
    assert cert['generation_sha256']==data['generation_sha256'];m=len(NAMES)
    sums=[F(0)]*m;obj=F(0)
    for s in cert['support']:
        r=rows[s['id']];f=feats(r);w=F(s['weight']);assert w>=0
        assert list(map(str,f))==s['features'];assert r['representative']==s['representative']
        sums[0]+=w
        for j in range(1,m):sums[j]+=w*f[j]
        obj+=w*f[0]
    assert sums==[F(1)]+[F(0)]*(m-1)
    dual=[F(cert['dual'][n]) for n in NAMES];assert obj==F(cert['objective'])==dual[0]
    tight=0;minimum=None
    for r in rows:
        f=feats(r);slack=f[0]-dual[0]-sum((dual[j]*f[j] for j in range(1,m)),F(0))
        assert slack>=0,(r['id'],slack)
        tight+=slack==0;minimum=slack if minimum is None else min(slack,minimum)
    out=dict(all_passed=True,objective=str(obj),target_gap=str(obj-F(4,11)),
             support_size=len(cert['support']),primal_sums=list(map(str,sums)),
             columns=len(rows),tight_columns=tight,minimum_dual_slack=str(minimum),regenerated=False)
    if regenerate:
        sha=hashlib.sha256();count=0;known={(r['denominator'],tuple(r['numerators'])) for r in rows}
        # Integer-scaled inequality avoids repeatedly constructing all54 Fractions.
        D=1
        for x in dual:D=D*x.denominator//gcd(D,x.denominator)
        yd=[x.numerator*(D//x.denominator) for x in dual]
        for rep in all_full_views():
            den,nums=numerator_features(rep);assert(den,nums) in known
            assert D*nums[0]-yd[0]*den-sum(yd[j]*nums[j] for j in range(1,m))>=0
            sha.update((json.dumps(rep,sort_keys=True,separators=(',',':'))+'\n').encode());count+=1
        assert count==data['counts']['full_views'];assert sha.hexdigest()==data['generation_sha256']
        out.update(regenerated=True,full_views_checked=count,generation_sha256=sha.hexdigest())
    dump(HERE/('verification_v2.json' if regenerate else 'verification_v2_columns.json'),out)
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['enumerate','solve','verify']);a=p.parse_args()
    {'enumerate':enumerate_views,'solve':solve,'verify':verify}[a.command]()
