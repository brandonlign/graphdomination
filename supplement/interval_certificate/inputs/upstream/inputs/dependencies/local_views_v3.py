#!/usr/bin/env python3
"""V3 enumeration, numerical search and exact primal/dual replay."""
from __future__ import annotations
if not __debug__:
    raise RuntimeError('V3 exact checks refuse optimized Python (-O/-OO)')
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
from functools import reduce
from math import gcd,lcm
from collections import Counter
import argparse,gzip,hashlib,json,time
import v2_features_snapshot as v2
HERE=Path(__file__).resolve().parent
STARS=[(d,x,T) for d in (2,3) for x in (0,1) for T in combinations_with_replacement(range(4),d)]
SID={s:i for i,s in enumerate(STARS)}
CLASSES=[(d,h,k)for d in(2,3)for h in(2,3)for k in range(1,h+1)]
PATTERNS=[(2,r)for r in range(3)]+[(3,r)for r in range(4)]
BOUNDS=list(map(F,['1/3','9/25','5/13','79/207','69/197','15/47','13/45']))
EQNAMES=v2.NAMES+[f'S_{d}_{x}_'+''.join(map(str,T)) for d,x,T in STARS]
GNAMES=[f'{typ}_{d}_{h}_{k}'for d,h,k in CLASSES for typ in('X','B')]+[f'L_{d}_{r}'for d,r in PATTERNS]
NEQ=len(EQNAMES);NG=len(GNAMES)
assert NEQ==114 and NG==27

def dump(path,obj):Path(path).write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')

def feature_key(rep):
    oldden,old=v2.numerator_features(rep)
    d=rep['d'];ds=rep['neighbor_degrees'];sets=rep['feasible_masks'];z=len(sets)
    r=sum(k==3 for k in ds)if d==2 else sum(k==2 for k in ds)
    pattern=PATTERNS.index((d,r));L=BOUNDS[pattern]
    D=lcm(oldden,6*z,4,L.denominator)
    nums=[a*(D//oldden)for a in old]
    star=[0]*60
    # Neighbor-star symbol templates consist of root, inside neighbors, fixed boundary.
    inside={u:[] for u in range(1,d+1)};outside={u:[]for u in range(1,d+1)}
    for a,b in rep['neighbor_edges']:inside[a].append(b);inside[b].append(a)
    for T,h,(x,b) in zip(rep['blocks'],rep['boundary_degrees'],rep['flags']):
        for u in T:outside[u].append(2*(h-2)+x)
    for S in sets:
        bits=[int(bool(S&(1<<u)))for u in range(d+1)]
        symbols=[2*(ds[u-1]-2)+bits[u]for u in range(1,d+1)]
        star[SID[(d,bits[0],tuple(sorted(symbols)))]]+=6
        for u in range(1,d+1):
            neighbors=[2*(d-2)+bits[0]]+[symbols[w-1]for w in inside[u]]+outside[u]
            assert len(neighbors)==ds[u-1]
            star[SID[(ds[u-1],bits[u],tuple(sorted(neighbors)))]]-=6//ds[u-1]
    nums += [a*(D//(6*z))for a in star]
    assert len(nums)==NEQ
    for rd,h,k in CLASSES:
        chosen=[(x,b) for T,hd,(x,b)in zip(rep['blocks'],rep['boundary_degrees'],rep['flags'])if d==rd and hd==h and len(T)==k]
        N=len(chosen);X=sum(x for x,b in chosen);B=sum(b for x,b in chosen)
        nums += [N*(D//2)-X*D,X*(D//(2**(h-k)))-B*D]
    rootcount=sum(bool(S&1)for S in sets)
    for p,Lp in enumerate(BOUNDS):nums.append(rootcount*(D//z)-Lp.numerator*(D//Lp.denominator)if p==pattern else 0)
    assert len(nums)==NEQ+NG
    g=reduce(gcd,nums,D)
    return d,D//g,tuple(n//g for n in nums)

def enumerate_views():
    start=time.time();columns={};counts=Counter();sha=hashlib.sha256()
    for rep in v2.all_full_views():
        key=feature_key(rep)
        if key not in columns:columns[key]=dict(representative=rep,full_views=0)
        columns[key]['full_views']+=1;counts['full_views']+=1;counts[f'degree_{rep["d"]}']+=1
        sha.update((json.dumps(rep,sort_keys=True,separators=(',',':'))+'\n').encode())
    rows=[dict(id=i,d=d,denominator=den,numerators=nums,**info)for i,((d,den,nums),info)in enumerate(sorted(columns.items()))]
    out=dict(schema='V3: c +113 zero features +27 nonnegative-mean features; root d retained for secondary',
             equality_names=EQNAMES,inequality_names=GNAMES,star_order=STARS,classes=CLASSES,patterns=PATTERNS,bounds=list(map(str,BOUNDS)),
             columns=len(rows),counts=dict(counts),generation_sha256=sha.hexdigest(),rows=rows)
    with gzip.open(HERE/'views_v3.json.gz','wt')as f:json.dump(out,f,sort_keys=True)
    summary={k:v for k,v in out.items()if k!='rows'};summary['elapsed_seconds']=time.time()-start
    dump(HERE/'enumeration_v3.json',summary);print(json.dumps({k:v for k,v in summary.items()if k in('columns','counts','generation_sha256','elapsed_seconds')},indent=2))

def load():
    with gzip.open(HERE/'views_v3.json.gz','rt')as f:return json.load(f)

def feats(row):return tuple(F(x,row['denominator'])for x in row['numerators'])

def matrices(rows,secondary=False):
    C=[feats(r)for r in rows]
    A=[[F(1)]*len(rows)]+[[f[j]for f in C]for j in range(1,NEQ)]
    b=[F(1)]+[F(0)]*(NEQ-1)
    if secondary:A.append([f[0]for f in C]);b.append(F(4,11))
    G=[[f[NEQ+j]for f in C]for j in range(NG)]
    obj=[-F(r['d']==3)for r in rows]if secondary else[f[0]for f in C]
    return C,A,b,G,obj

def evaluate_certificate(cert,data,regenerate=False):
    assert cert['problem'] in ('primary','secondary')
    assert cert['generation_sha256']==data['generation_sha256']
    assert data['equality_names']==EQNAMES and data['inequality_names']==GNAMES
    rows=data['rows'];secondary=cert['problem']=='secondary';C,A,b,G,obj=matrices(rows,secondary)
    y=list(map(F,cert['equality_dual']));mu=list(map(F,cert['inequality_dual']))
    assert len(y)==len(A)and len(mu)==NG and all(x>=0 for x in mu)
    eq=[F(0)]*len(A);ineq=[F(0)]*NG;value=F(0)
    for item in cert['support']:
        i=item['id'];w=F(item['weight']);assert w>=0
        assert item['representative']==rows[i]['representative'];assert item['features']==list(map(str,C[i]))
        for j in range(len(A)):eq[j]+=w*A[j][i]
        for j in range(NG):ineq[j]+=w*G[j][i]
        value+=w*obj[i]
    assert eq==b and all(x>=0 for x in ineq)
    dual=sum((u*v for u,v in zip(y,b)),F(0))
    assert value==F(cert['objective'])==dual
    tight=0;mint=None
    for i in range(len(rows)):
        slack=obj[i]-sum((y[j]*A[j][i]for j in range(len(A))),F(0))-sum((mu[j]*G[j][i]for j in range(NG)),F(0))
        assert slack>=0,(i,slack)
        tight+=slack==0;mint=slack if mint is None else min(mint,slack)
    out=dict(all_passed=True,problem=cert['problem'],objective=str(value),dual_objective=str(dual),support_size=len(cert['support']),
             equality_sums=list(map(str,eq)),inequality_sums=list(map(str,ineq)),inequality_dual_nonnegative=True,
             columns=len(rows),tight_columns=tight,minimum_slack=str(mint),regenerated=False)
    if regenerate:
        H=1
        for x in y+mu:H=lcm(H,x.denominator)
        yi=[int(x*H)for x in y];mi=[int(x*H)for x in mu]
        known={(r['d'],r['denominator'],tuple(r['numerators']))for r in rows}
        count=0;sha=hashlib.sha256()
        for rep in v2.all_full_views():
            d,den,n=feature_key(rep);assert(d,den,n)in known
            target=-den*int(d==3)if secondary else n[0]
            rhs=yi[0]*den+sum(yi[j]*n[j]for j in range(1,NEQ))+sum(mi[j]*n[NEQ+j]for j in range(NG))
            if secondary:rhs+=yi[-1]*n[0]
            assert H*target>=rhs
            sha.update((json.dumps(rep,sort_keys=True,separators=(',',':'))+'\n').encode());count+=1
        assert count==data['counts']['full_views']and sha.hexdigest()==data['generation_sha256']
        out.update(regenerated=True,full_views_checked=count,generation_sha256=sha.hexdigest())
    return out

def solve(secondary=False):
    import verify_association_source
    verify_association_source.verify()
    gate=HERE/'audit'/'validation_v3.json'
    assert gate.exists()and json.loads(gate.read_text()).get('all_passed')is True,'V3 validation gate missing/not passed'
    if secondary:
        first=json.loads((HERE/'certificate_primary_v3.json').read_text());evaluate_certificate(first,load())
        assert F(first['objective'])==F(4,11),'secondary only authorized at primary optimum4/11'
    import numpy as np
    import sympy as sp
    from scipy.optimize import linprog
    data=load();rows=data['rows'];C,A,b,G,obj=matrices(rows,secondary);kind='secondary'if secondary else'primary'
    ans=linprog(np.array(list(map(float,obj))),A_eq=np.array([[float(x)for x in row]for row in A]),b_eq=np.array(list(map(float,b))),
                A_ub=-np.array([[float(x)for x in row]for row in G]),b_ub=np.zeros(NG),bounds=(0,None),method='highs')
    raw=dict(success=bool(ans.success),status=int(ans.status),message=ans.message,objective=float(ans.fun)if ans.success else None)
    if ans.success:
        support=[i for i,w in enumerate(ans.x)if w>1e-9]
        raw.update(support=[dict(id=i,weight=float(ans.x[i]))for i in support],equality_dual=ans.eqlin.marginals.tolist(),
                   inequality_dual=(-ans.ineqlin.marginals).tolist(),equality_residual=ans.eqlin.residual.tolist(),inequality_slack=ans.ineqlin.residual.tolist())
    dump(HERE/f'lp_{kind}_v3_numerical.json',raw);print(json.dumps({k:v for k,v in raw.items()if k in('success','message','objective','support')},indent=2))
    if not ans.success:return
    rational=lambda x:sp.Rational(x.numerator,x.denominator)
    binding=[j for j,x in enumerate(ans.ineqlin.residual)if abs(x)<1e-8]
    M=sp.Matrix([[rational(row[i])for i in support]for row in A]+[[rational(G[j][i])for i in support]for j in binding])
    rhs=sp.Matrix([rational(x)for x in b]+[0]*len(binding))
    weights,params=M.gauss_jordan_solve(rhs);assert params.rows==0,'nonunique primal recovery'
    weights=[F(str(x))for x in weights]
    yfloat=raw['equality_dual'];mufloat=raw['inequality_dual'];qfloat=yfloat+mufloat
    active=[j for j,x in enumerate(qfloat)if abs(x)>1e-9]
    combined=A+G
    tight=[i for i in range(len(rows))if abs(float(obj[i])-sum(qfloat[j]*float(combined[j][i])for j in range(len(qfloat))))<1e-8]
    D=sp.Matrix([[rational(combined[j][i])for j in active]for i in tight]);db=sp.Matrix([rational(obj[i])for i in tight])
    sol,free=D.gauss_jordan_solve(db)
    base=dict(problem=kind,generation_sha256=data['generation_sha256'],objective=str(sum((w*obj[i]for i,w in zip(support,weights)),F(0))),
              support=[dict(id=i,weight=str(w),features=list(map(str,C[i])),representative=rows[i]['representative'])for i,w in zip(support,weights)])
    attempts=[]
    for lim in (10**6,10**9,10**12):
        subs={}
        for par in free:
            loc=next(j for j,s in enumerate(sol)if s==par)
            frac=F(qfloat[active[loc]]).limit_denominator(lim);subs[par]=rational(frac)
        exact=[F(0)]*len(qfloat)
        for j,s in zip(active,sol):exact[j]=F(str(s.subs(subs)))
        cert=dict(base,equality_dual=list(map(str,exact[:len(A)])),inequality_dual=list(map(str,exact[len(A):])))
        path=HERE/f'certificate_{kind}_v3_attempt_{lim}.json';dump(path,cert)
        try:out=evaluate_certificate(cert,data)
        except Exception as e:
            attempts.append(dict(limit=lim,passed=False,error=repr(e)));continue
        attempts.append(dict(limit=lim,passed=True))
        dump(HERE/f'certificate_{kind}_v3.json',cert);dump(HERE/f'verification_{kind}_v3_columns.json',out)
        recovery=dict(problem=kind,tight_columns_used=len(tight),active_dual_coordinates=len(active),exact_rank=len(active)-free.rows,free_parameters=free.rows,attempts=attempts)
        dump(HERE/f'recovery_{kind}_v3.json',recovery)
        print(json.dumps(dict(exact_objective=cert['objective'],support_size=len(support),dual_nonzero=len(active),tight_columns=out['tight_columns'],free_parameters=free.rows),indent=2))
        return
    dump(HERE/f'recovery_{kind}_v3_failed.json',dict(attempts=attempts));raise RuntimeError('All exact dual recoveries failed; candidates preserved')

def verify(secondary=False):
    import verify_association_source
    verify_association_source.verify()
    kind='secondary'if secondary else'primary';cert=json.loads((HERE/f'certificate_{kind}_v3.json').read_text())
    out=evaluate_certificate(cert,load(),regenerate=True);dump(HERE/f'verification_{kind}_v3.json',out)
    print(json.dumps({k:v for k,v in out.items()if k not in('equality_sums','inequality_sums')},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['enumerate','solve','verify','secondary','verify-secondary']);a=p.parse_args()
    {'enumerate':enumerate_views,'solve':solve,'verify':verify,'secondary':lambda:solve(True),'verify-secondary':lambda:verify(True)}[a.command]()
