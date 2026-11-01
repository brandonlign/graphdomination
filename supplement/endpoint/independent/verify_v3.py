#!/usr/bin/env python3
"""Independent exact V3 primary/optimum-face certificate replay."""
import sys
if not __debug__: raise RuntimeError('Optimized Python is forbidden for V3 exact checks')
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import gzip,json,argparse
from features_v3_independent import feature_vector,canonical,SIGS,CLASSES,PATTERNS
from enumerate_v2_independent import COORDS,family
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
EQNAMES=['normalization','m2','m3','q2','q3','j']+[f'R_{t}_{a}_{b}_{c}' for t,a,b,c in COORDS]+[f'S_{d}_{x}_'+''.join(map(str,T)) for d,x,T in SIGS]
GNAMES=[f'{typ}_{d}_{h}_{k}' for d,h,k in CLASSES for typ in ('X','B')]+[f'L_{d}_{r}' for d,r,L in PATTERNS]

def load(p):
    with gzip.open(p,'rt') as f:return json.load(f)

def key(r):return tuple(r['numerators'])+(r['denominator'],)

def repkey(R):
    d=R['d'];deg=tuple(R['neighbor_degrees']);E=frozenset(tuple(e) for e in R['neighbor_edges'])
    assert d in (2,3) and len(deg)==d and all(h in (2,3) for h in deg)
    assert all(1<=u<v<=d for u,v in E) and len(E)==len(R['neighbor_edges'])
    blocks=tuple(sorted((tuple(T),h,x,b) for T,h,(x,b) in zip(R['blocks'],R['boundary_degrees'],R['flags'])))
    assert len(blocks)==len(R['blocks'])==len(R['boundary_degrees'])==len(R['flags'])
    for T,h,x,b in blocks:
        assert T and len(T)==len(set(T)) and all(1<=u<=d for u in T)
        assert h in (2,3) and h>=len(T) and (x,b) in ((0,0),(1,0),(1,1))
        if h==len(T):assert x==b
    for u in range(1,d+1):assert deg[u-1]==1+sum(u in e for e in E)+sum(u in T for T,h,x,b in blocks)
    masks=sorted(sum(1<<u for u in S) for S in family(d,E,blocks))
    assert masks==sorted(R['feasible_masks']) and len(masks)==len(R['feasible_masks'])
    return canonical(feature_vector(d,deg,E,blocks)),int(d==3)

def run(kind):
    independent=load(HERE/'independent_columns_v3.json.gz');primary=load(ROOT/'views_v3.json.gz')
    assert primary['equality_names']==EQNAMES and primary['inequality_names']==GNAMES
    assert primary['star_order']==[[d,x,list(T)] for d,x,T in SIGS]
    assert primary['classes']==[list(C) for C in CLASSES]
    assert primary['patterns']==[[d,r] for d,r,L in PATTERNS]
    assert primary['bounds']==[str(L) for d,r,L in PATTERNS]
    assert all(r['d'] in (2,3) and r['id']==i for i,r in enumerate(primary['rows']))
    left={(key(r),r['degree3']):r['full_views'] for r in independent['rows']}
    right={(key(r),int(r['d']==3)):r['full_views'] for r in primary['rows']}
    assert len(primary['rows'])==len(right)==len(independent['rows'])==len(left)==9853
    assert left==right and sum(left.values())==231816
    path=ROOT/f'certificate_{kind}_v3.json';cert=json.loads(path.read_text())
    assert cert['problem']==kind and cert['generation_sha256']==primary['generation_sha256']
    secondary=kind=='secondary';y=list(map(F,cert['equality_dual']));mu=list(map(F,cert['inequality_dual']))
    assert len(y)==114+int(secondary) and len(mu)==27 and all(x>=0 for x in mu)
    H=lcm(*(q.denominator for q in y+mu));yi=[int(q*H) for q in y];mi=[int(q*H) for q in mu]
    tightcols=0;tightfull=0;minslack=None
    for (K,d3),multiplicity in left.items():
        nums=K[:-1];D=K[-1];objective=-d3*D if secondary else nums[0]
        rhs=yi[0]*D+sum(yi[j]*nums[j] for j in range(1,114))+sum(mi[j]*nums[114+j] for j in range(27))
        if secondary:rhs+=yi[-1]*nums[0]
        slack=H*objective-rhs;assert slack>=0,(K,d3,slack)
        tightcols+=int(slack==0);tightfull+=multiplicity*int(slack==0)
        q=F(slack,H*D);minslack=q if minslack is None else min(minslack,q)
    eq=[F(0)]*len(y);g=[F(0)]*27;objective=F(0);root3=F(0);support=[]
    for S in cert['support']:
        i=S['id'];K,d3=repkey(S['representative']);assert (K,d3) in left
        row=primary['rows'][i];assert K==key(row) and d3==int(row['d']==3)
        f=tuple(F(n,K[-1]) for n in K[:-1]);assert list(map(str,f))==S['features']
        w=F(S['weight']);assert w>=0;eq[0]+=w
        for j in range(1,114):eq[j]+=w*f[j]
        if secondary:eq[-1]+=w*f[0]
        for j in range(27):g[j]+=w*f[114+j]
        objective+=w*(-d3 if secondary else f[0]);root3+=w*d3
        support.append({'id':i,'weight':str(w),'root_degree':3 if d3 else 2,'occupancy':str(f[0])})
    assert eq==[F(1)]+[F(0)]*113+([F(4,11)] if secondary else []) and all(x>=0 for x in g)
    dual=y[0]+(F(4,11)*y[-1] if secondary else F(0))
    assert objective==F(cert['objective'])==dual
    if secondary:assert objective==-root3
    out={'all_passed':True,'problem':kind,'independent_columns_match':True,'full_multiplicities_match':True,
         'full_views_covered':sum(left.values()),'columns_checked':len(left),'inequality_multipliers_nonnegative':True,
         'exact_objective':str(objective),'dual_objective':str(dual),'primary_occupancy':str(eq[-1] if secondary else objective),
         'root_degree3_mass':str(root3),'support':support,'equality_sums':list(map(str,eq)),
         'inequality_sums':list(map(str,g)),'tight_columns':tightcols,'tight_full_views':tightfull,
         'minimum_dual_slack':str(minslack),'generation_sha256':primary['generation_sha256']}
    out['nonzero_dual_by_family']={'normalization':int(y[0]!=0),'V1':sum(q!=0 for q in y[1:6]),'reversals':sum(q!=0 for q in y[6:54]),'stars':sum(q!=0 for q in y[54:114]),'exterior':sum(q!=0 for q in mu[:20]),'association':sum(q!=0 for q in mu[20:]),'face':int(secondary and y[-1]!=0)}
    if secondary:
        assert y[-1]<0 and y[0]+F(4,11)*y[-1]==0
        out['face_multiplier']=str(y[-1]);out['degree3_gap_coefficient']=str(-1/y[-1])
    (HERE/f'certificate_audit_{kind}_v3.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('equality_sums','inequality_sums')},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('kind',choices=('primary','secondary'));a=p.parse_args();run(a.kind)
