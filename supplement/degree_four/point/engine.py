#!/usr/bin/env python3
"""Exact z=1 full degree-four evaluator. Replay needs standard library only."""
from __future__ import annotations
if not __debug__: raise RuntimeError('Exact checks refuse Python -O/-OO')
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
from functools import reduce
from math import gcd,lcm,prod
from collections import Counter,defaultdict
import argparse,gzip,hashlib,importlib.util,json,sys,time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
GEN_SHA='af0020c5efce7b021f901e798013b34aa7967c019c4e5c98b8352e0e2985f021'
class ExactCheckFailure(RuntimeError): pass
spec=importlib.util.spec_from_file_location('view_grammar',HERE/'inputs/enumerate_validate.py')
grammar=importlib.util.module_from_spec(spec);spec.loader.exec_module(grammar)
def require(ok,msg):
    if not ok: raise ExactCheckFailure(msg)
def dump(path,x): Path(path).write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
def canonical(x): return json.dumps(x,sort_keys=True,separators=(',',':'))
def low(symbol): return symbol<2
def edgeok(a,b): return low(a) or low(b)
REV=tuple((t,a,b,c) for t in (0,1) for a in range(6) for b in range(6) for c in range(a+1,6)
          if edgeok(a,b) and edgeok(b,c) and (not t or edgeok(a,c)))
RID={s:i for i,s in enumerate(REV)}
STARS=tuple((d,x,T) for d in (2,3,4) for x in (0,1)
            for T in combinations_with_replacement(range(6 if d==2 else 2),d))
SID={s:i for i,s in enumerate(STARS)}
CLASSES=tuple((d,h,k) for d in (2,3,4) for h in (2,3,4) for k in range(1,min(d,h)+1))
CID={s:i for i,s in enumerate(CLASSES)}
PATTERNS=tuple((d,T) for d in (2,3,4) for T in combinations_with_replacement((2,3,4) if d==2 else (2,),d))
BOUNDS=[]
for d,ds in PATTERNS:
    aa=[1-F(1,2**(h-1)) for h in ds]
    T=prod((1+a)/2 for a in aa)-F(1,2**d)*prod(aa)
    require(T>0,'Nonpositive association constant')
    BOUNDS.append(T/(1+T))
EQNAMES=([f'm{k}' for k in (2,3,4)]+[f'q{k}' for k in (2,3,4)]+['J23','J24']+
         [f'R_{t}_{a}_{b}_{c}' for t,a,b,c in REV]+
         [f'S_{d}_{x}_'+','.join(map(str,T)) for d,x,T in STARS]+
         [f'BR_{kind}_{d}_{h}_{k}' for d,h,k in CLASSES for kind in ('X','B')])
GNAMES=([f'E{kind}_{d}_{h}_{k}' for d,h,k in CLASSES for kind in ('X','B')]+
        [f'L_{d}_'+','.join(map(str,T)) for d,T in PATTERNS])
NE=170;NG=54;NF=225
require((len(REV),len(STARS),len(CLASSES),len(PATTERNS),len(EQNAMES),len(GNAMES))==(56,60,23,8,NE,NG),'Raw dimensions')

def feature_key(v):
    d=v['d'];ds=v['degrees'];masks=v['masks'];W=len(masks);n=d+1
    require(W>0,'Empty local family')
    occ=[sum(bool(S>>i&1) for S in masks) for i in range(n)]
    nums=[12*occ[0]]
    for k in (2,3,4): nums.append(12*W*(d==k)-(12//k)*W*ds.count(k))
    for k in (2,3,4): nums.append(12*occ[0]*(d==k)-(12//k)*sum(occ[i] for i in range(1,n) if ds[i-1]==k))
    for h in (3,4):
        nums.append(12*sum(sum(bool(S&1) and bool(S>>i&1) for S in masks)*((d==2 and ds[i-1]==h)-(d==h and ds[i-1]==2)) for i in range(1,n)))
    paths=[];inside={i:[] for i in range(1,n)};outside={i:[] for i in range(1,n)}
    for a,b in v['edges']:
        inside[a].append(b);inside[b].append(a)
        paths.extend(((a,b,ds[b-1],0,1),(b,a,ds[a-1],0,1)))
    for T,(h,x,b) in zip(v['blocks'],v['choices']):
        for i in T:
            outside[i].append(2*(h-2)+x);paths.append((i,0,h,x,0))
    require(len(paths)==sum(h-1 for h in ds),'Complete nonbacktracking paths')
    rev=[0]*len(REV);stars=[0]*len(STARS)
    for S in masks:
        bits=[int(bool(S>>i&1)) for i in range(n)]
        a=2*(d-2)+bits[0];symbols=[2*(h-2)+bits[i] for i,h in enumerate(ds,1)]
        for i,j,h,x,t in paths:
            b=symbols[i-1];c=2*(h-2)+(bits[j] if j else x)
            if a<c: rev[RID[t,a,b,c]]+=12
            elif c<a: rev[RID[t,c,b,a]]-=12
        stars[SID[d,bits[0],tuple(sorted(symbols))]]+=12
        for i in range(1,n):
            U=[a]+[symbols[j-1] for j in inside[i]]+outside[i]
            require(len(U)==ds[i-1],'Neighbor star completeness')
            stars[SID[ds[i-1],bits[i],tuple(sorted(U))]]-=12//ds[i-1]
    nums.extend(rev);nums.extend(stars)
    br=[0]*46;allnb=(1<<n)-2
    for T,(h,x,b) in zip(v['blocks'],v['choices']):
        k=len(T);j=2*CID[d,h,k];r=2*CID[h,d,k]
        br[j]+=12*W*x;br[j+1]+=12*W*b;br[r]-=12*occ[0]
        required=1|(allnb^sum(1<<i for i in T))
        br[r+1]-=12*sum(S&required==required for S in masks)
    nums.extend(br);require(len(nums)==171,'Transport vector length')
    D=lcm(12*W,8,*(L.denominator for L in BOUNDS))
    nums=[x*(D//(12*W)) for x in nums]
    for rd,h,k in CLASSES:
        chosen=[(x,b) for T,(hd,x,b) in zip(v['blocks'],v['choices']) if d==rd and hd==h and len(T)==k]
        X=sum(x for x,b in chosen);B=sum(b for x,b in chosen)
        nums.extend((len(chosen)*(D//2)-X*D,X*(D//(2**(h-k)))-B*D))
    pattern=(d,tuple(sorted(ds)))
    for P,L in zip(PATTERNS,BOUNDS): nums.append(occ[0]*(D//W)-L.numerator*(D//L.denominator) if P==pattern else 0)
    require(len(nums)==NF,'Complete feature length')
    g=reduce(gcd,nums,D)
    return d,D//g,tuple(x//g for x in nums)

def load():
    with gzip.open(HERE/'columns.json.gz','rt') as f:return json.load(f)
def vector(row): return tuple(F(x,row['denominator']) for x in row['numerators'])
def rowkey(row): return row['d'],row['denominator'],tuple(row['numerators'])
def enumerate_columns():
    start=time.monotonic();col={};sha=hashlib.sha256();keys=[];views=[]
    with gzip.open(HERE/'inputs/views.jsonl.gz','rt') as source:
        for i,v in enumerate(grammar.all_views()):
            raw=canonical(v)+'\n';sha.update(raw.encode());require(raw==next(source),'Fixed view mismatch')
            k=feature_key(v);keys.append(k);views.append(v)
            if k not in col:col[k]={'representative':v,'source_indices':[]}
            col[k]['source_indices'].append(i)
        require(source.read()=='','Fixed extra views')
    require(len(keys)==17553 and sha.hexdigest()==GEN_SHA,'Grammar coverage digest')
    rows=[]
    for i,((d,D,nums),extra) in enumerate(sorted(col.items())):
        rows.append(dict(id=i,d=d,denominator=D,numerators=nums,full_views=len(extra['source_indices']),**extra))
    ids={rowkey(r):r['id'] for r in rows};mapping=[ids[k] for k in keys]
    zeroeq=[EQNAMES[j] for j in range(NE) if all(r['numerators'][j+1]==0 for r in rows)]
    zerog=[GNAMES[j] for j in range(NG) if all(r['numerators'][171+j]==0 for r in rows)]
    metadata=dict(activity='1',q='1/2',target='4/11',generation_sha256=GEN_SHA,full_views=len(keys),columns=len(rows),
                  feature_length=NF,raw_transports=NE,normalization_plus_transports=NE+1,raw_inequalities=NG,
                  equality_names=EQNAMES,inequality_names=GNAMES,reversals=REV,stars=STARS,classes=CLASSES,patterns=PATTERNS,
                  association_bounds=list(map(str,BOUNDS)),zero_transport_names=zeroeq,zero_inequality_names=zerog)
    data=dict(metadata,rows=rows,full_to_column=mapping)
    with gzip.open(HERE/'columns.json.gz','wt') as f:json.dump(data,f,sort_keys=True)
    dump(HERE/'enumeration.json',dict(metadata,elapsed_seconds=time.monotonic()-start))
    print(json.dumps({k:metadata[k] for k in ('columns','full_views','feature_length','zero_transport_names','zero_inequality_names','association_bounds')},indent=2))
    primary_fixtures(data,views)

def primary_fixtures(data,views):
    known={canonical(grammar.key_of(v)):i for i,v in enumerate(views)};reports={}
    for name,(n,edges) in grammar.fixtures().items():
        adj=[set() for _ in range(n)]
        for a,b in edges:adj[a].add(b);adj[b].add(a)
        closed=[{i}|A for i,A in enumerate(adj)]
        families=[{i for i in range(n) if S>>i&1} for S in range(1<<n) if all(any(not(S>>u&1) for u in C) for C in closed)]
        sums=[F(0)]*NF;masses=defaultdict(F);exposures=0
        for root in range(n):
            exterior=sorted(set(range(n))-({root}|adj[root]));mass=F(0)
            for B in range(1<<len(exterior)):
                E={u for j,u in enumerate(exterior) if B>>j&1}
                completed=[I for I in families if I.intersection(exterior)==E]
                if not completed:continue
                v,local=grammar.actual_view(adj,root,E);actual=sorted(sum(1<<local[u] for u in I if u in local) for I in completed)
                require(actual==list(v['masks']),'Primary fixture masks')
                raw=known[canonical(grammar.key_of(v))];cid=data['full_to_column'][raw];f=vector(data['rows'][cid])
                p=F(len(completed),n*len(families));mass+=p;exposures+=1;masses[cid]+=p
                for j,x in enumerate(f):sums[j]+=p*x
            require(mass==F(1,n),'Fixture root exposure mass')
        require(sum(masses.values())==1,'Fixture all mass')
        require(all(x==0 for x in sums[1:171]),'Fixture transport mean')
        require(all(x>=0 for x in sums[171:]),'Fixture averaged inequality')
        direct=F(sum(map(len,families)),n*len(families));require(sums[0]==direct,'Fixture occupancy')
        reports[name]=dict(vertices=n,feasible_sets=len(families),exposures=exposures,objective=str(direct),feature_means=list(map(str,sums)),masses={str(k):str(v) for k,v in masses.items()})
        if name=='C4':
            require(direct==F(4,11),'C4 target')
            cert=dict(activity='1',generation_sha256=GEN_SHA,objective='4/11',support=[dict(id=i,weight=str(w)) for i,w in sorted(masses.items())])
            dump(HERE/'c4_primal.json',cert);dump(HERE/'c4_primal_verification.json',verify_primal(cert,data))
    dump(HERE/'primary_fixtures.json',dict(all_passed=True,fixtures=reports,total_exposures=sum(r['exposures'] for r in reports.values())))

def verify_primal(cert,data):
    require(cert['activity']=='1' and cert['generation_sha256']==GEN_SHA,'Primal scope')
    sums=[F(0)]*NF;total=F(0);seen=set()
    for item in cert['support']:
        i=item['id'];w=F(item['weight']);require(0<=i<len(data['rows']) and i not in seen,'Primal id');seen.add(i)
        require(w>=0,'Negative primal weight');total+=w
        for j,x in enumerate(vector(data['rows'][i])):sums[j]+=w*x
    require(total==1,'Primal normalization');require(all(x==0 for x in sums[1:171]),'Primal transport failure')
    require(all(x>=0 for x in sums[171:]),'Primal inequality failure')
    require(sums[0]==F(cert['objective']),'Primal objective failure')
    return dict(all_passed=True,objective=str(sums[0]),support_size=len(seen),feature_means=list(map(str,sums)),total_mass=str(total))

def verify_dual(cert,data):
    require(cert['activity']=='1' and cert['generation_sha256']==GEN_SHA,'Dual scope')
    y=list(map(F,cert['equality_dual']));mu=list(map(F,cert['inequality_dual']))
    require(len(y)==171 and len(mu)==54,'Dual dimension');require(all(x>=0 for x in mu),'Negative inequality multiplier')
    require(y[0]==F(cert['objective']),'Dual objective mismatch')
    H=lcm(*(x.denominator for x in y+mu));yi=[int(H*x) for x in y];mi=[int(H*x) for x in mu]
    minimum=None;tight=0
    for row in data['rows']:
        D=row['denominator'];nums=row['numerators']
        s=H*nums[0]-yi[0]*D-sum(a*b for a,b in zip(yi[1:],nums[1:171]))-sum(a*b for a,b in zip(mi,nums[171:]))
        require(s>=0,f'Negative dual column slack {row["id"]}: {s}')
        tight+=s==0;slack=F(s,H*D);minimum=slack if minimum is None else min(minimum,slack)
    return dict(all_passed=True,objective=str(y[0]),columns=len(data['rows']),tight_columns=tight,minimum_slack=str(minimum),dual_nonzero=sum(x!=0 for x in y+mu))

def replay():
    data=load();require(data['equality_names']==EQNAMES and data['inequality_names']==GNAMES,'Metadata row order')
    seen=defaultdict(list);sha=hashlib.sha256();tight=0;full_count=0
    for i,v in enumerate(grammar.all_views()):
        sha.update((canonical(v)+'\n').encode());key=feature_key(v);cid=data['full_to_column'][i]
        require(key==rowkey(data['rows'][cid]),f'Full regenerated feature mismatch {i}');seen[cid].append(i);full_count+=1
    require(full_count==17553 and sha.hexdigest()==GEN_SHA,'Full replay coverage')
    require(len(seen)==len(data['rows']),'Column coverage')
    for r in data['rows']:
        require(seen[r['id']]==r['source_indices'] and len(seen[r['id']])==r['full_views'],'Multiplicity mismatch')
    result=dict(all_passed=True,full_views=full_count,columns=len(data['rows']),generation_sha256=sha.hexdigest(),witnesses={})
    for name in ('c4_primal','primal','certificate'):
        path=HERE/(name+'.json')
        if path.exists():
            cert=json.loads(path.read_text());v=verify_primal(cert,data)
            if 'equality_dual' in cert:v['dual']=verify_dual(cert,data)
            result['witnesses'][name]=v
    dump(HERE/'primary_replay.json',result);print(json.dumps({k:v for k,v in result.items() if k!='witnesses'},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['enumerate','replay']);a=p.parse_args()
    {'enumerate':enumerate_columns,'replay':replay}[a.command]()
