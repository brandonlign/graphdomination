#!/usr/bin/env python3
from __future__ import annotations
if not __debug__: raise RuntimeError('Symbolic interval proof refuses optimized Python')
import sys,json,gzip,hashlib,time
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'inputs/dependencies'))
import v2_features_snapshot as structural
import local_views_v3 as u
NAMES=['normalization']+u.EQNAMES[1:]+u.GNAMES
SUPPORT=[x for x in json.loads((HERE/'inputs/SUPPORT_COMPARISON.json').read_text())['union_coordinates'] if x!='normalization']
assert len(SUPPORT)==74
IND=[NAMES.index(x)for x in SUPPORT]
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def trim(v):
 v=list(v)
 while len(v)>1 and v[-1]==0:v.pop()
 return tuple(v)
def add(a,b):return trim([(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))])
def scale(a,b):return trim([x*b for x in a])
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)
def peval(p,z):return sum((x*z**i for i,x in enumerate(p)),F(0))
@lru_cache(None)
def association():
 import sympy as s
 z=s.Symbol('z');q=z/(1+z);out=[]
 for d,r in u.PATTERNS:
  ds=[3]*r+[2]*(d-r)if d==2 else[2]*r+[3]*(d-r)
  aa=[1-q**(h-1)for h in ds]
  T=s.prod(1-q+q*a for a in aa)-q**d*s.prod(aa)
  N,D=s.fraction(s.cancel(z*T/(1+z*T)))
  np=trim([int(s.expand(N).coeff(z,i))for i in range(s.degree(N,z)+1)])
  dp=trim([int(s.expand(D).coeff(z,i))for i in range(s.degree(D,z)+1)])
  if dp[0]<0:np=scale(np,-1);dp=scale(dp,-1)
  assert dp[0]>0 and all(x>=0 for x in dp)
  out.append((np,dp))
 return tuple(out)
def feature_key(rep):
 d=rep['d'];ds=rep['neighbor_degrees'];sets=rep['feasible_masks'];n=d+1
 Z=[0]*4;occ=[[0]*4 for _ in range(n)];joint=[[0]*4 for _ in range(d)]
 for S in sets:
  k=S.bit_count();assert k<4;Z[k]+=1
  for a in range(n):occ[a][k]+=bool(S&(1<<a))
  for a in range(1,n):joint[a-1][k]+=bool(S&1 and S&(1<<a))
 nums=[]
 nums.append([6*x for x in occ[0]])
 for h in(2,3):nums.append([(6*int(d==h)-(6//h)*sum(k==h for k in ds))*w for w in Z])
 for h in(2,3):nums.append([6*occ[0][k]*int(d==h)-(6//h)*sum(occ[a][k]for a in range(1,n)if ds[a-1]==h)for k in range(4)])
 nums.append([6*sum((int(d==2 and ds[a]==3)-int(d==3 and ds[a]==2))*joint[a][k]for a in range(d))for k in range(4)])
 paths=[];inside={a:[]for a in range(1,n)};outside={a:[]for a in range(1,n)}
 for a,b in rep['neighbor_edges']:
  paths.extend([(a,b,ds[b-1],0,1),(b,a,ds[a-1],0,1)]);inside[a].append(b);inside[b].append(a)
 for T,h,(x,b)in zip(rep['blocks'],rep['boundary_degrees'],rep['flags']):
  for a in T:paths.append((a,0,h,x,0));outside[a].append(2*(h-2)+x)
 assert len(paths)==sum(h-1 for h in ds)
 rev=[[0]*4 for _ in range(48)];star=[[0]*4 for _ in range(60)]
 for S in sets:
  k=S.bit_count();bits=[int(bool(S&(1<<a)))for a in range(n)];a=2*(d-2)+bits[0];symbols=[2*(ds[v-1]-2)+bits[v]for v in range(1,n)]
  for v,w,h,x,t in paths:
   b=symbols[v-1];c=2*(h-2)+(bits[w]if w else x)
   if a<c:rev[structural.RID[t,a,b,c]][k]+=6
   elif c<a:rev[structural.RID[t,c,b,a]][k]-=6
  star[u.SID[(d,bits[0],tuple(sorted(symbols)))]][k]+=6
  for v in range(1,n):
   ns=[a]+[symbols[w-1]for w in inside[v]]+outside[v];assert len(ns)==ds[v-1]
   star[u.SID[(ds[v-1],bits[v],tuple(sorted(ns)))]][k]-=6//ds[v-1]
 nums+=rev+star;assert len(nums)==114
 ext=[]
 for rd,h,k in u.CLASSES:
  chosen=[(x,b)for T,hd,(x,b)in zip(rep['blocks'],rep['boundary_degrees'],rep['flags'])if d==rd and hd==h and len(T)==k]
  ext.append((len(chosen),sum(x for x,b in chosen),sum(b for x,b in chosen)))
 r=sum(h==(3 if d==2 else 2)for h in ds)
 return d,tuple(Z),tuple(tuple(v)for v in nums),tuple(ext),u.PATTERNS.index((d,r))
def key_from_row(row):return row['d'],tuple(row['Z']),tuple(tuple(p)for p in row['features6']),tuple(tuple(t)for t in row['exterior']),row['pattern']
def features_at(key,z):
 d,Z,nums,ext,p=key;W=peval(Z,z);out=[peval(v,z)/(6*W)for v in nums];q=z/(1+z)
 for (rd,h,k),(N,X,B)in zip(u.CLASSES,ext):out +=[q*N-X,q**(h-k)*X-B]
 for j,(N,D)in enumerate(association()):out.append(out[0]-peval(N,z)/peval(D,z)if p==j else F(0))
 return tuple(out)
def cleared(row):
 d,Z,nums,ext,p=key_from_row(row);O=tuple(x//6 for x in nums[0]);A=(1,4,6)
 target=scale(add(mul(A,O),scale(mul((0,1,3),Z),-1)),6);polys=[]
 for name,idx in zip(SUPPORT,IND):
  if idx<114:poly=nums[idx]
  elif idx<134:
   ei=(idx-114)//2;rd,h,k=u.CLASSES[ei];N,X,B=ext[ei]
   if (idx-114)%2==0:ex=add(scale((0,1,1),N),scale((1,2,1),-X))
   else:ex=add(scale(mul((0,)* (h-k)+(1,), (1,)if h-k==2 else (1,1)if h-k==1 else(1,2,1)),X),scale((1,2,1),-B))
   poly=scale(mul(ex,Z),6)
  else:
   j=idx-134
   if p==j:
    N,D=association()[j];poly=scale(add(mul(D,O),scale(mul(N,Z),-1)),6)
   else:poly=(0,)
  polys.append(trim(poly))
 return target,polys

def enumerate_symbolic():
 start=time.time();cols={};sha=hashlib.sha256();counts=Counter()
 for rep in structural.all_full_views():
  key=feature_key(rep)
  if key not in cols:cols[key]=dict(representative=rep,full_views=0)
  cols[key]['full_views']+=1;counts['full_views']+=1;counts['degree_'+str(rep['d'])]+=1
  sha.update((json.dumps(rep,sort_keys=True,separators=(',',':'))+'\n').encode())
  if counts['full_views']%50000==0:print(counts['full_views'],time.time()-start,flush=True)
 assert counts['full_views']==231816
 expected=json.loads((HERE/'inputs/dependencies/enumeration_v3.json').read_text())
 assert sha.hexdigest()==expected['generation_sha256']
 rows=[dict(id=i,d=k[0],Z=k[1],features6=k[2],exterior=k[3],pattern=k[4],**v)for i,(k,v)in enumerate(sorted(cols.items()))]
 out=dict(rows=rows,columns=len(rows),counts=dict(counts),generation_sha256=sha.hexdigest(),support=SUPPORT,association=[dict(N=N,D=D,pattern=p)for(N,D),p in zip(association(),u.PATTERNS)],elapsed_seconds=time.time()-start)
 with gzip.open(HERE/'symbolic_columns.json.gz','wt')as f:json.dump(out,f,separators=(',',':'))
 dump(HERE/'enumeration.json',{k:v for k,v in out.items()if k!='rows'})
 print('DONE',len(rows),time.time()-start,flush=True)
def load():
 with gzip.open(HERE/'symbolic_columns.json.gz','rt')as f:return json.load(f)
if __name__=='__main__':enumerate_symbolic()
