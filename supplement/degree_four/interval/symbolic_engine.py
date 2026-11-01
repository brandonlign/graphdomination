#!/usr/bin/env python3
"""Complete z-weighted degree-four features; optimizer-free exact enumeration."""
from __future__ import annotations
if not __debug__: raise RuntimeError('Exact polynomial checks refuse Python -O/-OO')
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb
from collections import defaultdict,Counter
import sys,json,gzip,hashlib,time,importlib.util
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
GEN_SHA='af0020c5efce7b021f901e798013b34aa7967c019c4e5c98b8352e0e2985f021'
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def dump(path,x):Path(path).write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def trim(a):
 a=list(a)
 while len(a)>1 and a[-1]==0:a.pop()
 return tuple(a or [0])
def add(a,b):return trim([(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([c*v for v in a])
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)
def shift(a,k):return (0,)*k+tuple(a) if any(a)else(0,)
def oneplus(k):return tuple(comb(k,i)for i in range(k+1))
def low(s):return s<2
def edgeok(a,b):return low(a)or low(b)
REV=tuple((t,a,b,c)for t in(0,1)for a in range(6)for b in range(6)for c in range(a+1,6)if edgeok(a,b)and edgeok(b,c)and(not t or edgeok(a,c)))
STARS=tuple((d,x,T)for d in(2,3,4)for x in(0,1)for T in combinations_with_replacement(range(6 if d==2 else 2),d))
CLASSES=tuple((d,h,k)for d in(2,3,4)for h in(2,3,4)for k in range(1,min(d,h)+1))
PATTERNS=tuple((d,T)for d in(2,3,4)for T in combinations_with_replacement((2,3,4)if d==2 else(2,),d))
RID={v:i for i,v in enumerate(REV)};SID={v:i for i,v in enumerate(STARS)};CID={v:i for i,v in enumerate(CLASSES)}
EQNAMES=([f'm{k}'for k in(2,3,4)]+[f'q{k}'for k in(2,3,4)]+['J23','J24']+[f'R_{t}_{a}_{b}_{c}'for t,a,b,c in REV]+[f'S_{d}_{x}_'+','.join(map(str,T))for d,x,T in STARS]+[f'BR_{kind}_{d}_{h}_{k}'for d,h,k in CLASSES for kind in('X','B')])
GNAMES=([f'E{kind}_{d}_{h}_{k}'for d,h,k in CLASSES for kind in('X','B')]+[f'L_{d}_'+','.join(map(str,T))for d,T in PATTERNS])
SUPPORT=EQNAMES+GNAMES
require((len(REV),len(STARS),len(CLASSES),len(PATTERNS),len(SUPPORT))==(56,60,23,8,224),'Dimensions')
spec=importlib.util.spec_from_file_location('view_grammar',HERE/'inputs/enumerate_validate.py');grammar=importlib.util.module_from_spec(spec);spec.loader.exec_module(grammar)

def association():
 out=[]
 for d,hs in PATTERNS:
  aa=[h-1 for h in hs];M=d+sum(aa)-min(aa);R=(0,)
  for mask in range((1<<d)-1):
   selected=[i for i in range(d)if mask>>i&1];e=M-d-sum(aa[i]for i in selected)
   require(e>=0,'Association exponent')
   term=shift(oneplus(e),len(selected))
   for i in selected:term=mul(term,oneplus(aa[i])[:-1]) # (1+z)^a-z^a
   R=add(R,term)
  N=shift(R,1);D=add(oneplus(M),N)
  require(R[0]==1 and D[0]==1 and all(x>=0 for x in R+N+D),'Association coefficients')
  out.append(dict(M=M,R=R,N=N,D=D))
 return out
ASSOCIATION=association()

def feature(v):
 d=v['d'];hs=v['degrees'];n=d+1
 Z=[0]*n;O=[0]*n;polys=[[0]*n for _ in range(170)]
 internal={i:[]for i in range(1,n)};external={i:[]for i in range(1,n)};walks=[]
 for i,j in v['edges']:
  internal[i].append(j);internal[j].append(i);walks.extend([(i,j,hs[j-1],0,1),(j,i,hs[i-1],0,1)])
 ext=[[0,0,0]for _ in CLASSES]
 for T,(h,x,b)in zip(v['blocks'],v['choices']):
  j=CID[d,h,len(T)];ext[j][0]+=1;ext[j][1]+=x;ext[j][2]+=b
  for i in T:external[i].append(2*(h-2)+x);walks.append((i,0,h,x,0))
 require(len(walks)==sum(h-1 for h in hs),'Walk completeness')
 for mask in v['masks']:
  k=mask.bit_count();require(k<=d,'Root clause degree bound');Z[k]+=1
  bits=[(mask>>i)&1 for i in range(n)];O[k]+=bits[0]
  f=[12*(d==h)-(12//h)*hs.count(h)for h in(2,3,4)]
  f += [12*bits[0]*(d==h)-(12//h)*sum(bits[i]for i in range(1,n)if hs[i-1]==h)for h in(2,3,4)]
  f += [12*bits[0]*sum(bits[i]*((d==2 and hs[i-1]==h)-(d==h and hs[i-1]==2))for i in range(1,n))for h in(3,4)]
  rev=[0]*56;star=[0]*60;a=2*(d-2)+bits[0];symbols=[2*(h-2)+bits[i]for i,h in enumerate(hs,1)]
  for i,j,h,x,t in walks:
   mid=symbols[i-1];c=2*(h-2)+(bits[j]if j else x)
   if a<c:rev[RID[t,a,mid,c]]+=12
   elif c<a:rev[RID[t,c,mid,a]]-=12
  star[SID[d,bits[0],tuple(sorted(symbols))]]+=12
  for i in range(1,n):
   U=[a]+[symbols[j-1]for j in internal[i]]+external[i];require(len(U)==hs[i-1],'Neighbor star degree')
   star[SID[hs[i-1],bits[i],tuple(sorted(U))]]-=12//hs[i-1]
  br=[0]*46
  for T,(h,x,b)in zip(v['blocks'],v['choices']):
   j=2*CID[d,h,len(T)];r=2*CID[h,d,len(T)]
   br[j]+=12*x;br[j+1]+=12*b;br[r]-=12*bits[0]
   br[r+1]-=12*bits[0]*all(bits[i]for i in range(1,n)if i not in T)
  f+=rev+star+br;require(len(f)==170,'Feature dimensions')
  for j,value in enumerate(f):polys[j][k]+=value
 require(Z[0]==1,'Positive local constant')
 return dict(d=d,Z=trim(Z),O=trim(O),F=tuple(map(trim,polys)),exterior=tuple(map(tuple,ext)),pattern=PATTERNS.index((d,tuple(sorted(hs)))))

def key(row):return(row['d'],tuple(row['Z']),tuple(row['O']),tuple(map(tuple,row['F'])),tuple(map(tuple,row['exterior'])),row['pattern'])
def cleared(row):
 Z=row['Z'];O=row['O'];target=scale(add(mul((1,4,6),O),scale(mul((0,1,3),Z),-1)),12);ps=[tuple(p)for p in row['F']]
 for(d,h,k),(n,x,b)in zip(CLASSES,row['exterior']):
  qx=mul(oneplus(2),(-x,n-x));qb=add(scale(shift(oneplus(3-h+k),h-k),x),scale(oneplus(3),-b))
  ps.extend((scale(mul(Z,qx),12),scale(mul(Z,qb),12)))
 for j,p in enumerate(ASSOCIATION):ps.append(scale(add(mul(p['D'],O),scale(mul(p['N'],Z),-1)),12)if j==row['pattern']else(0,))
 require(len(ps)==224,'Cleared dimension');return target,ps

def load():
 with gzip.open(HERE/'symbolic_columns.json.gz','rt')as f:return json.load(f)
def c4_column_ids(data):return data['c4_columns']
def point_vector(row):
 W=sum(row['Z']);a=[Q(sum(row['O']),W)]+[Q(sum(f),12*W)for f in row['F']]
 for(d,h,k),(n,x,b)in zip(CLASSES,row['exterior']):a.extend([Q(n,2)-x,Q(x,2**(h-k))-b])
 for j,p in enumerate(ASSOCIATION):a.append(a[0]-Q(sum(p['N']),sum(p['D']))if j==row['pattern']else Q(0))
 return a

def enumerate_columns():
 start=time.monotonic();h=hashlib.sha256();found={};rawkeys=[];views=[];point= None
 with gzip.open(HERE/'inputs/columns.json.gz','rt')as f:point=json.load(f)
 with gzip.open(HERE/'inputs/views.jsonl.gz','rt')as src:
  for i,v in enumerate(grammar.all_views()):
   raw=canonical(v)+'\n';h.update(raw.encode());require(raw==next(src),'Fixed grammar differs')
   row=feature(v);k=key(row);rawkeys.append(k);views.append(v)
   old=point['rows'][point['full_to_column'][i]];require(point_vector(row)==[Q(x,old['denominator'])for x in old['numerators']],f'Point validation {i}')
   if k not in found:found[k]=dict(row,representative=v,source_indices=[])
   found[k]['source_indices'].append(i)
  require(src.read()=='','Extra source views')
 require(len(rawkeys)==17553 and h.hexdigest()==GEN_SHA,'Universe digest')
 rows=[]
 for i,k in enumerate(sorted(found)):
  row=found[k];rows.append(dict(row,id=i,full_views=len(row['source_indices'])))
 ids={key(r):r['id']for r in rows};mapping=[ids[k]for k in rawkeys]
 c4=[]
 for r in rows:
  v=r['representative']
  if v['d']==2 and tuple(v['degrees'])==(2,2)and not v['edges']and tuple(map(tuple,v['blocks']))==((1,2),)and tuple(v['choices'][0])in((2,0,0),(2,1,1)):c4.append(r['id'])
 require(len(c4)==2,'C4 columns')
 meta=dict(generation_sha256=GEN_SHA,full_views=17553,columns=len(rows),transports=170,inequalities=54,equality_names=EQNAMES,inequality_names=GNAMES,support=SUPPORT,association=ASSOCIATION,patterns=PATTERNS,classes=CLASSES,c4_columns=c4,
 zero_transport_names=[EQNAMES[j]for j in range(170)if all(not any(r['F'][j])for r in rows)],
 zero_point_transports=point['zero_transport_names'])
 data=dict(meta,rows=rows,full_to_column=mapping)
 with gzip.open(HERE/'symbolic_columns.json.gz','wt')as f:json.dump(data,f,sort_keys=True)
 dump(HERE/'enumeration.json',dict(meta,elapsed_seconds=time.monotonic()-start,max_cleared_degree=max(len(p)-1 for r in rows for p in [cleared(r)[0]]+cleared(r)[1])))
 print(json.dumps({k:meta[k]for k in('columns','full_views','c4_columns','zero_transport_names')},indent=2),flush=True)
 primary_fixtures(data,views)

def primary_fixtures(data,views):
 known={canonical(grammar.key_of(v)):i for i,v in enumerate(views)};reports={}
 for name,(n,edges)in grammar.fixtures().items():
  adj=[set()for _ in range(n)]
  for a,b in edges:adj[a].add(b);adj[b].add(a)
  closed=[{v}|A for v,A in enumerate(adj)];family=[{i for i in range(n)if S>>i&1}for S in range(1<<n)if all(any(not(S>>u&1)for u in C)for C in closed)]
  Zg=trim([sum(len(I)==k for I in family)for k in range(n+1)]);den=scale(Zg,n)
  total=(0,);occ=(0,);trans=[(0,)]*170;inequalities=[(0,)]*54;exposures=0;exposure_records=[]
  for root in range(n):
   outer=sorted(set(range(n))-({root}|adj[root]));root_mass=(0,)
   for B in range(1<<len(outer)):
    E={v for j,v in enumerate(outer)if B>>j&1};completed=[I for I in family if I.intersection(outer)==E]
    if not completed:continue
    v,local=grammar.actual_view(adj,root,E);actual=sorted(sum(1<<local[u]for u in I if u in local)for I in completed)
    require(actual==list(v['masks']),'Fixture full mask equality')
    raw=known[canonical(grammar.key_of(v))];r=data['rows'][data['full_to_column'][raw]];require(key(feature(v))==key(r),'Fixture map')
    exp=len(E);mass=shift(r['Z'],exp);root_mass=add(root_mass,mass);occ=add(occ,shift(r['O'],exp));exposures+=1
    direct=trim([sum(len(I)==k for I in completed)for k in range(n+1)]);require(mass==direct,'Polynomial exposure identity')
    for j,Fj in enumerate(r['F']):trans[j]=add(trans[j],shift(Fj,exp))
    for j,((d,h,k),(nc,x,b))in enumerate(zip(CLASSES,r['exterior'])):
     inequalities[2*j]=add(inequalities[2*j],shift(mul(r['Z'],(-x,nc-x)),exp))
     inequalities[2*j+1]=add(inequalities[2*j+1],shift(mul(r['Z'],add(shift((x,),h-k),scale(oneplus(h-k),-b))),exp))
    j=r['pattern'];a=ASSOCIATION[j];inequalities[46+j]=add(inequalities[46+j],shift(add(mul(a['D'],r['O']),scale(mul(a['N'],r['Z']),-1)),exp))
    exposure_records.append(dict(root=root,exterior=sorted(E),raw=raw,column=r['id'],numerator=mass))
   require(root_mass==Zg,'Root mass polynomial');total=add(total,root_mass)
  require(total==den and all(not any(t)for t in trans),'Fixture polynomial transport or normalization')
  direct=trim([k*Zg[k]for k in range(len(Zg))]);require(occ==direct,'Fixture occupancy polynomial')
  require(all(all(c>=0 for c in p)for p in inequalities[:46]),'Deletion coefficient proof control')
  reports[name]=dict(vertices=n,feasible_sets=len(family),exposures=exposures,Z=Zg,occupancy_numerator=occ,transport_numerators=trans,inequality_numerators=inequalities,exposure_records=exposure_records)
  if name=='C4':
   require(Zg==(1,4,6)and occ==(0,4,12),'C4 target identity')
   dump(HERE/'c4_calibration.json',dict(all_passed=True,Z=Zg,occupancy_numerator=occ,columns=data['c4_columns'],views=[{k:data['rows'][j][k]for k in('Z','O','representative')}for j in data['c4_columns']]))
 dump(HERE/'primary_fixtures.json',dict(all_passed=True,total_exposures=sum(r['exposures']for r in reports.values()),fixtures=reports))
 print('primary fixture polynomial identities passed',flush=True)

def full_replay():
 data=load();seen=defaultdict(list);h=hashlib.sha256()
 require(data['support']==SUPPORT and data['generation_sha256']==GEN_SHA,'Replay scope')
 for i,v in enumerate(grammar.all_views()):
  h.update((canonical(v)+'\n').encode());r=data['rows'][data['full_to_column'][i]];require(key(feature(v))==key(r),f'Full feature mismatch {i}');seen[r['id']].append(i)
 require(i+1==17553 and h.hexdigest()==GEN_SHA,'Replay grammar')
 for r in data['rows']:require(r['source_indices']==seen[r['id']] and r['full_views']==len(seen[r['id']]),'Replay multiplicity')
 out=dict(all_passed=True,full_views=i+1,columns=len(data['rows']));dump(HERE/'primary_full_replay.json',out);return out
if __name__=='__main__':
 if sys.argv[1]=='enumerate':enumerate_columns()
 elif sys.argv[1]=='replay':print(full_replay())
 else:raise RuntimeError('Unknown command')
