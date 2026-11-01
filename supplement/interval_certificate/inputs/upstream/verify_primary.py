#!/usr/bin/env python3
"""Optimizer-free exact integer Bernstein replay with primary full regeneration."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Primary interval exact replay refuses optimized Python')
import sys,json,hashlib,time,math
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import symbolic_engine as s

def verify(path,regenerate=True):
 start=time.time();cert=json.loads(Path(path).read_text());data=s.load();a,b=F(cert['a']),F(cert['b']);d=cert['degree'];M=cert['bernstein_degree']
 assert cert['support']==s.SUPPORT and cert['generation_sha256']==data['generation_sha256']
 allowed=[(F(0),F(1),2,4),(F(0),F(1),4,8),(F(0),F(1),6,12),(F(0),F(1,2),6,12),(F(1,2),F(1),6,12),(F(1,4),F(1,2),6,12),(F(1,2),F(3,4),6,12),(F(3,4),F(1),6,12)]
 assert(a,b,d,cert['elevation'])in allowed and M==10+d+cert['elevation']
 vals=[[F(x)for x in row]for row in cert['bernstein_polynomials']];assert len(vals)==74 and all(len(v)==d+1 for v in vals)
 H=math.lcm(*(v.denominator for row in vals for v in row));assert H>0
 P=[];mult=0
 for name,idx,row in zip(s.SUPPORT,s.IND,vals):
  if idx>=114:assert all(x>=0 for x in row);mult+=len(row)
  beta=[int(H*x)for x in row];power=[0]*(d+1)
  for i,v in enumerate(beta):
   for j in range(d-i+1):power[i+j]+=v*math.comb(d,i)*math.comb(d-i,j)*(-1)**j
  P.append(s.trim(power))
 B=math.lcm(a.denominator,b.denominator);av=int(a*B);delta=int((b-a)*B);assert B>0 and delta>0
 @lru_cache(None)
 def change(p):
  assert len(p)<=11
  return s.trim([sum(v*math.comb(j,k)*av**(j-k)*delta**k*B**(10-j)for j,v in enumerate(p)if j>=k)for k in range(11)])
 combs=[math.comb(M,j)for j in range(11+d)];L=math.lcm(*combs)
 matrix=[[math.comb(i,j)*(L//combs[j])if j<=i else 0 for j in range(11+d)]for i in range(M+1)]
 for N,D in s.association():assert D[0]>0 and all(x>=0 for x in D)
 known={};sha=hashlib.sha256();checks=0;zeros=0;fullzeros=0;c4=0
 for row in data['rows']:
  key=s.key_from_row(row);assert row['Z'][0]==1 and all(x>=0 for x in row['Z']);target,polys=s.cleared(row)
  slack=s.scale(change(tuple(target)),H)
  for p,q in zip(P,polys):slack=s.add(slack,s.scale(s.mul(p,change(tuple(q))),-1))
  bern=[sum(v*mat[j]for j,v in enumerate(slack))for mat in matrix];assert all(x>=0 for x in bern),(row['id'],min(bern))
  checks+=len(bern);iszero=not any(slack);zeros+=iszero;fullzeros+=iszero*row['full_views'];known[key]=row['full_views']
  rep=row['representative']
  if rep['d']==2 and rep['neighbor_degrees']==[2,2]and rep['blocks']==[[1,2]]and rep['boundary_degrees']==[2]and rep['neighbor_edges']==[]:assert iszero;c4+=1
  sha.update(('|'.join(map(str,bern))+'\n').encode())
 assert c4==2 and len(known)==9867
 full=0;generation=hashlib.sha256()
 if regenerate:
  seen={}
  for rep in s.structural.all_full_views():
   key=s.feature_key(rep);assert key in known;seen[key]=seen.get(key,0)+1;full+=1
   generation.update((json.dumps(rep,sort_keys=True,separators=(',',':'))+'\n').encode())
  assert seen==known and full==231816 and generation.hexdigest()==data['generation_sha256']
 out=dict(all_passed=True,certificate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),interval=[str(a),str(b)],columns=len(known),full_views=full,regenerated=regenerate,
 slack_coefficients_checked=checks,multiplier_coefficients_checked=mult,zero_slack_columns=zeros,zero_slack_full_views=fullzeros,c4_identical_zero=2,
 all_denominators_positive=True,integer_bernstein_sequence_sha256=sha.hexdigest(),elapsed_seconds=time.time()-start,optimizer_used=False)
 dest=Path(path).parent/'primary_exact_replay.json';s.dump(dest,out);print(json.dumps(out,indent=2),flush=True);return out
if __name__=='__main__':verify(sys.argv[1])
