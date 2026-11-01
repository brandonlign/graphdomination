#!/usr/bin/env python3
"""Standalone exact example checks: no imports from the certificate engines.
Only Python's standard library; writes outputs to this addendum directory.
"""
from __future__ import annotations
import hashlib,json
from pathlib import Path
from itertools import combinations
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def trim(p):
 p=list(p)
 while len(p)>1 and p[-1]==0:p.pop()
 return p
def add(a,b):return trim([(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return trim(out)
def polynomial(sets,weight=lambda s:1):
 out=[0]*(max(map(len,sets),default=0)+1)
 for s in sets:out[len(s)]+=weight(s)
 return trim(out)
def subsets(vs):
 vs=tuple(vs)
 return [frozenset(c)for k in range(len(vs)+1)for c in combinations(vs,k)]
def complete_bipartite(left,right):
 return {v:frozenset(right if v in left else left)for v in left+right}
def feasible(g,I):return all(not({v}|g[v])<=I for v in g)
def sigma(g,v,I):return 2*(len(g[v])-2)+int(v in I)
def reversal(g,v,I,t,a,b,c):
 ans=0
 for u in g[v]:
  for w in g[u]-{v}:
   if int(w in g[v])!=t:continue
   sy=(sigma(g,v,I),sigma(g,u,I),sigma(g,w,I))
   ans+=int(sy==(a,b,c))-int(sy==(c,b,a))
 return ans

g=complete_bipartite(('a','ap'),('b','c','d'))
Fs=[I for I in subsets(g)if feasible(g,I)]
Z=polynomial(Fs);O={v:polynomial(Fs,lambda I:int(v in I))for v in g}
need(Z==[1,5,10,7],'K2,3 partition')
need(O['a']==O['ap']==[0,1,4,3],'cubic marginal')
need(O['b']==O['c']==O['d']==[0,1,4,5],'quadratic marginal')
A=[1,4,6];target=[0,1,3]
need(add(mul(O['a'],A),scale(mul(target,Z),-1))==[0,0,0,0,-1,-3],'cubic difference identity')
need(add(mul(O['b'],A),scale(mul(target,Z),-1))==[0,0,0,2,7,9],'degree2 difference identity')
agg=add(scale(O['a'],2),scale(O['b'],3))
need(agg==[0,5,20,21],'global occupancy numerator')
need(add(mul(agg,A),scale(mul(target,Z),-5))==[0,0,0,6,19,21],'global difference identity')
need(Q(sum(O['a']),sum(Z))==Q(8,23) and Q(sum(O['b']),sum(Z))==Q(10,23),'unit marginals')
need(Q(sum(agg),5*sum(Z))==Q(2,5),'unit density')

exposures=[];sumZ=[0];sumO=[0];sumR=[0]
for E in subsets(('c','d')):
 local=[S for S in subsets(('b','a','ap'))if feasible(g,S|E)]
 zc=polynomial(local);oc=polynomial(local,lambda S:int('b'in S))
 r=polynomial(local,lambda S:reversal(g,'b',S|E,0,0,2,1))
 q2=oc
 q3=polynomial(local,lambda S:-Q(int('a'in S)+int('ap'in S),3))
 j=polynomial(local,lambda S:int('b'in S)*(int('a'in S)+int('ap'in S)))
 k=len(E)
 need(zc==[[1,3,3],[1,3,2],[1,3]][k],'exposure Z')
 need(oc==[[0,1,2],[0,1,2],[0,1]][k],'exposure O')
 if E=={'c'}:
  need(set(local)==set(map(frozenset,[(),('b',),('a',),('ap',),('b','a'),('b','ap')])),'explicit local feasible list')
  need(r==[2,0,-2],'worked reversal polynomial')
  need(q3==[Q(0),Q(-2,3),Q(-2,3)],'worked q3')
  need(j==[0,0,2],'worked j')
 shift=[0]*k+[1]
 sumZ=add(sumZ,mul(shift,zc));sumO=add(sumO,mul(shift,oc));sumR=add(sumR,mul(shift,r))
 exposures.append(dict(E=sorted(E),feasible_inside=[sorted(S)for S in local],Z=zc,O=oc,R_0_0_2_1_numerator=r,q2_numerator=q2,q3_numerator=list(map(str,q3)),j_numerator=j,conditional_mass_at_one=str(Q(sum(zc),sum(Z))),root_exposure_mass_at_one=str(Q(sum(zc),5*sum(Z)))))
need(sumZ==Z and sumO==O['b'],'exposure reconstruction')
for I in Fs:
 for k in(2,3):
  m=sum(Q(len(g[v])==k)-sum((Q(1,k)for u in g[v]if len(g[u])==k),Q(0))for v in g)
  q=sum(Q(len(g[v])==k and v in I)-sum((Q(u in I,k)for u in g[v]if len(g[u])==k),Q(0))for v in g)
  need(m==q==0,'transport root sums')
 need(sum(reversal(g,v,I,0,0,2,1)for v in g)==0,'selected reversal root sum')

families=[]
for name,H in [('K4',{v:frozenset(set(range(4))-{v})for v in range(4)}),('K3,3',complete_bipartite((0,1,2),(3,4,5)))]:
 edges=sorted((u,v)for u in H for v in H[u]if u<v);N=len(H)
 S={('old',v):set()for v in H}
 for i,(u,v)in enumerate(edges):
  S[('new',i)]={('old',u),('old',v)}
  S[('old',u)].add(('new',i));S[('old',v)].add(('new',i))
 n2=sum(len(a)==2 for a in S.values());n3=sum(len(a)==3 for a in S.values())
 need(n2==3*N//2 and n3==N and Q(n2,len(S))==Q(3,5),'subdivision profile')
 need(all(u[0]!=v[0]for u in S for v in S[u]),'subdivision bipartition')
 # Every cycle balances bipartition counts, so different part sizes exclude a spanning 2-factor.
 need(n2!=n3 and Q(n2,len(S))>Q(7,187),'profile threshold and part imbalance')
 families.append(dict(cubic_graph=name,N=N,n2=n2,n3=n3,n=len(S),bipartition_sizes=[N,len(edges)],degree_two_fraction='3/5',no_spanning_2_factor_by_part_imbalance=True))
report=dict(all_passed=True,arithmetic='integer polynomial coefficients and fractions.Fraction',independent_of_certificate_implementations=True,global_feasible_sets=len(Fs),Z=Z,vertex_numerators=O,global_occupancy_numerator=agg,exposures=exposures,subdivision_examples=families,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(ROOT/'EXAMPLE_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
