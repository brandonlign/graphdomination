#!/usr/bin/env python3
"""Ordinary-set controls for the all-length singleton-exposure argument.
This checks C3,...,C20; the generic reasoning is stated in the article.
No local-view generator, numerical optimizer or certificate evaluator imported.
"""
if not __debug__:raise RuntimeError('Cycle exposure check refuses optimized Python')
from pathlib import Path
from itertools import combinations
import json

def require(c,m):
 if not c:raise RuntimeError(m)
def sets(S):
 for k in range(len(S)+1):
  for T in combinations(sorted(S),k):yield set(T)
def feasible(adj,I):return all(not ({u}|adj[u])<=I for u in adj)
def view(adj,v,E):
 ns=sorted(adj[v]);R={v}|set(ns);boundary=set().union(*(adj[u]for u in ns))-R
 edges=[(u,w)for u,w in combinations(ns,2)if w in adj[u]]
 blocks=[dict(vertex=w,T=sorted(adj[w]&set(ns)),h=len(adj[w]),x=int(w in E),b=int(({w}|(adj[w]-R))<=E))for w in sorted(boundary)]
 return R,edges,blocks
records=[]
for m in range(3,21):
 adj={u:{(u-1)%m,(u+1)%m}for u in range(m)};E=set()if m==3 else{2};R,edges,blocks=view(adj,0,E)
 require(feasible(adj,E),'Singleton not feasible');inside=[S for S in sets(R)if feasible(adj,E|S)]
 require(set()in inside,'Missing empty completion')
 zero=not edges and (all(b['x']==b['b']==0 for b in blocks)or(len(blocks)==1 and blocks[0]['h']==2 and len(blocks[0]['T'])==2 and blocks[0]['x']==blocks[0]['b']==1))
 if m==3:require(edges and not zero,'C3 root-neighbor adjacency not detected')
 elif m==4:require(zero,'C4 must have shared-boundary zero view')
 else:
  require(not zero and any(b['vertex']==2 and b['x']==1 and b['b']==0 for b in blocks),'Singleton does not leave both zero supports')
  require(3 not in R and 3 not in E,'Absent exterior neighbor explaining b=0 missing')
  Z=[sum(len(S)==k for S in inside)for k in range(3)];O=[sum(len(S)==k and 0 in S for S in inside)for k in range(3)]
  require(Z==[1,3,2]and O==[0,1,1],'Singleton completion polynomial differs')
 records.append(dict(m=m,exterior=sorted(E),neighbor_edges=edges,boundary=blocks,in_zero_shape=zero,feasible_inside=[sorted(S)for S in inside]))
out=dict(all_passed=True,tested_cycles=list(range(3,21)),project_generator_imported=False,
 generic_argument='For every m>=5, vertices 0,1,2,3 are distinct; 3 is not in N[0]. I={2} is feasible because every closed neighborhood has three vertices. At root 0, boundary vertex 2 has x=1,b=0 because its outside closed part contains absent vertex 3. Both verified zero supports exclude that flag. The actual empty-inside completion has positive weight at every z>0. C3 instead has adjacent root neighbors, also excluded.',records=records)
Path(__file__).with_name('CYCLE_EXPOSURES.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='records'},indent=2))
