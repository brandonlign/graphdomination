#!/usr/bin/env python3
"""Independent full-view iff classification for the final low-bit certificate."""
if not __debug__:raise RuntimeError('Independent structural zero replay refuses optimized Python')
import sys,json,gzip,hashlib,time
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent;SRC=ROOT/'inputs/upstream'
sys.modules['scipy']=None;sys.modules['scipy.optimize']=None
sys.path.insert(0,str(SRC/'audit'))
import symbolic_features_independent as ind
import verify_interval_independent as alg

def predicate(view):
 d,deg,edges,blocks=view
 if d!=2 or tuple(deg)!=(2,2) or edges:return False
 if all(x==0 and b==0 for T,h,x,b in blocks):return True
 return blocks==(((1,2),2,1,1),)

def main():
 start=time.time();cp=ROOT/'subcubic_certificate.json';digest=hashlib.sha256(cp.read_bytes()).hexdigest()
 report=json.loads((ROOT/'audit/final_average_1000/zero_support.json').read_text())
 assert report['all_passed'] and report['certificate_sha256']==digest
 zero={alg.tup(r['key'])for r in report['zero_records']}
 with gzip.open(SRC/'audit/symbolic_columns_independent.json.gz','rt')as f:records=json.load(f)['rows']
 expected={alg.tup(r['key']):r['full_views']for r in records};indices={alg.tup(r['key']):i for i,r in enumerate(records)}
 cache={};seen=Counter();matches=[];n=0
 for view in ind.enumerate_inputs():
  n+=1
  if view not in cache:cache[view]=ind.symbolic(*view)
  key=cache[view];seen[key]+=1;iszero=key in zero;pred=predicate(view)
  if iszero!=pred:raise RuntimeError('Structural zero iff failure at '+repr(view))
  if pred:
   d,deg,edges,blocks=view;matches.append(dict(independent_column=indices[key],d=d,neighbor_degrees=deg,neighbor_edges=sorted(edges),blocks=[dict(T=T,degree=h,x=x,b=b)for T,h,x,b in blocks]))
 assert n==231816 and len(cache)==52609 and len(seen)==9867 and dict(seen)==expected
 assert len(matches)==7 and len(zero)==6
 out=dict(all_passed=True,certificate='subcubic_certificate.json',certificate_sha256=digest,independently_regenerated_full_views=n,independent_combinatorial_records=len(cache),independent_symbolic_columns=len(seen),all_symbolic_tuples_and_multiplicities_matched=True,structural_iff_tested_on_every_full_view=True,zero_columns=sorted(indices[k]for k in zero),zero_full_views=len(matches),zero_views=matches,optimizer_used=False,elapsed_seconds=time.time()-start)
 dest=ROOT/'audit/final_average_1000/zero_structure_full_replay.json';dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='zero_views'},indent=2))
if __name__=='__main__':main()
