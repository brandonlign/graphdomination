#!/usr/bin/env python3
if not __debug__:raise RuntimeError('Structure check refuses Python -O/-OO')
import gzip,json
from pathlib import Path
B=Path(__file__).resolve().parent
v=[json.loads(line) for line in gzip.open(B/'inputs/views.jsonl.gz','rt')]
z=json.loads((B/'ZERO_SUPPORT.json').read_text())
def zero_shape(a):
 if a['d']!=2 or a['degrees']!=[2,2] or a['edges']:return False
 if a['blocks']==[[1,2]]:
  return a['choices'] in ([[2,0,0]],[[2,1,1]])
 return a['blocks']==[[1],[2]] and all(h in (2,3,4) and x==b==0 for h,x,b in a['choices'])
actual=[i for i,a in enumerate(v) if zero_shape(a)]
if actual!=z['raw_indices']:raise RuntimeError('Structural zero predicate differs from arithmetic')
empty4=[i for i,a in enumerate(v) if a['d']==4 and all(x==0 and b==0 for h,x,b in a['choices'])]
r=json.loads((B/'ENDPOINT_EMPTY_DEGREE4.json').read_text())
if empty4!=[x['raw_index'] for x in r['views']]:raise RuntimeError('Endpoint empty coverage differs')
empty_boundary=[i for i in empty4 if not v[i]['blocks']]
if len(empty_boundary)!=3:raise RuntimeError('Perfect matching empty-boundary views omitted')
out=dict(all_passed=True,existing_full_views_tested=len(v),no_new_enumeration=True,zero_predicate_indices=actual,
         degree4_empty_views=len(empty4),degree4_empty_boundary_indices=empty_boundary,
         description='Two shared-boundary C4 views; nine ordered separate-boundary all-unselected views')
(B/'STRUCTURAL_CHECKS.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps(out,indent=2))
