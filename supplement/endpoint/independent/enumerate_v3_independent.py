#!/usr/bin/env python3
"""Independent exhaustive V3 common-denominator columns over full degree views."""
import sys
if not __debug__: raise RuntimeError('Optimized Python is forbidden for V3 exact checks')
from itertools import product,combinations
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
import json,gzip,time
from base_views_v1_copy import subsets
from enumerate_v2_independent import partitions
from features_v3_independent import feature_vector,canonical,SIGS,CLASSES,PATTERNS
HERE=Path(__file__).resolve().parent

def run():
    start=time.time();cache={};mult=Counter();full=0;bydegree=Counter();degree3={}
    for d in (2,3):
        for E in subsets(combinations(range(1,d+1),2)):
            for deg in product((2,3),repeat=d):
                amounts=[deg[u-1]-1-sum(u in e for e in E) for u in range(1,d+1)]
                if min(amounts)<0:continue
                stubs=tuple((u,j) for u,a in enumerate(amounts,1) for j in range(a))
                for P in partitions(stubs):
                    Ts=tuple(tuple(u for u,j in B) for B in P)
                    opts=[[(h,x,b) for h in (2,3) if h>=len(T) for x,b in ((0,0),(1,0),(1,1)) if h!=len(T) or x==b] for T in Ts]
                    for choice in product(*opts):
                        blocks=tuple(sorted((T,*c) for T,c in zip(Ts,choice)));K=(d,deg,E,blocks)
                        if K not in cache:cache[K]=canonical(feature_vector(*K))
                        col=cache[K];mult[col]+=1;full+=1;bydegree[d]+=1
                        if col in degree3:assert degree3[col]==int(d==3)
                        degree3[col]=int(d==3)
    assert full==231816
    actual=json.loads((HERE/'graph_views_v3.json').read_text())
    for A in actual:
        col=canonical(tuple(map(F,A['features'])));assert col in mult
        assert degree3[col]==int(A['view']['d']==3)
    data={'all_passed':True,'full_views':full,'full_views_by_degree':dict(bydegree),
          'distinct_combinatorial_views':len(cache),'columns':len(mult),'actual_exposures_covered':len(actual),
          'schema':'c,53 V2 zero features,60 star zero features,20 exterior inequalities,7 association inequalities',
          'rows':[{'denominator':col[-1],'numerators':list(col[:-1]),'full_views':n,'degree3':degree3[col]} for col,n in sorted(mult.items())]}
    with gzip.open(HERE/'independent_columns_v3.json.gz','wt') as f:json.dump(data,f,sort_keys=True)
    out={k:v for k,v in data.items() if k!='rows'};out['elapsed_seconds']=time.time()-start
    (HERE/'independent_enumeration_v3.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':run()
