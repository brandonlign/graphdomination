#!/usr/bin/env python3
"""Compare independently derived entire cleared polynomials, not sampled signs."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Independent clearing comparison refuses optimized Python')
import sys,json,gzip,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
import compare_symbolic as c
import verify_interval_independent as v
import symbolic_features_independent as s
sys.path.insert(0,str(ROOT))
import symbolic_engine as primary

def run():
    start=time.time();data=primary.load();assoc=s.independent_association();maxdeg=0
    for row in data['rows']:
        K=c.primary_key(row);left=primary.cleared(row);right=v.cleared(K,data['support'],assoc)
        c.require(left==right,'Entire cleared-polynomial discrepancy at '+str(row['id']))
        maxdeg=max(maxdeg,max(len(x)-1 for x in [right[0]]+right[1]))
    report=dict(all_passed=True,complete_cleared_polynomials_compared=len(data['rows'])*(len(data['support'])+1),columns=len(data['rows']),maximum_degree=maxdeg,
                method='Independent formulas compared coefficientwise to primary normalization; no activity sampling.',elapsed_seconds=time.time()-start)
    (HERE/'clearing_comparison.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':run()
