#!/usr/bin/env python3
"""Complete symbolic-key/multiplicity comparison; no point-activity anchors."""
if not __debug__:raise RuntimeError('Complete key comparison refuses optimized Python')
from pathlib import Path
from collections import Counter
import json,gzip
import compare_symbolic as c
import symbolic_features_independent as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
a=c.load(HERE/'symbolic_columns_independent.json.gz')['rows'];b=c.load(ROOT/'symbolic_columns.json.gz')
ic=Counter({c.tuplify(r['key']):r['full_views']for r in a});pc=Counter({c.primary_key(r):r['full_views']for r in b['rows']})
c.require(len(pc)==len(b['rows']) and ic==pc,'Complete key/multiplicity mismatch')
c.require(sum(ic.values())==231816 and len(ic)==9867,'Wrong complete universe')
c.require(c.tuplify(s.independent_association())==c.tuplify(b['association']),'Association polynomials mismatch')
out=dict(all_passed=True,full_views=231816,complete_keys=9867,sampled_activity_anchors_used=False)
(HERE/'complete_key_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
