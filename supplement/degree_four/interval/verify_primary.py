#!/usr/bin/env python3
"""Full-view polynomial certificate replay with no optimizer imports."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Primary certificate replay refuses Python -O/-OO')
from pathlib import Path
from collections import defaultdict
import sys,json,hashlib,gzip
import symbolic_engine as s
from solver_math import Problem

def full_views(data):
 s.require(data['generation_sha256']==s.GEN_SHA and data['support']==s.SUPPORT,'Full replay metadata')
 s.require(data['equality_names']==s.EQNAMES and data['inequality_names']==s.GNAMES,'Full row ordering')
 s.require(data['association']==json.loads(json.dumps(s.ASSOCIATION)),'Association denominator metadata')
 seen=defaultdict(list);sha=hashlib.sha256();degree_counts=defaultdict(int)
 for i,v in enumerate(s.grammar.all_views()):
  sha.update((s.canonical(v)+'\n').encode());cid=data['full_to_column'][i]
  s.require(0<=cid<len(data['rows']),'Mapping range');r=data['rows'][cid]
  s.require(r['id']==cid and s.key(s.feature(v))==s.key(r),f'Full polynomial feature mismatch at {i}')
  seen[cid].append(i);degree_counts[v['d']]+=1
 s.require(i+1==17553 and sha.hexdigest()==s.GEN_SHA and len(data['full_to_column'])==17553,'Full grammar universe')
 s.require(len(seen)==len(data['rows']),'All symbolic columns exposed')
 for r in data['rows']:s.require(r['source_indices']==seen[r['id']] and r['full_views']==len(seen[r['id']]),'Multiplicity disagreement')
 return dict(full_views=i+1,columns=len(data['rows']),degree_counts=dict(degree_counts),generation_sha256=sha.hexdigest())

def verify(path,full=True,data=None):
 if data is None:data=s.load()
 cert=json.loads(Path(path).read_text());report={}
 if full:report['full_regeneration']=full_views(data)
 p=Problem(cert['stage'],data);report['exact']=p.verify(cert)
 s.require(report['exact']['all_passed'],'Certificate exact Bernstein or calibration failure')
 report['all_passed']=True;return report
if __name__=='__main__':
 path=Path(sys.argv[1]);out=verify(path);dest=Path(sys.argv[2])if len(sys.argv)>2 else s.HERE/'primary_certificate_replay.json';s.dump(dest,out);print(json.dumps(out,sort_keys=True))
