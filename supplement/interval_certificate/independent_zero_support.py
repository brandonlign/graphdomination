#!/usr/bin/env python3
"""Replacement-certificate equality premise from independent exact polynomials."""
if not __debug__:raise RuntimeError('Independent zero support refuses optimized Python')
import sys,json,gzip,hashlib,math,argparse
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'inputs/upstream'
sys.modules['scipy']=None;sys.modules['scipy.optimize']=None
sys.path.insert(0,str(SOURCE/'audit'))
import symbolic_features_independent as ind
import verify_interval_independent as alg

def verify(candidate,sign_report,output):
 cp=Path(candidate);digest=hashlib.sha256(cp.read_bytes()).hexdigest();proof=json.loads(Path(sign_report).read_text())
 if not proof['all_passed']or proof['certificate_sha256']!=digest:raise RuntimeError('Matching independently verified signs required')
 cert=json.loads(cp.read_text());names=cert['support'];assoc=ind.independent_association();beta=[[F(x)for x in row]for row in cert['bernstein_polynomials']]
 H=math.lcm(*(x.denominator for row in beta for x in row));P=[tuple(int(H*x)for x in (a,2*(b-a),a-2*b+c))for a,b,c in beta]
 zero=[];cubic_nonzero=0;views=0
 with gzip.open(SOURCE/'audit/symbolic_columns_independent.json.gz','rt')as f:rows=json.load(f)['rows']
 for i,row in enumerate(rows):
  K=alg.tup(row['key']);tar,fs=alg.cleared(K,names,assoc);slack=alg.scale(tar,H)
  for p,q in zip(P,fs):slack=alg.sub(slack,alg.times(p,q))
  if not any(slack):zero.append(dict(independent_column=i,root_degree=K[0],full_views=row['full_views'],key=K))
  elif K[0]==3:cubic_nonzero+=1
  views+=row['full_views']
 if len(rows)!=9867 or views!=231816:raise RuntimeError('Universe coverage mismatch')
 degrees=sorted({r['root_degree']for r in zero});strict=all(r['root_degree']==2 for r in zero)
 out=dict(all_passed=True,certificate_sha256=digest,independent_columns=len(rows),full_views=views,zero_columns=len(zero),zero_full_views=sum(r['full_views']for r in zero),zero_root_degrees=degrees,nonzero_cubic_columns=cubic_nonzero,all_cubic_slacks_strict_on_open_interval=strict,zero_records=zero,optimizer_used=False,equality_consequence='Interior equality forces 2-regularity; existing strict cycle theorem and independent z=1 equality theorem then give exactly C4 unions.'if strict else'No cubic-strictness equality conclusion from this candidate; original independent equality theorem remains available.')
 Path(output).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='zero_records'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('certificate');p.add_argument('sign_report');p.add_argument('output');a=p.parse_args();verify(a.certificate,a.sign_report,a.output)
