#!/usr/bin/env python3
"""Small exact analytic controls; no optimizer, source imports, or graph census."""
if not __debug__:
    raise RuntimeError('Exact checks refuse optimized Python')
from pathlib import Path
import gzip, hashlib, json
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent

def require(ok, text):
    if not ok: raise RuntimeError(text)
def norm(p):
    q=list(p)
    while len(q)>1 and q[-1]==0:q.pop()
    return tuple(q or [0])
def add(a,b):return norm([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def neg(a):return tuple(-x for x in a)
def scale(a,n):return norm([n*x for x in a])
def mul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):r[i+j]+=x*y
    return norm(r)
def deriv(a):return norm([i*a[i] for i in range(1,len(a))])
A=(1,4,6); B=(1,3); z=(0,1); z2=(0,0,1)
Z={1:(1,),2:(1,2),3:(1,3,3)}
for n in range(4,7):Z[n]=add(add(Z[n-1],mul(z,Z[n-2])),mul(z2,Z[n-3]))
H={n:add(mul(A,deriv(Z[n])),neg(scale(mul(B,Z[n]),n))) for n in range(3,7)}
expected={3:(0,0,6,9),4:(0,),5:(0,0,0,5,15),6:(0,0,0,6,30,18)}
require(H==expected,'Cycle base identities fail')
# Linear forms in three unrestricted symbols Z_(n-6), Z_(n-5), Z_(n-4).
# Thus the following identity is symbolic, not a finite index check.
def ladd(a,b):return tuple(add(x,y) for x,y in zip(a,b))
def lmul(p,a):return tuple(mul(p,x) for x in a)
U={0:((1,),(0,),(0,)),1:((0,),(1,),(0,)),2:((0,),(0,),(1,))}
for j in range(3,6):U[j]=ladd(ladd(U[j-1],lmul(z,U[j-2])),lmul(z2,U[j-3]))
T=ladd(ladd(lmul((-1,-3),U[5]),lmul((1,2),U[4])),lmul((0,2,5,3),U[3]))
force=ladd(lmul((0,0,0,1,3),U[1]),lmul((0,0,0,0,2,3),U[0]))
require(T==force,'General positive forcing identity fails')
# This is a reader-only cross-check, not independent regeneration of the LP.
p=ROOT/'interval'
d=json.load(gzip.open(p/'symbolic_columns.json.gz','rt'))
r=json.loads((p/'primary_certificate_replay.json').read_text())
tight=r['exact']['tight_columns']
full=[]
for j in tight:
    row=d['rows'][j]
    require(row['d']==2 and row['representative']['d']==2,'Reported zero support has a high root')
    require(row['full_views']==len(row['source_indices']),'Reported source multiplicity mismatch')
    full.extend(row['source_indices'])
require(len(set(full))==len(full),'Duplicate reported full-view support')
require(len(tight)==8 and len(full)==11,'Unexpected reported zero support size')
report={'all_passed':True,'cycle_base_polynomials':{str(n):list(H[n]) for n in H},'symbolic_forcing_linear_form':[list(x) for x in T], 'reported_tight_columns':tight,'reported_full_view_indices':sorted(full),'reported_full_view_count':len(full),'reported_root_degrees':[2], 'reader_only_limit':'Zero-support data were read, not independently reconstructed in this analytic audit; the full-view equality checker supplies independent arithmetic/full-view reconstruction. No source code was imported.'}
(HERE/'analytic_checks.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
