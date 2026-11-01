#!/usr/bin/env python3
"""Independent direct endpoint evaluation from the stored symbolic columns.
No imports of earlier evaluators; this reuses symbolic features, not grammar.
"""
if not __debug__:raise RuntimeError('Exact checks refuse optimized Python')
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import gzip, hashlib, json
R=Path(__file__).resolve().parents[2]/'interval'
O=Path(__file__).resolve().parent
c=json.loads((R/'certificate.json').read_text())
d=json.load(gzip.open(R/'symbolic_columns.json.gz','rt'))
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
require(c['support']==d['equality_names']+d['inequality_names'],'Support ordering differs')
P=dict(zip(c['support'],(Q(a[-1]) for a in c['bernstein_polynomials'])))
require(P['L_3_2,2,2']>0,'Degree-three endpoint multiplier differs')
require(P['L_4_2,2,2,2']==0,'Degree-four endpoint multiplier differs')
zeros=[];zero4=[];empty4=[];empty4full=0;fullzeros=0
for row in d['rows']:
    Z=sum(row['Z']);occ=Q(sum(row['O']),Z)
    slack=occ-Q(4,11)
    for name,F in zip(d['equality_names'],row['F']):
        slack-=P[name]*Q(sum(F),12*11*Z)
    for (root,h,k),(N,X,B) in zip(d['classes'],row['exterior']):
        slack-=Q(8,11)*P[f'EX_{root}_{h}_{k}']*(Q(N,2)-X)
        slack-=Q(8,11)*P[f'EB_{root}_{h}_{k}']*(Q(X,2**(h-k))-B)
    hs=row['representative']['degrees'];r=row['d'];M=r+sum(h-1 for h in hs)-min(h-1 for h in hs)
    T=Q(0)
    for bits in product((0,1),repeat=r):
        if all(bits):continue
        term=Q(1,2**r)
        for h,x in zip(hs,bits):
            if x:term*=1-Q(1,2**(h-1))
        T+=term
    N=T*2**M;D=2**M+N
    name='L_'+str(r)+'_'+','.join(map(str,sorted(hs)))
    slack-=Q(D,11)*P[name]*(occ-Q(N,D))
    H=12*11*Z*slack
    require(H>=0,'Negative endpoint slack')
    if H==0:
        zeros.append(row['id']);fullzeros+=row['full_views']
        if r==4:zero4.append(row['id'])
    if r==4 and all(X==0 for _,X,_ in row['exterior']):
        # Includes rows with no boundary vertices: all(empty) is deliberately True.
        require(H>0,'Degree-four empty-exterior support is endpoint-tight')
        empty4.append((row['id'],str(H)));empty4full+=row['full_views']
require(len(zeros)==80 and fullzeros==335 and len(zero4)==18,'Endpoint support counts differ')
report={'all_passed':True,'endpoint_zero_columns':zeros,'endpoint_zero_full_view_count':fullzeros,'endpoint_zero_degree_four_columns':zero4,'degree_four_all_x_zero_columns':empty4,'degree_four_all_x_zero_full_view_count':empty4full,'degree_three_P_at_one':str(P['L_3_2,2,2']),'degree_three_association_D_at_one':45,'degree_three_mu_at_one':str(Q(45,11)*P['L_3_2,2,2']),'empty_boundary_covered':True,'implementation_limit':'Independent direct un-cleared endpoint arithmetic, reusing the primary symbolic feature table; no grammar reconstruction here. The full-view equality checker supplies independently reconstructed full views.'}
(O/'endpoint_checks.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='degree_four_all_x_zero_columns'},indent=2,sort_keys=True))
print('All-x-zero degree-four column count:',len(empty4))
