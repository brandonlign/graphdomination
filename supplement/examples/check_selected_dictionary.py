#!/usr/bin/env python3
"""Read-only selected low-bit inspection; no imported proof engine or optimizer.
Usage: python check_selected_dictionary.py /path/to/extracted/exact/supplement
Outputs only in this checker directory; all rationals are parsed exactly.
"""
from pathlib import Path
from itertools import combinations_with_replacement
from fractions import Fraction as F
from collections import Counter
import argparse,gzip,hashlib,json
from math import lcm
ROOT=Path(__file__).resolve().parent
CERT_SHA='c6f20f9d56dd89013c8e1b186e1f0614368b7585c3243fa12ddc8b64fbfa0205'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def dump(name,x):(ROOT/name).write_text(json.dumps(x,indent=2)+'\n')
p=argparse.ArgumentParser();p.add_argument('supplement',type=Path);args=p.parse_args();base=args.supplement
folder=base/'interval_certificate/inputs/upstream';certpath=base/'interval_certificate/subcubic_certificate.json'
need(hashlib.sha256(certpath.read_bytes()).hexdigest()==CERT_SHA,'Not the selected low-bit certificate')
c=json.loads(certpath.read_text());support=c['support'];beta=[[F(x)for x in row]for row in c['bernstein_polynomials']]
need(c['degree']==2 and c['bernstein_degree']==16 and len(support)==74 and all(len(row)==3 for row in beta),'selected degree/shape')
rev=[(t,a,b,d)for t in(0,1)for a in range(4)for b in range(4)for d in range(a+1,4)]
stars=[(d,x,T)for d in(2,3)for x in(0,1)for T in combinations_with_replacement(range(4),d)]
classes=[(d,h,k)for d in(2,3)for h in(2,3)for k in range(1,h+1)]
patterns=[(2,r)for r in range(3)]+[(3,r)for r in range(4)]
names=['normalization','m2','m3','q2','q3','j']+[f'R_{t}_{a}_{b}_{d}'for t,a,b,d in rev]+[f'S_{d}_{x}_'+''.join(map(str,T))for d,x,T in stars]+[f'{typ}_{d}_{h}_{k}'for d,h,k in classes for typ in('X','B')]+[f'L_{d}_{r}'for d,r in patterns]
need(len(names)==141 and len(set(names))==141,'dictionary length')
source=json.loads((folder/'inputs/dependencies/enumeration_v3.json').read_text())
need(names==source['equality_names']+source['inequality_names'],'dictionary does not match the coordinate order')
mult=dict(zip(support,beta));rows=[];counts=Counter()
for i,name in enumerate(names):
 family='objective'if i==0 else 'basic_transport'if i<6 else 'reversal_transport'if i<54 else 'star_transport'if i<114 else 'exterior'if i<134 else 'association'
 selected=name in mult;active=selected and any(mult[name]);counts[family+'_available']+=1
 counts[family+'_selected']+=selected;counts[family+'_active']+=active
 row=dict(feature_index=i,code_name=name,family=family,selected=selected,active_nonzero_polynomial=bool(active),multiplier_bernstein=list(map(str,mult[name]))if selected else None)
 if i==0:row['warning']='Stored features6[0]=6 O_C; objective is O_C/Z_C. This file key is NOT the constant normalization feature. The LP normalization is separately sum(pi)=1; its prescribed dual right side is c4.'
 if 6<=i<54:row['indices']=list(rev[i-6])
 if 54<=i<114:
  d,x,T=stars[i-54];row['indices']=[d,x,list(T)]
 if 114<=i<134:row['indices']=list(classes[(i-114)//2]);row['kind']='q N-X'if (i-114)%2==0 else 'q^(h-k) X-B'
 if i>=134:row['indices']=list(patterns[i-134])
 rows.append(row)
active_ineq=[r['code_name']for r in rows if r['feature_index']>=114 and r['active_nonzero_polynomial']]
selected_ineq=[r['code_name']for r in rows if r['feature_index']>=114 and r['selected']]
zero_selected=[name for name in support if not any(mult[name])]
need(set(active_ineq)=={'B_3_3_1','L_2_1','L_2_2','L_3_0','L_3_1','L_3_2','L_3_3'},'active inequality set')
need(set(zero_selected)=={'X_3_2_1','X_3_3_1'},'selected zero set')
need(len([r for r in rows[1:114]if r['selected']])==65 and len([r for r in rows[1:114]if r['active_nonzero_polynomial']])==65,'65 selected active transports')
need(len(selected_ineq)==9 and all(v>=0 for name in selected_ineq for v in mult[name]),'inequality signs')
need(sum(v!=0 for row in beta for v in row)==205,'nonzero rational coefficient total')
with gzip.open(folder/'symbolic_columns.json.gz','rt')as f:data=json.load(f)
need(len(data['rows'])==9867 and sum(r['full_views']for r in data['rows'])==231816,'symbolic data counts')
need(all(len(r['features6'])==114 and len(r['exterior'])==10 for r in data['rows']),'dictionary shape')
need(len(data['association'])==7,'association ordering length')
nonzero=[v for row in beta for v in row if v]
metrics=dict(nonzero_coefficients=len(nonzero),total_integer_bits=sum(abs(v.numerator).bit_length()+v.denominator.bit_length()for v in nonzero),maximum_integer_bits=max(max(abs(v.numerator).bit_length(),v.denominator.bit_length())for v in nonzero),maximum_decimal_digits=max(max(len(str(abs(v.numerator))),len(str(v.denominator)))for v in nonzero),common_denominator_bits=lcm(*(v.denominator for v in nonzero)).bit_length(),certificate_json_bytes=certpath.stat().st_size)
need(metrics==dict(nonzero_coefficients=205,total_integer_bits=9389,maximum_integer_bits=240,maximum_decimal_digits=73,common_denominator_bits=539,certificate_json_bytes=7827),'selected exact encoding metrics')
dump('FEATURE_DICTIONARY.json',dict(selected_certificate_sha256=CERT_SHA,coordinate_count=len(rows),zero_based_index=True,rows=rows))
dump('SELECTED_SUPPORT_CHECKS.json',dict(all_passed=True,read_only_selected_certificate=True,not_a_new_certificate_or_full_replay=True,certificate_sha256=CERT_SHA,encoding_metrics=metrics,counts=dict(counts),active_inequalities=active_ineq,selected_inequalities=selected_ineq,selected_zero_polynomials=zero_selected,selected_polynomials=74,active_polynomials=72,nonzero_rational_coefficients=205,slack_degree=16,multiplier_degree=2,columns=9867,full_views=231816,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
print(json.dumps(dict(counts=counts,active_inequalities=active_ineq,all_passed=True),indent=2))
