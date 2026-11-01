#!/usr/bin/env python3
"""Exact consequences of read sources; no graph census, LP, or claim of novelty."""
from fractions import Fraction
from pathlib import Path
import json
if not __debug__: raise RuntimeError('Run without -O')
def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 while len(c)>1 and c[-1]==0:c.pop()
 return c
def neg(a):return [-x for x in a]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return add(c,[0])
def power(a,n):
 c=[1]
 for _ in range(n):c=mul(c,a)
 return c
def mono(n,c=1):return [0]*n+[c]
q=[0,1]
one=[1]
# A(z)(1-q)^4 = (1-q)^2(1+2q+3q^2), z=q/(1-q).
base=mul(power([1,-1],2),[1,2,3])
assert base==[1,0,0,-4,3]
# F_3^4 / A >= 1 for 0<q<=1/2.
d3=sub(power(sub(one,mono(4)),4),base)
d3_rhs=mul(mono(3),add(add([4,-7],mul(mono(5,2),[3,0,0,0,-2])),mono(13)))
assert d3==d3_rhs
# (F_2 F_4)^2 / A >= 1 for 0<q<=1/2.
pair=sub(mul(power(sub(one,mono(3)),2),power(sub(one,mono(5)),2)),base)
bracket=mul([1,-2],[2,1])
bracket=add(bracket,mono(3))
bracket=add(bracket,mul(mono(5,2),[2,0,0,-1,0,-1]))
bracket=add(bracket,mono(7));bracket=add(bracket,mono(13))
assert pair==mul(mono(3),bracket)
# Every displayed bracket is positive on 0<q<=1/2 by the written inequalities.
# Wagner Thm 4 expression uses factors of sizes s=2,...,5.
wagner={s:{'left':(2**(s-1)+1)**5,'right':17**s} for s in range(2,6)}
assert all(v['left']<=v['right'] for v in wagner.values())
assert 17**4<11**5
# Every edge-minimal profile: 2*n2 >= 3*n3+4*n4.
low_gap=11*4**4-7**4
cubic_block_gap=11**5*4**12*8**8-7**12*15**8
assert low_gap>0 and cubic_block_gap==9858121268852468831>0
# Explicit uncovered profile: once-subdivision of a 4-regular core.
profile_mean=(2*Fraction(5,7)+Fraction(18,31))/3
assert profile_mean==Fraction(436,651)>Fraction(7,11)
product_base=Fraction(7,4)**2*Fraction(31,16)
assert product_base==Fraction(1519,256)
gap=11**3*256**4-1519**4
assert gap==392686686655>0
out={'status':'PASS','degree3_product_numerator':d3,'degree2_degree4_pair_numerator':pair,
'wagner_factor_checks':wagner,'wagner_target_comparison':{'17^4':17**4,'11^5':11**5},
'all_edge_minimal_product_failure':{'degree2_gap':low_gap,'degree3_plus_1.5_degree2_gap':cubic_block_gap,'degree4_plus_2_degree2_gap':gap},
'once_subdivided_4regular':{'mean_bound':str(profile_mean),'product_base':str(product_base),'integer_gap':gap},
'limits':'Symbolic comparison checks only.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
