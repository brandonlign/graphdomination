#!/usr/bin/env python3
"""Exact symbolic unit tests for verifier's affine and Bernstein integer scaling."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Independent Bernstein unit tests refuse optimized Python')
from fractions import Fraction as F
from math import comb,lcm
from pathlib import Path
import sympy as sp,json
HERE=Path(__file__).resolve().parent

def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def run():
    t,z=sp.symbols('t z');base=10;count=0
    stages=[(F(0),F(1),2,4),(F(0),F(1),4,8)]+[(a,b,6,12)for a,b in((F(0),F(1)),(F(0),F(1,2)),(F(1,2),F(1)),(F(1,4),F(1,2)),(F(1,2),F(3,4)),(F(3,4),F(1)))]
    rational=lambda x:sp.Rational(x.numerator,x.denominator)
    for a,b,d,e in stages:
        beta=[F((-1)**i*(i+2),i+3)for i in range(d+1)];H=lcm(*(x.denominator for x in beta));bi=[int(x*H)for x in beta]
        mono=[sum(bi[i]*comb(d,i)*comb(d-i,r-i)*(-1)**(r-i)for i in range(r+1))for r in range(d+1)]
        literal=H*sum(rational(x)*comb(d,i)*t**i*(1-t)**(d-i)for i,x in enumerate(beta))
        require(sp.expand(literal-sum(x*t**i for i,x in enumerate(mono)))==0,'Bernstein polynomial expansion failed');count+=1
        poly=[(-1)**j*(j+1)for j in range(base+1)];B=lcm(a.denominator,b.denominator);A=int(a*B);D=int((b-a)*B)
        transform=[sum(v*comb(j,k)*A**(j-k)*D**k*B**(base-j)for j,v in enumerate(poly)if j>=k)for k in range(base+1)]
        literal=B**base*sum(x*(rational(a)+(rational(b)-rational(a))*t)**i for i,x in enumerate(poly))
        require(sp.expand(literal-sum(x*t**i for i,x in enumerate(transform)))==0,'Affine scaling failed');count+=1
        M=base+d+e;val=[(-1)**j*(j+2)for j in range(base+d+1)];div=[comb(M,j)for j in range(len(val))];Q=lcm(*div)
        bern=[sum(x*comb(k,j)*(Q//div[j])for j,x in enumerate(val)if j<=k)for k in range(M+1)]
        literal=sum(x*comb(M,k)*t**k*(1-t)**(M-k)for k,x in enumerate(bern))
        require(sp.expand(literal-Q*sum(x*t**j for j,x in enumerate(val)))==0,'Bernstein recovery failed');count+=1
    out=dict(all_passed=True,exact_symbolic_identities=count,stage_combinations=len(stages),method='SymPy expansion, no activity evaluation or floating-point test')
    (HERE/'bernstein_unit_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':run()
