#!/usr/bin/env python3
"""Independent all-q certificate replay from 13 elementary cleared terms.
Uses only standard-library integer sparse polynomial arithmetic. No SymPy,
primary-code import, graph/parameter search, or numerical evaluation.
"""
from pathlib import Path
from math import comb
import hashlib,json
if not __debug__: raise RuntimeError('Optimized Python is not accepted.')
HERE=Path(__file__).resolve().parent
PRIMARY=HERE.parent/'upper_interval'

def plus(*polys):
    z={}
    for p in polys:
        for mon,c in p.items(): z[mon]=z.get(mon,0)+c
    return {mon:c for mon,c in z.items() if c}

def times(*polys):
    z={(0,0):1}
    for p in polys:
        n={}
        for (i,j),a in z.items():
            for (k,l),b in p.items():
                mon=(i+k,j+l); n[mon]=n.get(mon,0)+a*b
        z={mon:c for mon,c in n.items() if c}
    return z

def scale(p,c): return {mon:c*a for mon,a in p.items() if c*a}
def power(p,k): return times(*([p]*k))
def require(test,message):
    if not test: raise RuntimeError(message)

one={(0,0):1}; q={(1,0):1}; x={(0,1):1}
r=plus(one,x)
a=plus(scale(one,2),scale(q,-1))
V=plus(one,scale(q,2),scale(power(q,2),3))
W=plus(V,one); Vm1=plus(V,scale(one,-1))
S3=plus(one,q,power(q,2))
A=times(plus(r,scale(q,-1)),S3)
G=plus(A,one); H=plus(A,scale(one,-1))
Cnum=plus(times(plus(x,scale(power(q,2),-2)),power(r,2)),
          scale(times(x,power(q,2)),-2))
X2=plus(x,scale(power(q,2),2))
# Positive denominator D=5040*q^2*r^2*a^7*V*W^7*G^3.
# Clear C_L * 2*(q/a)^odd/odd one logarithm-series term at a time.
Lterms=[]
for odd in [1,3,5,7]:
    require(2520%odd==0,'Nonintegral L clearing coefficient')
    Lterms.append(scale(times(Cnum,power(q,odd),power(a,7-odd),
                              V,power(W,7),power(G,3)),2520//odd))
# Clear the four entropy terms p, p^2/2, 2*p*c, 2*p*c^3/3 separately.
Eterms=[
 scale(times(power(q,3),r,power(a,7),V,power(W,7),power(G,3)),5040),
 scale(times(power(q,4),power(a,7),V,power(W,7),power(G,3)),2520),
 scale(times(power(q,3),r,power(a,7),V,power(W,7),H,power(G,2)),10080),
 scale(times(power(q,3),r,power(a,7),V,power(W,7),power(H,3)),3360)]
# Clear -(d+2)/8 times each finite upper log(V) term and then its tail.
Vterms=[]
for odd in [1,3,5,7]:
    require(1260%odd==0,'Nonintegral V clearing coefficient')
    Vterms.append(scale(times(X2,power(r,2),power(a,7),V,
                              power(W,7-odd),power(Vm1,odd),power(G,3)), -1260//odd))
Vterms.append(scale(times(X2,power(r,2),power(a,7),power(Vm1,9),power(G,3)),-35))
P=plus(*Lterms,*Eterms,*Vterms)
# Apply x=3*q^2+s exactly by the binomial theorem.
shift={}
for (i,j),c in P.items():
    for sdegree in range(j+1):
        mon=(i+2*(j-sdegree),sdegree)
        shift[mon]=shift.get(mon,0)+c*comb(j,sdegree)*3**(j-sdegree)
shift={mon:c for mon,c in shift.items() if c}
k=max(i for i,j in shift)
# q=(1+2u)/(2(1+u)); clear [2(1+u)]^k, independently expand both binomials.
M={}
for (i,j),c in shift.items():
    for e in range(i+1):
        left=c*comb(i,e)*2**(e+k-i)
        for f in range(k-i+1):
            mon=(e+f,j)
            M[mon]=M.get(mon,0)+left*comb(k-i,f)
M={mon:c for mon,c in M.items() if c}
primary=json.loads((PRIMARY/'CERTIFICATE.json').read_text())
expected={(t['u_degree'],t['s_degree']):int(t['coefficient']) for t in primary['terms']}
require(len(expected)==len(primary['terms']),'Duplicate primary coefficient indices')
require(expected==M,'Transformed coefficient mismatch')
require(primary['clearing_power']==k,'Clearing-power mismatch')
numdegrees=[max(i for i,j in P),max(j for i,j in P)]
Mdegrees=[max(i for i,j in M),max(j for i,j in M)]
require(primary['numerator_degrees']==numdegrees,'Numerator-degree mismatch')
require(primary['transformed_degrees']==Mdegrees,'Transformed-degree mismatch')
negative=[dict(u_degree=i,s_degree=j,coefficient=str(c)) for (i,j),c in sorted(M.items()) if c<0]
passed=not negative and M.get((0,0),0)>0
require(primary['negative_count']==len(negative),'Negative-count mismatch')
require(int(primary['strict_constant'])==M.get((0,0),0),'Constant mismatch')
require(primary['passed']==passed,'Certificate-status mismatch')
terms=[dict(u_degree=i,s_degree=j,coefficient=str(c)) for (i,j),c in sorted(M.items())]
canonical=json.dumps(terms,sort_keys=True,separators=(',',':'))
out={'status':'PASS' if passed else 'MATCHED_BUT_CERTIFICATE_NOT_POSITIVE',
 'method':'Independent standard-library sparse integer expansion from 13 separately cleared rational terms',
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'primary_input_sha256':hashlib.sha256((PRIMARY/'CERTIFICATE.json').read_bytes()).hexdigest(),
 'numerator_degrees':numdegrees,'numerator_nonzero_terms':len(P),
 'clearing_power':k,'transformed_degrees':Mdegrees,
 'nonzero_terms':len(M),'negative_terms':negative,'strict_constant':str(M.get((0,0),0)),
 'all_transformed_coefficients_match':True,'all_nonzero_coefficients_strictly_positive':all(c>0 for c in M.values()),
 'coefficient_certificate_passed':passed,
 'canonical_transformed_sha256':hashlib.sha256(canonical.encode()).hexdigest(),
 'scope':'All real d>=3 and 1/2<=q<1, subject to the separately checked analytic bounds. q=1 excluded.',
 'terms':terms,
 'numerator_terms':[dict(q_degree=i,x_degree=j,coefficient=str(c)) for (i,j),c in sorted(P.items())]}
(HERE/'INDEPENDENT_EXACT_REPLAY.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({key:val for key,val in out.items() if key not in ['terms','numerator_terms']},indent=2))
