#!/usr/bin/env python3
"""Optimizer-free exact interval algebra and replay; no numerical-library imports."""
from __future__ import annotations
if not __debug__:
    raise RuntimeError('Exact interval algebra refuses Python -O/-OO')
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from pathlib import Path
import hashlib,json
import symbolic_engine as s
STAGES={'global_d2':(2,16),'global_d4':(4,20)}
LIMITS=(10**6,10**9,10**12)
class CheckFailure(RuntimeError): pass
class AffineInconsistent(RuntimeError): pass
def require(ok,message):
    if not ok: raise CheckFailure(message)
def dump(path,value):
    Path(path).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return tuple(p or [0])
def add(a,b):
    out=[F(0)]*max(len(a),len(b))
    for j,v in enumerate(a):out[j]+=v
    for j,v in enumerate(b):out[j]+=v
    return trim(out)
def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:out[i+j]+=x*y
    return trim(out)
@lru_cache(None)
def power_to_bernstein(poly,degree):
    require(len(poly)<=degree+1,'Polynomial exceeds declared Bernstein degree')
    return tuple(sum((F(v)*F(comb(i,k),comb(degree,k)) for k,v in enumerate(poly) if v and k<=i),F(0)) for i in range(degree+1))
@lru_cache(None)
def basis_power(i,d):
    return tuple([0]*i+[comb(d,i)*comb(d-i,h)*(-1)**h for h in range(d-i+1)])
def bernstein_to_power(values):
    out=(F(0),);d=len(values)-1
    for i,v in enumerate(values):
        if v:out=add(out,tuple(v*c for c in basis_power(i,d)))
    return out
class AffineSystem:
    """Deterministic sparse exact echelon elimination; inconsistency has one type."""
    def __init__(self,n):self.n=n;self.basis={};self.rows=0
    def add(self,entries,rhs):
        self.rows+=1
        row={int(j):F(x) for j,x in (entries.items() if isinstance(entries,dict) else enumerate(entries)) if x}
        require(all(0<=j<self.n for j in row),'Affine coordinate range')
        rhs=F(rhs)
        for pivot,(known,b) in sorted(self.basis.items()):
            if pivot not in row:continue
            q=row.pop(pivot);rhs-=q*b
            for j,x in known.items():
                if j==pivot:continue
                v=row.get(j,F(0))-q*x
                if v:row[j]=v
                else:row.pop(j,None)
        if not row:
            if rhs:raise AffineInconsistent(f'Exact row {self.rows}: zero LHS but RHS {rhs}')
            return
        pivot=min(row);q=row[pivot]
        self.basis[pivot]=({j:x/q for j,x in row.items()},rhs/q)
    def solve(self,free_values):
        require(len(free_values)==self.n,'Free approximation dimension')
        x=list(map(F,free_values))
        for pivot,(row,b) in sorted(self.basis.items(),reverse=True):
            x[pivot]=b-sum((a*x[j] for j,a in row.items() if j!=pivot),F(0))
        return x
    def affine_map(self):
        """Return x=particular+sum transform[:,free] * independent x_free."""
        free=[i for i in range(self.n) if i not in self.basis]
        index={v:i for i,v in enumerate(free)}
        base=[F(0)]*self.n;T=[({index[i]:F(1)} if i in index else {}) for i in range(self.n)]
        for p,(row,b) in sorted(self.basis.items(),reverse=True):
            val=b;tr={}
            for j,a in row.items():
                if j==p:continue
                val-=a*base[j]
                for k,c in T[j].items():
                    v=tr.get(k,F(0))-a*c
                    if v:tr[k]=v
                    else:tr.pop(k,None)
            base[p]=val;T[p]=tr
        return free,base,T
    def summary(self):return {'equations':self.rows,'variables':self.n,'rank':len(self.basis),'free_parameters':self.n-len(self.basis)}
class Problem:
    def __init__(self,stage,data):
        require(stage in STAGES,'Unknown stage')
        self.stage=stage;self.d,self.M=STAGES[stage];self.R=self.M-self.d
        self.data=data;self.rows=data['rows'];self.m=224;self.nv=self.m*(self.d+1)
        require(len(s.SUPPORT)==self.m,'Full support dimension')
        require(data['generation_sha256']==s.GEN_SHA,'Generation digest')
        self.cleared=[s.cleared(row) for row in self.rows]
        for target,features in self.cleared:
            require(len(features)==self.m,'Full cleared feature count')
            require(len(target)<=self.M+1 and all(len(p)<=self.R+1 for p in features),'Fixed slack degree insufficient')
        self.c4=list(data['c4_columns']) if 'c4_columns' in data else list(s.c4_column_ids(data))
        require(len(self.c4)==2 and len(set(self.c4))==2,'Exactly two C4 exposures')
        self.sign_indices=list(range(170*(self.d+1),self.nv))
    def exact_row(self,c,r):
        require(0<=c<len(self.rows) and 0<=r<=self.M,'Bernstein row range')
        target,features=self.cleared[c];row={}
        for j,poly in enumerate(features):
            if not any(poly):continue
            bb=power_to_bernstein(tuple(poly),self.R)
            for i in range(self.d+1):
                k=r-i
                if 0<=k<=self.R and bb[k]:
                    row[j*(self.d+1)+i]=bb[k]*F(comb(self.d,i)*comb(self.R,k),comb(self.M,r))
        return row,power_to_bernstein(tuple(target),self.M)[r]
    def power_row(self,c,k):
        target,features=self.cleared[c];row={}
        for j,poly in enumerate(features):
            if not any(poly):continue
            for i in range(self.d+1):
                basis=basis_power(i,self.d)
                val=sum((F(a)*poly[k-h] for h,a in enumerate(basis) if a and 0<=k-h<len(poly)),F(0))
                if val:row[j*(self.d+1)+i]=val
        return row,F(target[k] if k<len(target) else 0)
    def calibration(self):
        system=AffineSystem(self.nv);source=[]
        for c in self.c4:
            for k in range(self.M+1):
                row,rhs=self.power_row(c,k);system.add(row,rhs)
                source.append({'column':c,'power_coefficient':k,'entries':{str(j):str(v) for j,v in row.items()},'rhs':str(rhs)})
        return system,source
    def certificate(self,x):
        require(len(x)==self.nv,'Certificate coordinate dimension')
        return {'schema_version':'degree_four_weighted_interval_v1','stage':self.stage,'interval':['0','1'],
                'degree':self.d,'bernstein_degree':self.M,'generation_sha256':s.GEN_SHA,'support':list(s.SUPPORT),
                'clearing_factor':'12*A*Z_C','bernstein_polynomials':[[str(v) for v in x[j*(self.d+1):(j+1)*(self.d+1)]] for j in range(self.m)]}
    def parse(self,cert):
        require(cert['schema_version']=='degree_four_weighted_interval_v1','Certificate schema')
        require(cert['stage']==self.stage and cert['interval']==['0','1'],'Certificate stage/interval')
        require(cert['degree']==self.d and cert['bernstein_degree']==self.M,'Certificate degrees')
        require(cert['generation_sha256']==s.GEN_SHA and cert['support']==list(s.SUPPORT),'Certificate coverage/support')
        require(cert['clearing_factor']=='12*A*Z_C','Certificate denominator normalization')
        polynomials=cert['bernstein_polynomials']
        require(len(polynomials)==self.m and all(len(p)==self.d+1 for p in polynomials),'Coefficient shape')
        return [F(v) for p in polynomials for v in p]
    def verify(self,cert):
        x=self.parse(cert);failures=[];slack_hash=hashlib.sha256();minimum=None;tight=[]
        for j in self.sign_indices:
            if x[j]<0:failures.append({'kind':'negative_multiplier','coordinate':j,'value':str(x[j])})
        ps=[bernstein_to_power(x[j*(self.d+1):(j+1)*(self.d+1)]) for j in range(self.m)]
        for c,(target,features) in enumerate(self.cleared):
            slack=tuple(map(F,target))
            for p,q in zip(ps,features):
                if any(p) and any(q):slack=add(slack,tuple(-v for v in mul(p,q)))
            bern=power_to_bernstein(tuple(slack),self.M)
            for k,v in enumerate(bern):
                if v<0:failures.append({'kind':'negative_slack','column':c,'coefficient':k,'value':str(v)})
                if minimum is None or v<minimum:minimum=v
            if all(v==0 for v in bern):tight.append(c)
            if c in self.c4 and any(bern):failures.append({'kind':'nonzero_c4','column':c,'power_coefficients':list(map(str,slack))})
            slack_hash.update(('|'.join(map(str,bern))+'\n').encode())
        return {'all_passed':not failures,'failures':failures,'minimum_bernstein':str(minimum),
                'checked_multiplier_coefficients':len(self.sign_indices),'checked_slack_coefficients':len(self.rows)*(self.M+1),
                'columns':len(self.rows),'tight_columns':tight,'slack_certificate_sha256':slack_hash.hexdigest(),
                'maximum_absolute_bernstein_coefficient':str(max(map(abs,x),default=F(0))),
                'maximum_denominator':max((v.denominator for v in x),default=1)}
