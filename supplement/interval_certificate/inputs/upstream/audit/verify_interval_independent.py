#!/usr/bin/env python3
"""Independent optimizer-free Bernstein replay using independent symbolic columns.

Integer scaling retains the signs exactly. This checker imports neither the
primary symbolic engine nor the primary solve/recovery/check code.
"""
from __future__ import annotations
if not __debug__: raise RuntimeError('Independent interval certificate replay refuses optimized Python')
from pathlib import Path
from fractions import Fraction as F
from math import comb,lcm
import gzip,json,hashlib,time,argparse
import symbolic_features_independent as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
TNAME=['m2','m3','q2','q3','j']+[f'R_{t}_{a}_{b}_{c}'for t,a,b,c in s.COORDS]+[f'S_{d}_{x}_'+''.join(map(str,T))for d,x,T in s.SIGS]
ENAME=[f'{typ}_{d}_{h}_{k}'for d,h,k in s.CLASSES for typ in('X','B')]
LNAME=[f'L_{d}_{r}'for d,r in s.PATTERNS]
ALL=TNAME+ENAME+LNAME

def require(p,msg):
    if not p:raise RuntimeError(msg)
def tup(v):return tuple(tup(x)for x in v)if isinstance(v,list)else v
def plus(a,b):return s.trim([(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))])
def times(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):out[i+j]+=x*y
    return s.trim(out)
def scale(a,c):return s.trim([x*c for x in a])
def sub(a,b):return plus(a,scale(b,-1))
def ppow(a,n):
    ans=(1,)
    for _ in range(n):ans=times(ans,a)
    return ans

def cleared(K,names,assoc):
    d,Z,O,T,ext,pattern=K
    target=scale(sub(times((1,4,6),O),times((0,1,3),Z)),6)
    coords=dict(zip(TNAME,T))
    for (rd,h,k),(N,X,B),nx,nb in zip(s.CLASSES,ext,ENAME[::2],ENAME[1::2]):
        qnum=(0,1);qden=(1,1)
        # (1+z)^2(q N-X) and (1+z)^2(q^(h-k) X-B).
        polyX=sub(scale(times(qnum,qden),N),scale(ppow(qden,2),X))
        polyB=sub(scale(times(ppow(qnum,h-k),ppow(qden,2-h+k)),X),scale(ppow(qden,2),B))
        coords[nx]=scale(times(polyX,Z),6);coords[nb]=scale(times(polyB,Z),6)
    for j,(name,L)in enumerate(zip(LNAME,assoc)):
        coords[name]=scale(sub(times(L['D'],O),times(L['N'],Z)),6)if j==pattern else(0,)
    return target,[coords[n]for n in names]

def verify(path,outpath=None):
    start=time.time();C=json.loads(Path(path).read_text());a,b=F(C['a']),F(C['b']);degree=C['degree'];M=C['bernstein_degree'];names=C['support']
    expected_generation=json.loads((ROOT/'inputs/dependencies/enumeration_v3.json').read_text())['generation_sha256']
    require(C['generation_sha256']==expected_generation,'Certificate structural generation hash mismatch')
    frozen=json.loads((ROOT/'inputs/SUPPORT_COMPARISON.json').read_text())['union_coordinates']
    require(names==[x for x in frozen if x!='normalization'],'Support changed/reordered')
    require(0<=a<b<=1,'Claimed interval outside original domain')
    stages={(F(0),F(1),2,4),(F(0),F(1),4,8)}
    stages.update((left,right,6,12)for left,right in((F(0),F(1)),(F(0),F(1,2)),(F(1,2),F(1)),(F(1,4),F(1,2)),(F(1,2),F(3,4)),(F(3,4),F(1))))
    require((a,b,degree,C['elevation'])in stages,'Undeclared interval/degree/elevation combination')
    beta=[[F(x)for x in row]for row in C['bernstein_polynomials']]
    require(len(beta)==74 and all(len(v)==degree+1 for v in beta),'Coefficient shape mismatch')
    H=lcm(*(x.denominator for row in beta for x in row))
    require(H>0,'Coefficient common denominator nonpositive')
    P=[];multiplier_coefficients=0
    for name,row in zip(names,beta):
        if name in ENAME+LNAME:
            require(all(x>=0 for x in row),'Negative multiplier Bernstein coefficient: '+name)
            multiplier_coefficients+=len(row)
        bi=[int(x*H)for x in row]
        P.append(tuple(sum(bi[i]*comb(degree,i)*comb(degree-i,r-i)*(-1)**(r-i)for i in range(r+1))for r in range(degree+1)))
    assoc=s.independent_association()
    for L in assoc:require(L['D'][0]>0 and all(x>=0 for x in L['D']),'Association denominator sign failure')
    with gzip.open(HERE/'symbolic_columns_independent.json.gz','rt')as f:records=json.load(f)['rows']
    # Uniform z-degree scaling; this is exact, not a floating-point rescale.
    base=10
    B=lcm(a.denominator,b.denominator);A=int(a*B);D=int((b-a)*B)
    require(B>0 and D>0,'Affine change denominator/length not positive')
    change=[[comb(j,k)*A**(j-k)*D**k*B**(base-j)if k<=j else 0 for k in range(base+1)]for j in range(base+1)]
    def trans(poly):
        require(len(poly)<=base+1,'Unexpected cleared feature degree')
        return s.trim([sum(v*change[j][k]for j,v in enumerate(poly)if j>=k)for k in range(base+1)])
    require(M==base+degree+C['elevation'],'Incorrect frozen-stage Bernstein degree')
    divisors=[comb(M,j)for j in range(base+degree+1)];Q=lcm(*divisors)
    BM=[[comb(k,j)*(Q//divisors[j])if j<=k else 0 for j in range(base+degree+1)]for k in range(M+1)]
    zero=0;fullzero=0;checks=0;viewcount=0;sha=hashlib.sha256();c4matches=0
    for idx,record in enumerate(records):
        K=tup(record['key']);d,Z,O,T,ext,pattern=K
        require(Z[0]==1 and all(x>=0 for x in Z),'Local partition denominator sign failure')
        tar,fs=cleared(K,names,assoc)
        val=list(scale(trans(tar),H))+[0]*(base+degree+1-len(trans(tar)))
        for coeff,feature in zip(P,fs):
            if not any(feature)or not any(coeff):continue
            product=times(coeff,trans(feature))
            for j,x in enumerate(product):val[j]-=x
        val=s.trim(val)
        bs=[sum(x*BM[k][j]for j,x in enumerate(val))for k in range(M+1)]
        bad=next((k for k,x in enumerate(bs)if x<0),None)
        require(bad is None,'Negative full-view slack Bernstein coefficient at independent column '+str(idx)+' k='+str(bad))
        checks+=len(bs);viewcount+=record['full_views'];zero+=not any(val);fullzero+=(not any(val))*record['full_views']
        # Actual C4 pair: degree2, degree2 shared boundary; no other boundary class.
        c4ext=ext[s.CLASSES.index((2,2,2))]
        if d==2 and pattern==0 and all(not any(e)for j,e in enumerate(ext)if j!=s.CLASSES.index((2,2,2))) and c4ext in((1,0,0),(1,1,1)):
            if (c4ext==(1,0,0)and Z==(1,3,3)and O==(0,1,2))or(c4ext==(1,1,1)and Z==(1,3)and O==(0,1)):
                require(not any(val),'C4 slack not identically zero');c4matches+=1
        sha.update((json.dumps(list(val),separators=(',',':'))+'\n').encode())
    require(viewcount==231816 and len(records)==9867,'Independent coverage count mismatch')
    require(c4matches==2,'Did not identify exactly two C4 polynomial views')
    out=dict(all_passed=True,certificate=str(Path(path).resolve()),certificate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),interval=[str(a),str(b)],degree=degree,bernstein_degree=M,
             independent_full_views=viewcount,independent_columns=len(records),slack_bernstein_coefficients_checked=checks,multiplier_bernstein_coefficients_checked=multiplier_coefficients,
             zero_slack_columns=zero,zero_slack_full_views=fullzero,c4_identically_zero_views=c4matches,all_denominators_strictly_positive=True,
             integer_slack_polynomial_sequence_sha256=sha.hexdigest(),elapsed_seconds=time.time()-start,optimizer_used=False)
    dest=Path(outpath)if outpath else HERE/('certificate_'+Path(path).parent.name+'_audit.json')
    dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True);return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('certificate');p.add_argument('--out');args=p.parse_args();verify(args.certificate,args.out)
