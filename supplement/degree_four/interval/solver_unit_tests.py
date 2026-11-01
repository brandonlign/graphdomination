#!/usr/bin/env python3
"""Small exact tests that use no numerical optimizer."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Unit tests refuse Python -O/-OO')
from fractions import Fraction as F
from math import comb
import time,resource,signal,traceback
import symbolic_engine as s
import solver_math as m

def main():
    start=time.monotonic();tests=[]
    for degree in (0,1,2,4):
        for i in range(degree+1):
            m.require(m.power_to_bernstein(m.basis_power(i,degree),degree)==tuple(F(j==i) for j in range(degree+1)),'Basis identity')
        vals=[F(2*j-3,7) for j in range(degree+1)]
        m.require(m.power_to_bernstein(m.bernstein_to_power(vals),degree)==tuple(vals),'Basis roundtrip')
    tests.append({'test':'exact_Bernstein_power_roundtrips','passed':True})
    system=m.AffineSystem(4);system.add({0:2,1:3},7);system.add({1:2,3:5},11)
    free,base,T=system.affine_map();x=system.solve([F(1,3),F(2,7),F(3,11),F(4,13)])
    m.require(2*x[0]+3*x[1]==7 and 2*x[1]+5*x[3]==11,'Affine equation solution')
    m.require(x==[base[j]+sum((a*x[free[k]] for k,a in T[j].items()),F(0)) for j in range(4)],'Affine map')
    try:system.add({0:4,1:6},15)
    except m.AffineInconsistent:pass
    else:raise RuntimeError('Exact contradiction not rejected')
    fresh=m.AffineSystem(1)
    try:fresh.add({2:1},0)
    except m.CheckFailure:pass
    except m.AffineInconsistent:raise RuntimeError('Implementation-range error misclassified as inconsistency')
    else:raise RuntimeError('Bad affine coordinate not rejected')
    tests.append({'test':'consistent_affine_map_and_distinct_failure_classification','passed':True})
    data=s.load();m.require(data['full_views']==17553 and len(data['rows'])==1722,'Expected complete API dimensions')
    m.require(len(s.SUPPORT)==224 and len(s.EQNAMES)==170 and len(s.GNAMES)==54,'Full feature support')
    stages=[]
    for name in m.STAGES:
        p=m.Problem(name,data);cal,source=p.calibration();free,base,T=cal.affine_map()
        for eq in source:
            row={int(j):F(v) for j,v in eq['entries'].items()};rhs=F(eq['rhs'])
            m.require(sum((a*base[j] for j,a in row.items()),F(0))==rhs,'C4 affine particular')
            residual={}
            for j,a in row.items():
                for k,b in T[j].items():residual[k]=residual.get(k,F(0))+a*b
            m.require(all(v==0 for v in residual.values()),'C4 homogeneous map')
        # Distinct manual direct-power product, followed by the monomial elevation identity.
        coords=[F((j%9)-4,(j%5)+1) for j in range(p.nv)]
        selected=sorted(set([0,len(p.rows)//2,len(p.rows)-1]+p.c4));checks=0
        for c in selected:
            target,features=p.cleared[c];total=[F(0)]*(p.M+1)
            for j,Q in enumerate(features):
                for i in range(p.d+1):
                    beta=coords[j*(p.d+1)+i]
                    for h in range(p.d-i+1):
                        a=beta*comb(p.d,i)*comb(p.d-i,h)*(-1)**h
                        if a:
                            for k,q in enumerate(Q):
                                if q:total[i+h+k]+=a*q
            for r in range(p.M+1):
                expected=sum((a*F(comb(r,k),comb(p.M,k)) for k,a in enumerate(total) if a and k<=r),F(0))
                row,rhs=p.exact_row(c,r)
                observed=sum((a*coords[j] for j,a in row.items()),F(0))
                m.require(observed==expected,'Exact feature-product Bernstein row')
                manual_rhs=sum((F(a)*F(comb(r,k),comb(p.M,k)) for k,a in enumerate(target) if a and k<=r),F(0))
                m.require(rhs==manual_rhs,'Exact target Bernstein coefficient');checks+=1
        stages.append({'stage':name,'original_variables':p.nv,'columns':len(p.rows),'c4_columns':p.c4,
            'calibration_power_equations':len(source),'calibration_rank':len(cal.basis),'free_coordinates':len(free),
            'exact_product_rows_checked':checks})
    tests.append({'test':'actual_API_C4_affine_calibration_and_exact_product_rows','passed':True,'stages':stages})
    result={'all_passed':True,'optimizer_free':True,'tests':tests,'elapsed_seconds':time.monotonic()-start,
            'maximum_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    m.dump(s.HERE/'solver_unit_tests.json',result);print(result,flush=True)
if __name__=='__main__':
    try:resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    except (ValueError,OSError):pass  # macOS rejects RLIMIT_AS; this limit only guards memory use
    signal.signal(signal.SIGALRM,signal.SIG_DFL);signal.alarm(1800)
    try:main()
    except BaseException as err:
        m.dump(s.HERE/'solver_unit_tests_exception.json',{'type':type(err).__name__,'message':str(err),'traceback':traceback.format_exc()});raise
    finally:signal.alarm(0)
