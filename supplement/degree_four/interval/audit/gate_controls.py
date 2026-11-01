#!/usr/bin/env python3
"""Exact independent reference/unit/negative controls; no primary code import."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Independent controls refuse optimized Python')
from pathlib import Path
from fractions import Fraction as F
import copy,gzip,json,subprocess,sys,time
import independent_weighted as a
HERE=Path(__file__).resolve().parent

def main():
    start=time.monotonic();controls=[]
    def reject(name,fn):
        try:fn()
        except (RuntimeError,ValueError,KeyError,TypeError) as e:controls.append(dict(name=name,rejected=True,message=str(e)));return
        raise RuntimeError('Failed rejection control: '+name)
    source=[json.loads(x) for x in gzip.open(a.source_path(),'rt')]
    V=source[0]
    bad=copy.deepcopy(V);bad['d']=5;reject('invalid_root_degree',lambda:a.evaluate(bad))
    bad=copy.deepcopy(V);bad['masks']=bad['masks']+[2**(bad['d']+1)-1];reject('corrupt_local_family',lambda:a.evaluate(bad))
    B=next(copy.deepcopy(v) for v in source if v['blocks']);B['blocks'][0].append(B['blocks'][0][0]);reject('repeated_neighbor_in_boundary_block',lambda:a.evaluate(B))
    B=next(copy.deepcopy(v) for v in source if any(len(T)==h for T,(h,x,b) in zip(v['blocks'],v['choices'])))
    i=next(i for i,(T,(h,x,b)) in enumerate(zip(B['blocks'],B['choices'])) if len(T)==h);B['choices'][i][1:]=[1,0];reject('exhausted_boundary_flags',lambda:a.evaluate(B))
    # Exact basis reconstruction, rather than an activity grid.
    tests=0
    for M in (2,4,16,20):
        for k in range(M+1):
            p=(0,)*k+(1,);a.require(a.basis_to_power(a.bern(p,M))==p,'Monomial Bernstein reconstruction');tests+=1
        p=tuple(F((i%5)-2,i+1) for i in range(M+1));a.require(a.basis_to_power(a.bern(p,M))==a.trim(p),'Rational Bernstein reconstruction');tests+=1
    with gzip.open(HERE/'independent_columns.json.gz','rt') as f:ind=json.load(f)
    with gzip.open(a.ROOT/'inputs/columns.json.gz','rt') as f:point=json.load(f)
    computed=[]
    for row in ind['rows']:
        d,Z,F12,ext,p=a.key_from_row(row);zz=sum(Z);v=[F(sum(poly),12*zz) for poly in F12]
        for (rd,h,k),(N,X,B) in zip(a.CLASSES,ext):v.extend((F(N,2)-X,F(X,2**(h-k))-B))
        for j,ass in enumerate(a.ASSOCIATION):v.append(v[0]-F(sum(ass['N']),sum(ass['D'])) if p==j else F(0))
        computed.append(tuple(v))
    pointvals=[tuple(F(x,r['denominator']) for x in r['numerators']) for r in point['rows']]
    for raw,(i,j) in enumerate(zip(ind['full_to_column'],point['full_to_column'])):a.require(computed[i]==pointvals[j],f'Full point-reference recovery {raw}')
    expected=(F(1,3),F(9,25),F(19,51),F(5,13),F(21,53),F(11,27),F(13,45),F(5,21))
    a.require(tuple(F(sum(p['N']),sum(p['D'])) for p in a.ASSOCIATION)==expected,'Fresh degree-four association constants')
    for script in ('independent_weighted.py','gate_controls.py'):
        r=subprocess.run([sys.executable,'-O',str(HERE/script)],capture_output=True,text=True)
        (HERE/(script+'.optimized_rejection.log')).write_text(r.stdout+r.stderr)
        a.require(r.returncode!=0 and 'refuse' in r.stderr.lower(),'-O rejection '+script)
        controls.append(dict(name='optimized_'+script,rejected=True,returncode=r.returncode))
    out=dict(all_passed=True,controls=controls,bernstein_polynomial_identity_tests=tests,point_reference_full_feature_entries=17553*225,point_reference_only_not_interval_evidence=True,elapsed_seconds=time.monotonic()-start)
    a.dump(HERE/'gate_controls.json',out);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
