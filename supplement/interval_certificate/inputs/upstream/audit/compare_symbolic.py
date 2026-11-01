#!/usr/bin/env python3
"""Compare independent complete-polynomial columns with primary and old anchors."""
from __future__ import annotations
if not __debug__: raise RuntimeError('Independent symbolic comparison refuses optimized Python')
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from functools import reduce
from math import gcd,lcm
import gzip,json,hashlib,time
import symbolic_features_independent as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent

def require(p,msg):
    if not p:raise RuntimeError(msg)
def tuplify(v):
    if isinstance(v,(list,tuple)):return tuple(tuplify(x)for x in v)
    if isinstance(v,dict):return tuple(sorted((k,tuplify(x))for k,x in v.items()))
    return v
def primary_key(r):
    require(all(x%6==0 for x in r['features6'][0]),'Objective numerator not six-divisible')
    return r['d'],s.trim(r['Z']),s.trim([x//6 for x in r['features6'][0]]),tuple(s.trim(p)for p in r['features6'][1:]),tuplify(r['exterior']),r['pattern']
def load(path):
    with gzip.open(path,'rt')as f:return json.load(f)
def canonical(d,v):
    den=lcm(*(f.denominator for f in v));nums=tuple(int(f*den)for f in v)
    g=reduce(gcd,nums,den)
    return d,den//g,tuple(x//g for x in nums)
def run():
    start=time.time();ind=load(HERE/'symbolic_columns_independent.json.gz')['rows'];primary=load(ROOT/'symbolic_columns.json.gz')
    ic=Counter({tuplify(r['key']):r['full_views']for r in ind});pc=Counter({primary_key(r):r['full_views']for r in primary['rows']})
    require(len(pc)==len(primary['rows']),'Primary rows merge under independent key normalization')
    require(ic==pc,'Complete symbolic tuples or full-view multiplicities disagree')
    require(sum(ic.values())==231816,'Full-view count mismatch')
    own=s.independent_association()
    require(tuplify(own)==tuplify(primary['association']),'Association polynomials disagree')
    expected=json.loads((ROOT/'inputs/dependencies/enumeration_v3.json').read_text())['generation_sha256']
    require(primary['generation_sha256']==expected,'Primary structural hash changed')
    anchors=[]
    for z,tag in((F(1,4),'z_1_4'),(F(1,2),'z_1_2'),(F(1),'z_1')):
        at=Counter()
        for K,n in ic.items():at[canonical(K[0],s.eval_complete(K,z,own))]+=n
        old=load(ROOT/'inputs/anchors'/tag/'views.json.gz')
        want=Counter({(r['d'],r['denominator'],tuple(r['numerators'])):r['full_views']for r in old['rows']})
        require(at==want,'Anchor complete feature multiplicity discrepancy at '+str(z))
        anchors.append(dict(activity=str(z),columns=len(at),full_views=sum(at.values()),all_passed=True))
        print('Anchor passed',str(z),len(at),flush=True)
    out=dict(all_passed=True,full_views=sum(ic.values()),symbolic_columns=len(ic),complete_symbolic_multiplicities_match=True,all_seven_association_polynomials_match=True,
             primary_generation_sha256=expected,independent_enumerator='ordinary-set recursion; structural multiplicities compared without borrowing primary ordering',anchors=anchors,elapsed_seconds=time.time()-start)
    (HERE/'symbolic_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':run()
