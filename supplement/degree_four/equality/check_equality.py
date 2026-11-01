#!/usr/bin/env python3
"""Fixed-data equality audit. Standard library only; no source evaluator imports."""
from __future__ import annotations
if not __debug__:
    raise RuntimeError('Equality audit refuses Python -O/-OO')
from pathlib import Path
from fractions import Fraction as Q
from math import comb, lcm
from itertools import combinations
from collections import Counter
import copy, gzip, hashlib, json, time
BASE=Path(__file__).resolve().parent
INP=BASE/'inputs'

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def read(name):
    p=INP/name
    if name.endswith('.gz'):
        with gzip.open(p,'rt') as f:return json.load(f)
    return json.loads(p.read_text())
def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a or [0]
def plus(a,b):
    c=[0]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return trim(c)
def times(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)
def scale(a,k):return trim([k*x for x in a])
def minus(a,b):return plus(a,scale(b,-1))
def power(k):return [comb(k,j) for j in range(k+1)]
def shift(a,k):return trim([0]*k+list(a))
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def sh(x):return hashlib.sha256(canonical(x).encode()).hexdigest()
def dump(name,x): (BASE/name).write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')

def associations(patterns):
    # Product-minus-full-subset formula, expanded before exact polynomial division.
    ans=[]
    for d,hs in patterns:
        top=[1];excluded=[1]
        for h in hs:
            top=times(top,minus(power(h),shift([1],h)))
            excluded=times(excluded,minus(power(h-1),shift([1],h-1)))
        numerator=minus(top,shift(excluded,d));M=d+sum(h-1 for h in hs)-min(h-1 for h in hs)
        divisor=power(sum(hs)-M);rem=list(numerator);quot=[0]*max(1,len(rem)-len(divisor)+1)
        while len(rem)>=len(divisor) and any(rem):
            k=len(rem)-len(divisor);need(rem[-1]%divisor[-1]==0,'Association exact division')
            q=rem[-1]//divisor[-1];quot[k]=q;rem=minus(rem,shift(scale(divisor,q),k))
        need(not any(rem),'Association remainder');R=trim(quot);N=shift(R,1)
        ans.append(dict(M=M,R=R,N=N,D=plus(power(M),N)))
    return ans

def factors(row,classes,assoc):
    Z=row['Z'];O=row['O']
    target=scale(minus(times([1,4,6],O),times([0,1,3],Z)),12)
    out=[list(x) for x in row['F']]
    for (d,h,k),(N,X,B) in zip(classes,row['exterior']):
        out.append(scale(times(times(Z,power(2)),[-X,N-X]),12))
        a=h-k
        out.append(scale(times(times(Z,power(3-a)),minus(shift([X],a),scale(power(a),B))),12))
    for j,a in enumerate(assoc):
        out.append(scale(minus(times(a['D'],O),times(a['N'],Z)),12) if j==row['pattern'] else [0])
    need(len(out)==224,'Cleared factors dimension')
    return target,out

def certificate_powers(cert,support,digest):
    need(cert['schema_version']=='degree_four_weighted_interval_v1','Schema')
    need(cert['degree']==2 and cert['bernstein_degree']==16 and cert['stage']=='global_d2','Fixed degrees')
    need(cert['interval']==['0','1'] and cert['clearing_factor']=='12*A*Z_C','Interval/clearing')
    need(cert['generation_sha256']==digest and cert['support']==support,'Support/digest')
    vals=[[Q(s) for s in p] for p in cert['bernstein_polynomials']]
    need(len(vals)==224 and all(len(p)==3 for p in vals),'Multiplier dimensions')
    need(all(v>=0 for p in vals[170:] for v in p),'Negative inequality multiplier')
    den=lcm(*(v.denominator for p in vals for v in p))
    pows=[]
    for p in vals:
        b=[int(v*den) for v in p]
        pows.append([b[0],2*(b[1]-b[0]),b[0]-2*b[1]+b[2]])
    return den,pows

def slacks(rows,classes,assoc,den,pows):
    zeros=[];ends=[];witnesses=[];digest=hashlib.sha256();slack_powers=[]
    for c,row in enumerate(rows):
        target,ff=factors(row,classes,assoc);p=scale(target,den)
        for multiplier,f in zip(pows,ff):p=minus(p,times(multiplier,f))
        need(len(p)<=17,'Slack degree exceeds sixteen')
        bern=[sum((Q(p[k]*comb(i,k),den*comb(16,k)) for k in range(min(i,len(p)-1)+1)),Q()) for i in range(17)]
        need(all(b>=0 for b in bern),f'Negative slack coefficient in column {c}')
        if not any(p):zeros.append(c)
        else:witnesses.append(dict(column=c,positive_bernstein_index=next(i for i,b in enumerate(bern) if b>0)))
        if bern[-1]==0:ends.append(c)
        digest.update(('|'.join(map(str,bern))+'\n').encode());slack_powers.append(p)
    return zeros,ends,witnesses,digest.hexdigest(),slack_powers

def inspect_mapping(data,views):
    need(len(data['full_to_column'])==len(views)==17553,'Fixed full-view coverage')
    by=[[] for r in data['rows']]
    for i,(v,c) in enumerate(zip(views,data['full_to_column'])):
        need(0<=c<len(by),'Column mapping range');by[c].append(i)
        need(v['d']==data['rows'][c]['d'],f'Root-degree mapping {i}')
    for i,(r,indices) in enumerate(zip(data['rows'],by)):
        need(r['id']==i and r['source_indices']==indices and r['full_views']==len(indices),'Inverse mapping/multiplicity')
    return by

def reconstruct(v,support,classes,patterns):
    # Rebuild an explicit partial graph, all clauses, and configuration transports.
    d=v['d'];need(d==len(v['degrees']) and d in (2,3,4),'Root degree inconsistency')
    n=d+1;adj={j:set() for j in range(n+len(v['blocks']))};deg={0:d}
    for u,h in enumerate(v['degrees'],1):
        adj[0].add(u);adj[u].add(0);deg[u]=h
    for u,w in v['edges']:adj[u].add(w);adj[w].add(u)
    boundary=[]
    for j,(T,(h,x,b)) in enumerate(zip(v['blocks'],v['choices']),n):
        need(len(set(T))==len(T) and h>=len(T) and h in (2,3,4),'Boundary degree/incidences')
        need((x,b) in ((0,0),(1,0),(1,1)) and (h!=len(T) or x==b),'Boundary flags')
        deg[j]=h;boundary.append((j,set(T),h,x,b))
        for u in T:adj[u].add(j);adj[j].add(u)
    need(all(len(adj[u])==deg[u] for u in range(n)),'Actual inside degrees')
    need(all(deg[u]==2 or deg[w]==2 for u in adj for w in adj[u]),'High-high edge')
    stubs=[u for u in range(1,n) for _ in range(deg[u]-1-sum(u in e for e in v['edges']))]
    need(stubs==v['stubs'],'Stub reconstruction')
    need(sorted(j for P in v['parts'] for j in P)==list(range(len(stubs))),'Stub partition coverage')
    need([[stubs[j] for j in P] for P in v['parts']]==v['blocks'],'Boundary block reconstruction')
    clauses=[set(range(n))]
    for u in range(1,n):
        if all(x for w,T,h,x,b in boundary if u in T):clauses.append(({u}|adj[u])&set(range(n)))
    clauses.extend(T for w,T,h,x,b in boundary if b)
    masks=[]
    for mask in range(1<<n):
        S={i for i in range(n) if mask>>i&1}
        if not any(C<=S for C in clauses):masks.append(mask)
    need(masks==v['masks'],'Recorded mask family incorrect')
    Z=[0]*n;O=[0]*n;F={name:[0]*n for name in support[:170]}
    for mask in masks:
        S={i for i in range(n) if mask>>i&1};k=len(S);Z[k]+=1;O[k]+=int(0 in S)
        occ=S|{w for w,T,h,x,b in boundary if x};sym={u:2*(h-2)+int(u in occ) for u,h in deg.items()}
        val={name:0 for name in F}
        for h in (2,3,4):
            val[f'm{h}']=12*int(d==h)-sum(12//h for u in adj[0] if deg[u]==h)
            val[f'q{h}']=12*int(d==h and 0 in S)-sum(12//h for u in adj[0] if deg[u]==h and u in S)
        for h in (3,4):val[f'J2{h}']=12*sum(int(0 in S and u in S)*(int(d==2 and deg[u]==h)-int(d==h and deg[u]==2)) for u in adj[0])
        for u in adj[0]:
            for w in adj[u]-{0}:
                a,b,c=sym[0],sym[u],sym[w]
                if a!=c:val[f'R_{int(w in adj[0])}_{min(a,c)}_{b}_{max(a,c)}']+=12 if a<c else -12
        for u,weight in [(0,12)]+[(u,-12//deg[u]) for u in adj[0]]:
            key=f'S_{deg[u]}_{int(u in S)}_'+','.join(map(str,sorted(sym[w] for w in adj[u])))
            val[key]+=weight
        for w,T,h,x,b in boundary:
            kT=len(T)
            val[f'BR_X_{d}_{h}_{kT}']+=12*x;val[f'BR_B_{d}_{h}_{kT}']+=12*b
            val[f'BR_X_{h}_{d}_{kT}']-=12*int(0 in S)
            val[f'BR_B_{h}_{d}_{kT}']-=12*int(({0}|(adj[0]-T))<=S)
        for name,x in val.items():F[name][k]+=x
    ext=[]
    for dd,h,k in classes:
        bb=[(x,b) for w,T,hh,x,b in boundary if dd==d and hh==h and len(T)==k]
        ext.append([len(bb),sum(x for x,b in bb),sum(b for x,b in bb)])
    pattern=patterns.index([d,sorted(v['degrees'])])
    return dict(d=d,Z=trim(Z),O=trim(O),F=[trim(F[name]) for name in support[:170]],exterior=ext,pattern=pattern)

def same_feature(a,b):return all(a[k]==b[k] for k in ('d','Z','O','F','exterior','pattern'))
def verify_zero_claim(claim,actual):need(claim==actual,'Zero-support claim does not match arithmetic')
def expect_reject(name,fn):
    try:fn()
    except (RuntimeError,KeyError,IndexError) as e:return dict(name=name,rejected=True,reason=str(e))
    raise RuntimeError('Corruption not rejected: '+name)

def main():
    begin=time.monotonic();data=read('primary_columns.json.gz');ind=read('independent_columns.json.gz');cert=read('certificate.json')
    raw=gzip.open(INP/'views.jsonl.gz','rb').read();views=[json.loads(line) for line in raw.splitlines()]
    need(hashlib.sha256(raw).hexdigest()==data['generation_sha256'],'Raw universe digest')
    need(len(data['rows'])==len(ind['rows'])==1722,'Fixed column count')
    inverse=inspect_mapping(data,views);inspect_mapping(ind,views)
    need(data['full_to_column']==ind['full_to_column'],'Independent mapping mismatch')
    need(data['classes']==ind['classes'] and data['patterns']==ind['patterns'],'Independent schema mismatch')
    for a,b in zip(data['rows'],ind['rows']):
        need(a['d']==b['d'] and a['Z']==b['Z'] and a['exterior']==b['exterior'] and a['pattern']==b['pattern'],'Independent row structure')
        need(scale(a['O'],12)==b['features12'][0] and a['F']==b['features12'][1:],'Independent feature mismatch')
    assoc=associations(data['patterns']);need(assoc==data['association']==ind['association'],'Association formulas differ')
    den,pows=certificate_powers(cert,data['support'],data['generation_sha256'])
    zeros,end,witnesses,digest,sp=slacks(data['rows'],data['classes'],assoc,den,pows)
    primary=read('primary_replay.json');audit=read('independent_certificate_verification.json')
    need(digest==primary['exact']['slack_certificate_sha256']==audit['slack_bernstein_sha256'],'Existing slack replay mismatch')
    need(zeros==primary['exact']['tight_columns'] and len(zeros)==audit['tight_columns'],'Existing zero-support mismatch')
    zero_raw=[i for i,c in enumerate(data['full_to_column']) if c in zeros];records=[]
    for i in zero_raw:
        c=data['full_to_column'][i];v=views[i];r=reconstruct(v,data['support'],data['classes'],data['patterns'])
        need(same_feature(r,data['rows'][c]),f'Zero-view independent features {i}')
        target,ff=factors(r,data['classes'],assoc);actual=scale(target,den)
        for p,f in zip(pows,ff):actual=minus(actual,times(p,f))
        need(not any(actual),'Reconstructed zero-view slack nonzero')
        records.append(dict(raw_index=i,column=c,view=v,reconstructed_features=r,features_sha256=sh(r)))
    need(len(zero_raw)==audit['tight_full_views'],'Existing raw zero count mismatch')
    need(all(r['view']['d']==2 for r in records),'Zero view has high root')
    # Same finite arithmetic supplies separate endpoint data, never inferred from interior.
    endpoint_degrees=dict(Counter(str(data['rows'][c]['d']) for c in end))
    endpoint_raw=[i for i,c in enumerate(data['full_to_column']) if c in end]
    # Intentional in-memory corruptions; the sealed sources are not touched.
    bad=copy.deepcopy(data);bad['full_to_column'][zero_raw[0]]=next(i for i,r in enumerate(data['rows']) if r['d']==3)
    badmask=copy.deepcopy(views[zero_raw[0]]);badmask['masks'].remove(0)
    badd=copy.deepcopy(views[zero_raw[0]]);badd['d']=3
    badcert=copy.deepcopy(cert);badcert['bernstein_polynomials'][170][0]='-1'
    controls=[expect_reject('missing_zero_column',lambda:verify_zero_claim(zeros[1:],zeros)),
              expect_reject('spurious_zero_column',lambda:verify_zero_claim(sorted(zeros+[next(i for i in range(1722) if i not in zeros)]),zeros)),
              expect_reject('raw_mapping_root_degree',lambda:inspect_mapping(bad,views)),
              expect_reject('removed_empty_inside_mask',lambda:reconstruct(badmask,data['support'],data['classes'],data['patterns'])),
              expect_reject('changed_root_degree',lambda:reconstruct(badd,data['support'],data['classes'],data['patterns'])),
              expect_reject('negative_inequality_multiplier',lambda:certificate_powers(badcert,data['support'],data['generation_sha256']))]
    dump('ZERO_SUPPORT.json',dict(columns=zeros,raw_indices=zero_raw,full_views=records))
    dump('STRICTNESS_WITNESSES.json',dict(nonzero_slack_columns=len(witnesses),positive_bernstein_witnesses=witnesses))
    dump('ENDPOINT_SUPPORT.json',dict(activity='z=1',columns=end,raw_indices=endpoint_raw,column_root_degrees=endpoint_degrees,
                                   raw_root_degrees=dict(Counter(str(views[i]['d']) for i in endpoint_raw))))
    dump('CONTROLS.json',dict(all_rejected=True,controls=controls))
    report=dict(all_passed=True,optimizer_free=True,no_source_evaluator_imported=True,no_new_grammar_or_graph_census=True,
                columns=1722,full_views=17553,checked_slack_coefficients=1722*17,checked_inequality_multiplier_coefficients=54*3,
                slack_bernstein_sha256=digest,identically_zero_columns=zeros,identically_zero_full_views=zero_raw,
                zero_column_root_degrees=dict(Counter(str(data['rows'][c]['d']) for c in zeros)),
                zero_full_view_root_degrees=dict(Counter(str(views[i]['d']) for i in zero_raw)),
                reconstructed_zero_view_feature_polynomials=len(zero_raw)*172,nonzero_slack_columns=1722-len(zeros),
                endpoint_zero_columns=len(end),endpoint_zero_full_views=len(endpoint_raw),endpoint_column_root_degrees=endpoint_degrees,
                source_comparison='All primary and prior independent rows/maps agree; all zero views freshly reconstructed',
                elapsed_seconds=time.monotonic()-begin)
    dump('EXACT_CHECKS.json',report);print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
