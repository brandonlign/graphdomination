#!/usr/bin/env python3
"""Adversarial third implementation; read-only inputs, explicit fail-closed checks.
No project verifier, enumerator, or feature module is imported.
"""
from __future__ import annotations
import gzip, hashlib, itertools as it, json, math, sys, time
from collections import Counter, defaultdict
from fractions import Fraction as Q
from pathlib import Path

if not __debug__:
    raise RuntimeError('Challenge checks deliberately reject optimized Python')
ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

STARS = tuple((d, bit, ns) for d in (2, 3) for bit in (0, 1)
              for ns in it.combinations_with_replacement(range(4), d))
STAR_INDEX = {s: i for i, s in enumerate(STARS)}
WALKS = tuple((t, a, b, c) for t in (0, 1) for a in range(4)
              for b in range(4) for c in range(a + 1, 4))
WALK_INDEX = {w: i for i, w in enumerate(WALKS)}
CLASSES = tuple((d, h, k) for d in (2, 3) for h in (2, 3) for k in range(1, h + 1))
PATTERNS = ((2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2), (3, 3))
BOUNDS = tuple(map(Q, ('1/3', '9/25', '5/13', '79/207', '69/197', '15/47', '13/45')))
EQ_NAMES = ['normalization', 'm2', 'm3', 'q2', 'q3', 'j'] + [f'R_{t}_{a}_{b}_{c}' for t,a,b,c in WALKS] + [f'S_{d}_{b}_'+''.join(map(str,ns)) for d,b,ns in STARS]
G_NAMES = [f'{typ}_{d}_{h}_{k}' for d,h,k in CLASSES for typ in ('X','B')] + [f'L_{d}_{r}' for d,r in PATTERNS]

def symbol(deg, bit):
    return 2 * (deg - 2) + bit

def family_and_features(d, deg, edges, blocks):
    """blocks contain (tuple of inside neighbors, actual degree, x flag, b flag)."""
    adj = [set(range(1, d+1))] + [{0} for _ in range(d)]
    for a,b in edges:
        adj[a].add(b); adj[b].add(a)
    for u in range(1,d+1):
        require(len(adj[u]) + sum(u in T for T,h,x,b in blocks) == deg[u-1], 'stub count')
    family = []
    for bits in it.product((0,1), repeat=d+1):
        if all(bits):
            continue
        ok = True
        for u in range(1,d+1):
            if bits[u] and all(bits[w] for w in adj[u]) and all(x for T,h,x,b in blocks if u in T):
                ok = False; break
        if ok and any(b and all(bits[u] for u in T) for T,h,x,b in blocks):
            ok = False
        if ok:
            family.append(bits)
    z = len(family)
    require(z > 0, 'empty conditional family')
    occ = [sum(bits[u] for bits in family) for u in range(d+1)]
    joint = [sum(bits[0]*bits[u] for bits in family) for u in range(1,d+1)]
    # Every integer below has denominator 6*z until exterior/association rows.
    nums = [6*occ[0]]
    nums += [6*z*int(d==k) - (6//k)*z*sum(h==k for h in deg) for k in (2,3)]
    nums += [6*occ[0]*int(d==k) - (6//k)*sum(occ[u] for u in range(1,d+1) if deg[u-1]==k) for k in (2,3)]
    nums += [6*sum(joint[u-1]*(int(d==2 and deg[u-1]==3)-int(d==3 and deg[u-1]==2)) for u in range(1,d+1))]
    walknums = [0]*48; starnums=[0]*60
    for bits in family:
        syms = [symbol(d,bits[0])] + [symbol(h,bits[u]) for u,h in enumerate(deg,1)]
        rootstar = (d,bits[0],tuple(sorted(syms[u] for u in range(1,d+1))))
        starnums[STAR_INDEX[rootstar]] += 6
        for u in range(1,d+1):
            endpoints = [(syms[w],1) for w in adj[u] if w != 0]
            endpoints += [(symbol(h,x),0) for T,h,x,b in blocks if u in T]
            require(len(endpoints)==deg[u-1]-1, 'walk endpoints')
            for end,t in endpoints:
                a,c = syms[0],end
                if a != c:
                    ix = WALK_INDEX[t,min(a,c),syms[u],max(a,c)]
                    walknums[ix] += 6 if a<c else -6
            ns = [syms[0]] + [end for end,t in endpoints]
            star = (deg[u-1],bits[u],tuple(sorted(ns)))
            starnums[STAR_INDEX[star]] -= 6//deg[u-1]
    nums += walknums + starnums
    pattern = PATTERNS.index((d,sum(h==(3 if d==2 else 2) for h in deg)))
    D = math.lcm(6*z,4,BOUNDS[pattern].denominator)
    nums = [a*(D//(6*z)) for a in nums]
    for dr,h,k in CLASSES:
        chosen = [(x,b) for T,hh,x,b in blocks if d==dr and hh==h and len(T)==k]
        X=sum(x for x,b in chosen); B=sum(b for x,b in chosen)
        nums += [len(chosen)*(D//2)-X*D, X*(D//(2**(h-k)))-B*D]
    for i,L in enumerate(BOUNDS):
        nums += [occ[0]*(D//z)-L.numerator*(D//L.denominator) if i==pattern else 0]
    require(len(nums)==141, 'feature length')
    g=math.gcd(D,*nums)
    return tuple(family), (d,D//g,tuple(a//g for a in nums))

def block_partitions(stubs):
    if not stubs:
        yield (); return
    first,rest=stubs[0],stubs[1:]
    for extra_size in range(min(2,len(rest))+1):
        for extra in it.combinations(rest,extra_size):
            block=(first,)+extra
            if len({u for u,j in block}) != len(block):
                continue
            chosen=set(extra)
            for tail in block_partitions(tuple(s for s in rest if s not in chosen)):
                yield (block,)+tail

def all_structures():
    for d in (2,3):
        poss=tuple(it.combinations(range(1,d+1),2))
        for edge_bits in it.product((0,1),repeat=len(poss)):
            edges=tuple(e for e,b in zip(poss,edge_bits) if b)
            for deg in it.product((2,3),repeat=d):
                remaining=tuple(deg[u-1]-1-sum(u in e for e in edges) for u in range(1,d+1))
                if min(remaining)<0:
                    continue
                stubs=tuple((u,j) for u,a in enumerate(remaining,1) for j in range(a))
                for partition in block_partitions(stubs):
                    Ts=tuple(tuple(u for u,j in B) for B in partition)
                    opts=[tuple((h,x,b) for h in (2,3) if h>=len(T) for x,b in it.product((0,1),repeat=2)
                                if b<=x and (h!=len(T) or b==x)) for T in Ts]
                    for choice in it.product(*opts):
                        blocks=tuple(sorted((T,h,x,b) for T,(h,x,b) in zip(Ts,choice)))
                        yield d,deg,edges,blocks

def read_data():
    with gzip.open(ROOT/'views_v3.json.gz','rt') as f:
        data=json.load(f)
    require(data['equality_names']==EQ_NAMES and data['inequality_names']==G_NAMES,'coordinate metadata')
    rows=data['rows']; table={}
    for i,row in enumerate(rows):
        require(row['id']==i,'noncanonical row id')
        K=(row['d'],row['denominator'],tuple(row['numerators']))
        require(K not in table,'duplicate stored column')
        require(K[1]>0 and math.gcd(K[1],*K[2])==1,'noncanonical rational vector')
        table[K]=row['full_views']
        R=row['representative']
        blocks=tuple(sorted((tuple(T),h,x,b) for T,h,(x,b) in zip(R['blocks'],R['boundary_degrees'],R['flags'])))
        fam,K2=family_and_features(R['d'],tuple(R['neighbor_degrees']),tuple(map(tuple,R['neighbor_edges'])),blocks)
        require(K==K2,'stored feature/representative mismatch')
        require(sorted(sum(bit<<j for j,bit in enumerate(bits)) for bits in fam)==sorted(R['feasible_masks']),'representative family mismatch')
    return data,table

def verify_certificate(kind,data,table):
    cert=json.loads((ROOT/f'certificate_{kind}_v3.json').read_text())
    require(cert['problem']==kind and cert['generation_sha256']==data['generation_sha256'],'certificate identity')
    sec=(kind=='secondary'); y=tuple(map(Q,cert['equality_dual'])); mu=tuple(map(Q,cert['inequality_dual']))
    require(len(y)==114+int(sec) and len(mu)==27 and min(mu)>=0,'dual dimensions/signs')
    H=math.lcm(*(q.denominator for q in y+mu)); yi=[int(q*H) for q in y]; mi=[int(q*H) for q in mu]
    tight=0; min_slack=None
    for (d,D,N),multiplicity in table.items():
        target=-int(d==3)*D if sec else N[0]
        rhs=yi[0]*D+sum(yi[j]*N[j] for j in range(1,114))+sum(mi[j]*N[114+j] for j in range(27))
        if sec: rhs+=yi[-1]*N[0]
        slack=H*target-rhs
        require(slack>=0,f'{kind} dual negative slack')
        tight+=int(slack==0)
        s=Q(slack,H*D); min_slack=s if min_slack is None else min(s,min_slack)
    eq=[Q(0)]*len(y); gs=[Q(0)]*27; obj=Q(0); d3mass=Q(0)
    for s in cert['support']:
        w=Q(s['weight']); require(w>=0,'negative primal mass')
        row=data['rows'][s['id']]; f=tuple(Q(n,row['denominator']) for n in row['numerators'])
        require(list(map(str,f))==s['features'] and row['representative']==s['representative'],'support binding')
        eq[0]+=w
        for j in range(1,114): eq[j]+=w*f[j]
        if sec: eq[-1]+=w*f[0]
        for j in range(27): gs[j]+=w*f[114+j]
        obj+=w*(-int(row['d']==3) if sec else f[0]); d3mass+=w*int(row['d']==3)
    want=[Q(1)]+[Q(0)]*113+([Q(4,11)] if sec else [])
    require(eq==want and min(gs)>=0,'primal feasibility')
    dual=sum(a*b for a,b in zip(y,want)); target=Q(0) if sec else Q(4,11)
    require(obj==dual==Q(cert['objective'])==target,'primal dual objective mismatch')
    if sec:
        require(y[0]==-Q(4,11)*y[-1] and Q(0)<-y[-1]<72,'secondary quantitative implication')
    return dict(objective=str(obj),dual_objective=str(dual),tight_columns=tight,minimum_slack=str(min_slack),
                degree3_support_mass=str(d3mass),nonzero_inequality_multipliers=sum(q!=0 for q in mu),
                face_multiplier=str(y[-1]) if sec else None,normalization_multiplier=str(y[0]))

def fixture(name,n,edges,table):
    adj=[set() for _ in range(n)]
    for a,b in edges:
        require(a!=b and b not in adj[a],'fixture nonsimple')
        adj[a].add(b);adj[b].add(a)
    require(all(len(ns) in (2,3) for ns in adj),'fixture degree outside scope')
    closed=[{v}|adj[v] for v in range(n)]
    feasible=[frozenset(v for v in range(n) if bits[v]) for bits in it.product((0,1),repeat=n)
              if not any(all(bits[w] for w in C) for C in closed)]
    totals=[Q(0)]*141; totalweight=Q(0); exposures=0
    adjacent_boundary=shared_boundary=far_centers=0; pointwise_association_failures=0; pointwise_exterior_failures=0
    witnesses={}
    for v in range(n):
        local=(v,)+tuple(sorted(adj[v])); labels={w:i for i,w in enumerate(local)}; R=set(local)
        boundary=set().union(*(adj[u] for u in adj[v]))-R
        adjacent_boundary+=int(any(adj[w]&boundary for w in boundary))
        shared_boundary+=int(any(len(adj[w]&R)>1 for w in boundary))
        far_centers+=int(bool(set(range(n))-R-boundary))
        groups=defaultdict(list)
        for I in feasible: groups[I-R].append(I)
        for E,completions in groups.items():
            d=len(adj[v]);deg=tuple(len(adj[u]) for u in local[1:])
            inside_edges=tuple((labels[a],labels[b]) for a,b in it.combinations(local[1:],2) if b in adj[a])
            blocks=tuple(sorted((tuple(sorted(labels[u] for u in adj[w]&R)),len(adj[w]),int(w in E),int(({w}|(adj[w]-R))<=E)) for w in boundary))
            family,K=family_and_features(d,deg,inside_edges,blocks)
            actual={tuple(int(u in I) for u in local) for I in completions}
            require(set(family)==actual and len(family)==len(completions),'actual conditional family/law mismatch')
            require(K in table,'actual exposure missing from full universe')
            weight=Q(len(completions),n*len(feasible)); totalweight+=weight; exposures+=1
            features=tuple(Q(q,K[1]) for q in K[2])
            for j,f in enumerate(features):totals[j]+=weight*f
            pointwise_exterior_failures+=int(min(features[114:134])<0)
            pointwise_association_failures+=int(min(features[134:])<0)
            for label,indices in (('exterior',range(114,134)),('association',range(134,141))):
                bad=[dict(name=G_NAMES[j-114],value=str(features[j])) for j in indices if features[j]<0]
                if bad and label not in witnesses:
                    witnesses[label]=dict(root=v,exterior=sorted(E),root_occupancy=str(features[0]),negative_features=bad)
    p=Q(sum(len(I) for I in feasible),n*len(feasible))
    require(totalweight==1 and totals[0]==p,'root/exposure weights')
    require(all(t==0 for t in totals[1:114]),'actual equality feature nonzero')
    require(min(totals[114:])>=0,'actual averaged inequality negative')
    require(p>=Q(4,11),'fixture theorem bound')
    return dict(name=name,vertices=n,edges=sorted(map(list,edges)),admissible_sets=len(feasible),exposures=exposures,
                occupied_fraction=str(p),dominating_mean=str(n*(1-p)),degree3_fraction=str(Q(sum(len(ns)==3 for ns in adj),n)),
                roots_with_adjacent_boundary=adjacent_boundary,roots_with_shared_boundary=shared_boundary,
                roots_with_far_centers=far_centers,views_violating_pointwise_exterior_bound=pointwise_exterior_failures,
                views_violating_pointwise_association_bound=pointwise_association_failures,
                inequality_means=list(map(str,totals[114:])),pointwise_counterexamples=witnesses)

def scalar_checks():
    constants=[]
    for (d,r),want in zip(PATTERNS,BOUNDS):
        ds=([3]*r+[2]*(d-r)) if d==2 else ([2]*r+[3]*(d-r))
        aa=[1-Q(1,2**(h-1)) for h in ds]
        T=Q(sum(math.prod(aa[i] for i in range(d) if (mask>>i)&1) for mask in range((1<<d)-1)),1<<d)
        require(T/(1+T)==want,'association constant')
        constants.append(dict(pattern=[d,r],ratio=str(T),occupancy=str(want)))
    # Cycle induction needs only E3..E6 and positivity of the exact T identity.
    Z={1:1,2:3,3:7};W={1:1,2:4,3:12};cycle=[]
    for n in range(4,11):
        Z[n]=sum(Z[n-j] for j in (1,2,3));W[n]=Z[n]+sum(W[n-j] for j in (1,2,3))
    for n in range(3,11):
        E=7*n*Z[n]-11*W[n];require(E>=0 and (E==0)==(n==4),'cycle base/recurrence')
        if n>=7: require(-4*Z[n-1]+3*Z[n-2]+10*Z[n-3]==4*Z[n-5]+5*Z[n-6]>0,'cycle T identity')
        cycle.append(dict(n=n,Z=Z[n],W=W[n],E=E))
    return dict(association_constants=constants,cycle_base_checks=cycle)

def main():
    started=time.time();data,table=read_data();print('All stored representatives independently reconstructed.',flush=True)
    seen=Counter();cache={};counts=Counter()
    for structure in all_structures():
        if structure not in cache: cache[structure]=family_and_features(*structure)[1]
        K=cache[structure];require(K in table,'regenerated view not present')
        seen[K]+=1;counts[str(K[0])]+=1
    require(seen==table,'complete feature multiset/multiplicity mismatch')
    require(sum(seen.values())==231816 and len(seen)==9853,'unexpected universe size')
    print('Full universe and exact multiplicities matched.',flush=True)
    certificates={kind:verify_certificate(kind,data,table) for kind in ('primary','secondary')}
    c5=[(i,(i+1)%5) for i in range(5)]
    prism_subdivided=[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(1,4),(2,5),(0,6),(6,3)]
    c4_c5=[(i,(i+1)%4) for i in range(4)]+[(4+i,4+(i+1)%5) for i in range(5)]
    fixtures=[fixture('C5: adjacent distance-two vertices',5,c5,table),
              fixture('Triangular prism with one matching edge subdivided',7,prism_subdivided,table),
              fixture('C4 disjoint union C5: remote components',9,c4_c5,table)]
    inputs=['views_v3.json.gz','certificate_primary_v3.json','certificate_secondary_v3.json','local_views_v3.py','v1_universe_snapshot.py','v2_features_snapshot.py']
    result=dict(all_checks_passed=True,full_views=sum(seen.values()),columns=len(seen),full_views_by_degree=dict(counts),
                distinct_combinatorial_structures=len(cache),certificates=certificates,fixtures=fixtures,scalar_checks=scalar_checks(),
                input_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},elapsed_seconds=time.time()-started)
    (OUT/'check_exact_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('input_hashes','fixtures','scalar_checks')},indent=2))
if __name__=='__main__': main()
