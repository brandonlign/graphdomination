#!/usr/bin/env python3
"""Independent ordinary-set evaluator and optimizer-free audit; no primary imports."""
from __future__ import annotations
if not __debug__:
    raise RuntimeError('Independent exact audit refuses Python -O/-OO')
from fractions import Fraction as F
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement, product
from functools import reduce
from math import gcd, lcm
from pathlib import Path
import argparse, gzip, hashlib, json, time, copy

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'inputs/views.jsonl.gz'
DIGEST='af0020c5efce7b021f901e798013b34aa7967c019c4e5c98b8352e0e2985f021'
DEGREES=(2,3,4)
def require(ok,message):
    if not ok: raise RuntimeError(message)
def dump(path,obj): Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def clean(v): return json.dumps(v,sort_keys=True,separators=(',',':'))
def sigdegree(s): return 2+s//2
def edge_allowed(a,b): return a==2 or b==2
REV=tuple((t,a,b,c) for t in (0,1) for a in range(6) for b in range(6) for c in range(a+1,6)
          if edge_allowed(sigdegree(a),sigdegree(b)) and edge_allowed(sigdegree(b),sigdegree(c))
          and (not t or edge_allowed(sigdegree(a),sigdegree(c))))
STARS=tuple((d,x,T) for d in DEGREES for x in (0,1)
            for T in combinations_with_replacement(range(6) if d==2 else range(2),d))
CLASSES=tuple((d,h,k) for d in DEGREES for h in DEGREES for k in range(1,min(d,h)+1))
PATTERNS=tuple((2,T) for T in combinations_with_replacement(DEGREES,2))+((3,(2,2,2)),(4,(2,2,2,2)))
DESCRIPTORS=(('c',),)+tuple(('m',h) for h in DEGREES)+tuple(('q',h) for h in DEGREES)+tuple(('J',2,h) for h in (3,4))+tuple(('R',)+r for r in REV)+tuple(('S',)+s for s in STARS)+tuple(('BR',d,h,k,kind) for d,h,k in CLASSES for kind in ('X','B'))+tuple(('D',d,h,k,kind) for d,h,k in CLASSES for kind in ('X','B'))+tuple(('L',)+p for p in PATTERNS)
INDEX={r:i for i,r in enumerate(DESCRIPTORS)}
NEQ=171; NFEAT=225
require(len(REV)==56 and len(STARS)==60 and len(CLASSES)==23 and len(PATTERNS)==8 and len(DESCRIPTORS)==NFEAT,'Independent schema dimensions')

def subsets(elements):
    seq=tuple(elements)
    for size in range(len(seq)+1):
        for x in combinations(seq,size): yield frozenset(x)

def association(neighbor_degrees):
    # Direct proper-subset sum, independently of the product-minus-product formula.
    d=len(neighbor_degrees); T=F(0)
    for selected in subsets(range(d)):
        if len(selected)==d: continue
        prob=F(1,2**d)
        for j in selected: prob*=1-F(1,2**(neighbor_degrees[j]-1))
        T+=prob
    return T/(1+T)
BOUNDS={p:association(p[1]) for p in PATTERNS}

def validate_structure(V):
    d=V['d']; deg=V['degrees']; nb=set(range(1,d+1)); edges=V['edges']; parts=V['parts']; blocks=V['blocks']; choices=V['choices']; stubs=V['stubs']
    require(d in DEGREES and len(deg)==d and all(h in DEGREES for h in deg),'Root/neighbor degree')
    require(all(edge_allowed(d,h) for h in deg),'High-high root edge')
    ed={frozenset(e) for e in edges}
    require(len(ed)==len(edges) and all(len(e)==2 and e<=nb for e in ed),'Inside edges')
    require(all(edge_allowed(deg[a-1],deg[b-1]) for a,b in edges),'High-high inside edge')
    expected=[u for u in range(1,d+1) for _ in range(deg[u-1]-1-sum(u in e for e in ed))]
    require(list(stubs)==expected and all(deg[u-1]>=1+sum(u in e for e in ed) for u in nb),'Stub counts')
    require(len(parts)==len(blocks)==len(choices),'Boundary lengths')
    require(sorted(j for B in parts for j in B)==list(range(len(stubs))),'Partition missing/repeated stub')
    require(all(list(B)==sorted(B) and B for B in parts),'Part ordering')
    require([B[0] for B in parts]==sorted(B[0] for B in parts),'Block ordering')
    require([list(stubs[j] for j in P) for P in parts]==[list(T) for T in blocks],'Block incidence disagreement')
    for T,(h,x,b) in zip(blocks,choices):
        require(len(set(T))==len(T) and set(T)<=nb and T,'Repeated/invalid neighbor in block')
        require(h in DEGREES and h>=len(T),'Boundary degree below incidence')
        require(all(edge_allowed(h,deg[u-1]) for u in T),'High-high boundary edge')
        require((x,b) in ((0,0),(1,0),(1,1)) and (h!=len(T) or x==b),'Invalid boundary flag')

def inside_sets(V):
    validate_structure(V)
    d=V['d']; universe=set(range(d+1)); edges=list(map(set,V['edges'])); blocked=[universe]
    for u in range(1,d+1):
        exposed=[x for T,(h,x,b) in zip(V['blocks'],V['choices']) if u in T]
        if all(exposed):
            clause={0,u}
            for e in edges:
                if u in e: clause|=e
            blocked.append(clause)
    blocked.extend(set(T) for T,(h,x,b) in zip(V['blocks'],V['choices']) if b)
    return tuple(S for S in subsets(range(d+1)) if not any(Q<=S for Q in blocked))

def independent_views():
    def partitions(stubs):
        def go(rem):
            if not rem:
                yield ()
                return
            first=min(rem); tail=sorted(rem-{first})
            for count in range(min(3,len(tail))+1):
                for other in combinations(tail,count):
                    part=(first,)+other
                    if len({stubs[i] for i in part})<len(part): continue
                    for rest in go(rem-set(part)): yield (part,)+rest
        yield from go(set(range(len(stubs))))
    for d in DEGREES:
        for deg in product(DEGREES,repeat=d):
            if not all(edge_allowed(d,h) for h in deg): continue
            for eset in subsets(combinations(range(1,d+1),2)):
                edges=sorted(eset)
                if any(not edge_allowed(deg[a-1],deg[b-1]) for a,b in edges): continue
                rem=[deg[u-1]-1-sum(u in e for e in edges) for u in range(1,d+1)]
                if min(rem)<0: continue
                stubs=[u for u,n in enumerate(rem,1) for _ in range(n)]
                for parts in partitions(stubs):
                    blocks=[tuple(stubs[j] for j in P) for P in parts]
                    options=[]
                    for T in blocks:
                        options.append([(h,x,b) for h in DEGREES for x,b in ((0,0),(1,0),(1,1))
                                        if h>=len(T) and all(edge_allowed(h,deg[u-1]) for u in T) and (h!=len(T) or x==b)])
                    for choices in product(*options):
                        V=dict(d=d,degrees=list(deg),edges=[list(e) for e in edges],stubs=stubs,parts=[list(p) for p in parts],blocks=[list(t) for t in blocks],choices=[list(c) for c in choices])
                        V['masks']=sorted(sum(2**j for j in S) for S in inside_sets(V))
                        yield V

def feature_vector(V):
    """Direct conditional average, sets plus rational coefficients; no bitwise features."""
    allowed=inside_sets(V)
    expected=sorted(sum(2**j for j in S) for S in allowed)
    require(expected==V['masks'],'Recorded local masks disagree with independent clauses')
    d=V['d']; deg={0:d,**{u:h for u,h in enumerate(V['degrees'],1)}}
    local={u:set() for u in range(d+1)}
    for u in range(1,d+1): local[0].add(u); local[u].add(0)
    for a,b in V['edges']: local[a].add(b); local[b].add(a)
    boundary=[]
    for j,(T,(h,x,b)) in enumerate(zip(V['blocks'],V['choices']),d+1):
        deg[j]=h; boundary.append((j,set(T),h,x,b))
        for u in T: local[u].add(j)
    accum=defaultdict(F)
    for S in allowed:
        occupied=set(S)|{j for j,T,h,x,b in boundary if x}
        symbols={u:2*(h-2)+int(u in occupied) for u,h in deg.items()}
        accum[('c',)]+=int(0 in S)
        for h in DEGREES:
            accum[('m',h)]+=int(d==h)-sum((F(1,deg[u]) for u in local[0] if deg[u]==h),F(0))
            accum[('q',h)]+=int(d==h and 0 in S)-sum((F(int(u in S),deg[u]) for u in local[0] if deg[u]==h),F(0))
        for h in (3,4):
            accum[('J',2,h)]+=sum(int(0 in S and u in S)*((d==2 and deg[u]==h)-(d==h and deg[u]==2)) for u in local[0])
        for u in local[0]:
            for w in local[u]-{0}:
                a,b,c=symbols[0],symbols[u],symbols[w]
                if a==c: continue
                desc=('R',int(w in local[0]),min(a,c),b,max(a,c))
                require(desc in INDEX,'Unknown reversal encountered')
                accum[desc]+=1 if a<c else -1
        signature=(d,int(0 in S),tuple(sorted(symbols[u] for u in local[0])))
        accum[('S',)+signature]+=1
        for u in local[0]:
            star=(deg[u],int(u in S),tuple(sorted(symbols[w] for w in local[u])))
            accum[('S',)+star]-=F(1,deg[u])
        for j,T,h,x,b in boundary:
            k=len(T)
            accum[('BR',d,h,k,'X')]+=x
            accum[('BR',d,h,k,'B')]+=b
            accum[('BR',h,d,k,'X')]-=int(0 in S)
            accum[('BR',h,d,k,'B')]-=int(({0}|(local[0]-T))<=S)
            accum[('D',d,h,k,'X')]+=F(1,2)-x
            accum[('D',d,h,k,'B')]+=F(x,2**(h-k))-b
        pattern=(d,tuple(sorted(V['degrees'])))
        accum[('L',)+pattern]+=int(0 in S)-BOUNDS[pattern]
    require(all(k in INDEX for k in accum),'Unknown independent feature')
    return tuple(accum[k]/len(allowed) for k in DESCRIPTORS)

def feature_key(V):
    vec=feature_vector(V); den=lcm(*(x.denominator for x in vec)); nums=tuple(int(x*den) for x in vec)
    g=reduce(gcd,nums,den)
    return V['d'],den//g,tuple(x//g for x in nums)

def fixtures():
    out={f'C{n}':(n,[(u,(u+1)%n) for u in range(n)]) for n in (3,4,5)}
    out.update({f'K2,{k}':(k+2,[(u,w) for u in (0,1) for w in range(2,k+2)]) for k in (3,4)})
    out['two_triangles_one_root']=(5,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)])
    out['triangle_and_square_one_root']=(6,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,5),(5,0)])
    out['once_subdivision_K4']=(10,[(a,w) for w,e in enumerate(combinations(range(4),2),4) for a in e])
    out['three_once_one_twice']=(7,[(u,w) for u in (0,1) for w in (2,3,4)]+[(0,5),(5,6),(6,1)])
    out['two_once_two_returning_triangles']=(8,[(0,2),(2,1),(0,3),(3,1),(0,4),(4,5),(5,0),(1,6),(6,7),(7,1)])
    return out

def graph_features(adj,I,v):
    """Entirely graph-based implementation, not a call through a view."""
    d=len(adj[v]); a=defaultdict(F); sym=lambda w:2*(len(adj[w])-2)+int(w in I)
    a[('c',)]=int(v in I)
    for h in DEGREES:
        a[('m',h)]=int(d==h)-sum((F(1,len(adj[w])) for w in adj[v] if len(adj[w])==h),F(0))
        a[('q',h)]=int(d==h and v in I)-sum((F(int(w in I),len(adj[w])) for w in adj[v] if len(adj[w])==h),F(0))
    for h in (3,4):
        a[('J',2,h)]=int(v in I)*sum(int(w in I)*((d==2 and len(adj[w])==h)-(d==h and len(adj[w])==2)) for w in adj[v])
    for u in adj[v]:
        for w in adj[u]-{v}:
            if sym(v)==sym(w): continue
            a[('R',int(w in adj[v]),min(sym(v),sym(w)),sym(u),max(sym(v),sym(w)))]+=1 if sym(v)<sym(w) else -1
    sigma=lambda u:(len(adj[u]),int(u in I),tuple(sorted(sym(w) for w in adj[u])))
    a[('S',)+sigma(v)]+=1
    for u in adj[v]: a[('S',)+sigma(u)]-=F(1,len(adj[u]))
    R=adj[v]|{v}; outer=set().union(*(adj[u] for u in adj[v]))-R
    for w in outer:
        h=len(adj[w]); T=adj[v]&adj[w]; k=len(T); x=int(w in I); b=int(({w}|(adj[w]-R))<=I)
        a[('BR',d,h,k,'X')]+=x; a[('BR',d,h,k,'B')]+=b
        a[('BR',h,d,k,'X')]-=int(v in I)
        a[('BR',h,d,k,'B')]-=int(({v}|(adj[v]-T))<=I)
        a[('D',d,h,k,'X')]+=F(1,2)-x; a[('D',d,h,k,'B')]+=F(x,2**(h-k))-b
    pattern=(d,tuple(sorted(len(adj[w]) for w in adj[v])))
    a[('L',)+pattern]=int(v in I)-BOUNDS[pattern]
    require(all(k in INDEX for k in a),'Unknown graph feature')
    return a

def graph_view(adj,v,E):
    ordered=[v]+sorted(adj[v]); index={w:j for j,w in enumerate(ordered)}; R=set(ordered); d=len(adj[v])
    endpoints=[w for u in ordered[1:] for w in sorted(adj[u]-R)]
    stubs=[index[u] for u in ordered[1:] for w in sorted(adj[u]-R)]
    boundary=list(dict.fromkeys(endpoints)); parts=[[i for i,w in enumerate(endpoints) if w==b] for b in boundary]
    V=dict(d=d,degrees=[len(adj[u]) for u in ordered[1:]],edges=[[i,j] for i,j in combinations(range(1,d+1),2) if ordered[j] in adj[ordered[i]]],stubs=stubs,parts=parts,blocks=[[stubs[i] for i in P] for P in parts],choices=[[len(adj[w]),int(w in E),int(({w}|(adj[w]-R))<=E)] for w in boundary])
    V['masks']=sorted(sum(2**j for j in S) for S in inside_sets(V))
    return V,ordered

def check_fixtures(raw_lookup):
    reports={}; total=0
    for name,(n,edges) in fixtures().items():
        adj=[set() for _ in range(n)]
        for a,b in edges: adj[a].add(b);adj[b].add(a)
        require(all(len(a) in DEGREES for a in adj) and all(edge_allowed(len(adj[a]),len(adj[b])) for a,b in edges),'Invalid fixture')
        closed=[a|{v} for v,a in enumerate(adj)]
        family=[I for I in subsets(range(n)) if not any(N<=I for N in closed)]
        direct=defaultdict(F); exposure_average=[F(0)]*NFEAT; exposures=0; masses={}; repeated=0; records=[]
        for I in family:
            configuration_sum=defaultdict(F)
            for v in range(n):
                for key,value in graph_features(adj,I,v).items(): direct[key]+=value; configuration_sum[key]+=value
            require(all(configuration_sum[k]==0 for k in DESCRIPTORS[1:NEQ]),f'{name}: fixed-configuration transport')
        for v in range(n):
            R=adj[v]|{v}; grouped=defaultdict(list)
            for I in family: grouped[I-R].append(I)
            masssum=F(0)
            for E,completions in grouped.items():
                V,ordered=graph_view(adj,v,E); raw=clean(V)
                require(raw in raw_lookup,f'{name}: actual view missing')
                actual=sorted(sum(2**j for j,w in enumerate(ordered) if w in I) for I in completions)
                require(actual==V['masks'],f'{name}: exact graph conditional family')
                vec=feature_vector(V); mass=F(len(completions),n*len(family)); masssum+=mass
                for i,val in enumerate(vec): exposure_average[i]+=mass*val
                if len({tuple(T) for T in V['blocks']})<len(V['blocks']):repeated+=1
                if name=='C4':
                    tag='opposite_present' if E else 'opposite_absent'; masses[tag]=masses.get(tag,F(0))+mass
                records.append(dict(root=v,exterior=sorted(E),full_view_id=raw_lookup[raw],mass=str(mass)))
                exposures+=1
            require(masssum==F(1,n),f'{name}: root exposure mass')
        expected=[direct[k]/(n*len(family)) for k in DESCRIPTORS]
        require(exposure_average==expected,f'{name}: graph-to-view all-feature mean mismatch')
        require(all(x==0 for x in expected[1:NEQ]),f'{name}: mean transport nonzero')
        require(all(x>=0 for x in expected[NEQ:]),f'{name}: averaged inequality negative')
        require(expected[0]==F(sum(len(I) for I in family),n*len(family)),f'{name}: occupancy mean')
        if name=='C4':
            require(masses=={'opposite_absent':F(7,11),'opposite_present':F(4,11)} and expected[0]==F(4,11),'C4 exact exposure witness')
        reports[name]=dict(vertices=n,edges=sorted(sorted(e) for e in edges),full_sets=len(family),exposures=exposures,repeated_boundary_exposures=repeated,occupancy=str(expected[0]),all_feature_means=list(map(str,expected)),exposure_records=records,c4_masses={k:str(v) for k,v in masses.items()})
        total+=exposures
    require(total==1221,'Fixture exposure total')
    return dict(all_passed=True,fixtures=reports,total_exposures=total,transport_rows=170,inequality_rows=54)

def full_replay():
    start=time.monotonic(); rawbytes=gzip.open(SOURCE,'rb').read(); require(hashlib.sha256(rawbytes).hexdigest()==DIGEST,'Fixed grammar digest')
    source=[json.loads(line) for line in rawbytes.splitlines()]; lookup={clean(V):i for i,V in enumerate(source)}
    require(len(source)==len(lookup)==17553,'Fixed grammar cardinality')
    # Separate combination-of-subsets partition recursion regenerates every full description.
    regenerated={clean(V) for V in independent_views()}
    require(regenerated==set(lookup),'Independent grammar set differs')
    columns=defaultdict(list); keys=[]
    for i,V in enumerate(source):
        key=feature_key(V); columns[key].append(i); keys.append(key)
    ordered=sorted(columns); cid={key:i for i,key in enumerate(ordered)}
    rows=[dict(id=i,d=k[0],denominator=k[1],numerators=k[2],full_view_ids=columns[k],multiplicity=len(columns[k])) for i,k in enumerate(ordered)]
    bundle=dict(schema=[list(x) for x in DESCRIPTORS],rows=rows,full_to_column=[cid[k] for k in keys],generation_sha256=DIGEST,association_bounds={str(p):str(b) for p,b in BOUNDS.items()})
    with gzip.open(HERE/'independent_columns.json.gz','wt') as f: json.dump(bundle,f,sort_keys=True)
    fixture=check_fixtures(lookup); dump(HERE/'fixture_checks.json',fixture)
    summary=dict(all_passed=True,full_views=len(source),compressed_columns=len(rows),transport_rows=170,inequality_rows=54,independent_generation_exact=True,by_root=dict(Counter(str(V['d']) for V in source)),zero_transport_indices=[i for i in range(1,NEQ) if all(r['numerators'][i]==0 for r in rows)],zero_inequality_indices=[i-NEQ for i in range(NEQ,NFEAT) if all(r['numerators'][i]==0 for r in rows)],fixture_exposures=fixture['total_exposures'],elapsed_seconds=time.monotonic()-start,generation_sha256=DIGEST)
    dump(HERE/'full_replay.json',summary); print(json.dumps(summary,indent=2),flush=True)
    return bundle

def load_independent():
    with gzip.open(HERE/'independent_columns.json.gz','rt') as f:return json.load(f)

def compare_primary(data=None):
    if data is None:
        with gzip.open(ROOT/'columns.json.gz','rt') as f:data=json.load(f)
    independent=load_independent()
    require(data['activity']=='1' and data['q']=='1/2' and data['target']=='4/11','Primary activity/target')
    require(data['generation_sha256']==DIGEST,'Primary grammar digest')
    require(data['raw_transports']==170 and data['raw_inequalities']==54 and data['feature_length']==225,'Primary schema sizes')
    require(data['full_views']==17553 and data['columns']==1722,'Primary column dimensions')
    require(data['reversals']==[list(x) for x in REV],'Reversal schema')
    require(data['stars']==json.loads(json.dumps(STARS)),'Star schema')
    require(data['classes']==[list(x) for x in CLASSES],'Boundary class schema')
    require(data['patterns']==json.loads(json.dumps(PATTERNS)),'Association pattern schema')
    require(data['association_bounds']==[str(BOUNDS[p]) for p in PATTERNS],'Fresh degree-four association constants')
    require(data['full_to_column']==independent['full_to_column'],'Complete raw-to-column mapping')
    source=[json.loads(line) for line in gzip.open(SOURCE,'rt')]
    require(len(data['rows'])==len(independent['rows']),'Column count')
    for left,right in zip(data['rows'],independent['rows']):
        for key in ('id','d','denominator','numerators'):
            require(left[key]==right[key],f'Primary/independent {key} differs at column {right["id"]}')
        require(left['source_indices']==right['full_view_ids'],'Complete inverse mapping')
        require(left['full_views']==right['multiplicity'],'Column multiplicity')
        require(left['representative']==source[right['full_view_ids'][0]],'Primary representative not first recorded full view')
    out=dict(all_passed=True,full_views=17553,columns=len(data['rows']),complete_features_checked=17553*225,all_raw_and_inverse_mappings=True,all_multiplicities=True,generation_sha256=DIGEST)
    return out

def serialization_gate():
    out=compare_primary();dump(HERE/'primary_comparison.json',out);print(json.dumps(out,indent=2))

def verify_primal(cert,bundle=None):
    if bundle is None:bundle=load_independent()
    require(cert['activity']=='1' and cert['generation_sha256']==DIGEST,'Primal activity or universe')
    rows=bundle['rows']; averages=[F(0)]*NFEAT; mass=F(0)
    require(bool(cert['support']),'Empty primal support')
    for item in cert['support']:
        i=item['id'];weight=F(item['weight'])
        require(type(i) is int and 0<=i<len(rows),'Invalid support index')
        require(weight>=0,'Negative primal weight')
        row=rows[i];mass+=weight
        for j,numerator in enumerate(row['numerators']):averages[j]+=weight*F(numerator,row['denominator'])
    require(mass==1,'Primal normalization')
    require(all(x==0 for x in averages[1:NEQ]),'Primal transport equality')
    require(all(x>=0 for x in averages[NEQ:]),'Primal inequality sign')
    require(averages[0]==F(cert['objective']),'Primal objective mismatch')
    return dict(all_passed=True,objective=str(averages[0]),support_size=len(cert['support']),transport_means=list(map(str,averages[1:NEQ])),inequality_means=list(map(str,averages[NEQ:])))

def verify_witness(cert,bundle=None):
    if bundle is None:bundle=load_independent()
    result=verify_primal(cert,bundle)
    y=[F(x) for x in cert['equality_dual']]; mu=[F(x) for x in cert['inequality_dual']]
    require(len(y)==171 and len(mu)==54,'Dual dimensions')
    require(y[0]==F(cert['objective']),'Primal-dual objective mismatch')
    require(all(x>=0 for x in mu),'Negative inequality multiplier')
    scale=lcm(*(x.denominator for x in y+mu)); coeff=[int(x*scale) for x in y+mu]
    slacks=[]
    for row in bundle['rows']:
        n=row['numerators'];D=row['denominator']
        value=scale*n[0]-coeff[0]*D-sum(coeff[j]*n[j] for j in range(1,225))
        require(value>=0,f'Negative dual slack in column {row["id"]}')
        slacks.append(F(value,scale*D))
    full=[slacks[j] for j in bundle['full_to_column']]
    require(len(full)==17553 and all(s>=0 for s in full),'Full-view slack coverage')
    result.update(dual_objective=str(y[0]),all_inequality_multipliers_nonnegative=True,compressed_slacks=len(slacks),full_view_slacks=len(full),tight_columns=sum(x==0 for x in slacks),tight_full_views=sum(x==0 for x in full),minimum_slack=str(min(full)))
    return result

def validation_gate():
    comparison=compare_primary()
    fixture=json.loads((HERE/'fixture_checks.json').read_text()); require(fixture['all_passed'] and fixture['total_exposures']==1221,'Fixture gate')
    c4=json.loads((ROOT/'c4_primal.json').read_text()); c4check=verify_primal(c4)
    require(F(c4check['objective'])==F(4,11),'Actual C4 objective')
    result=dict(all_passed=True,comparison=comparison,fixtures=10,exposures=1221,c4=c4check,optimizer_imported=False)
    dump(HERE/'validation.json',result);print(json.dumps({k:v for k,v in result.items() if k!='c4'},indent=2))

def witness_gate():
    compare_primary()
    cert=json.loads((ROOT/'certificate.json').read_text()); out=verify_witness(cert)
    dump(HERE/'certificate_verification.json',out); print(json.dumps({k:v for k,v in out.items() if k not in ('transport_means','inequality_means')},indent=2))

def corruption_controls():
    with gzip.open(ROOT/'columns.json.gz','rt') as f: data=json.load(f)
    c4=json.loads((ROOT/'c4_primal.json').read_text()); outcomes=[]
    def rejects(name,fn):
        try: fn()
        except (RuntimeError,ValueError,KeyError,IndexError) as e:outcomes.append(dict(control=name,rejected=True,error=str(e)))
        else:raise RuntimeError('Corruption accepted: '+name)
    bad=copy.deepcopy(data);bad['rows'][0]['numerators'][0]+=1
    rejects('altered full-column numerator',lambda:compare_primary(bad))
    badmap=copy.deepcopy(data);badmap['full_to_column'][0]=(badmap['full_to_column'][0]+1)%1722
    rejects('altered full-view map',lambda:compare_primary(badmap))
    badmult=copy.deepcopy(data);badmult['rows'][0]['full_views']+=1
    rejects('altered multiplicity',lambda:compare_primary(badmult))
    badsupport=copy.deepcopy(c4);badsupport['support'][0]['weight']='-1'
    rejects('negative primal weight',lambda:verify_primal(badsupport))
    badobjective=copy.deepcopy(c4);badobjective['objective']='0'
    rejects('wrong primal objective',lambda:verify_primal(badobjective))
    source=[json.loads(line) for line in gzip.open(SOURCE,'rt')]
    V=copy.deepcopy(next(v for v in source if any(len(T)==h for T,(h,x,b) in zip(v['blocks'],v['choices']))))
    j=next(j for j,(T,(h,x,b)) in enumerate(zip(V['blocks'],V['choices'])) if len(T)==h)
    V['choices'][j][1:]=[1,0]
    rejects('exhausted boundary b!=x',lambda:feature_vector(V))
    maskbad=copy.deepcopy(source[0]);maskbad['masks']=maskbad['masks'][1:]
    rejects('missing feasible inside mask',lambda:feature_vector(maskbad))
    if (ROOT/'certificate.json').exists():
        cert=json.loads((ROOT/'certificate.json').read_text())
        baddim=copy.deepcopy(cert);baddim['equality_dual'].pop()
        rejects('missing transport multiplier',lambda:verify_witness(baddim))
        badsign=copy.deepcopy(cert);badsign['inequality_dual'][0]='-1'
        rejects('negative inequality multiplier',lambda:verify_witness(badsign))
        badbound=copy.deepcopy(cert);badbound['equality_dual'][0]=str(F(cert['objective'])+1000000);badbound['objective']=badbound['equality_dual'][0]
        rejects('forged high dual bound',lambda:verify_witness(badbound))
        # Change a nonzero dual coordinate while preserving objective/multiplier sign;
        # the preceding controls independently test dimensions and sign rejection.
    out=dict(all_passed=True,controls=outcomes);dump(HERE/'corruption_controls.json',out);print(json.dumps(out,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['replay','compare','validate','witness','controls']);args=parser.parse_args()
    {'replay':full_replay,'compare':serialization_gate,'validate':validation_gate,'witness':witness_gate,'controls':corruption_controls}[args.command]()
