#!/usr/bin/env python3
"""Independent set-based degree-four weighted features and exact audits.
No primary evaluator, numerical library, or optimizer is imported.
Structural routines descend from the previously independent point audit;
all activity polynomial sums, clearing, fixture identities and signs are new.
"""
from __future__ import annotations
if not __debug__: raise RuntimeError('Independent exact audit refuses Python -O/-OO')
from fractions import Fraction as F
from collections import Counter,defaultdict
from itertools import combinations,combinations_with_replacement,product
from math import comb,lcm
from pathlib import Path
import argparse,gzip,hashlib,json,time,copy
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DIGEST='af0020c5efce7b021f901e798013b34aa7967c019c4e5c98b8352e0e2985f021'
DEGREES=(2,3,4)
def require(ok,message):
    if not ok: raise RuntimeError(message)
def dump(path,obj): Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def clean(v): return json.dumps(v,sort_keys=True,separators=(',',':'))
def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return tuple(p or [0])
def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([x*c for x in a])
def mul(a,b):
    p=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):p[i+j]+=x*y
    return trim(p)
def shift(p,k):return trim((0,)*k+tuple(p))
def oneplus(k):return tuple(comb(k,i) for i in range(k+1))
def subtract(a,b):return add(a,scale(b,-1))
def bern(p,M):
    require(M>=len(p)-1,'Insufficient Bernstein degree')
    return tuple(sum((F(x*comb(i,k),comb(M,k)) for k,x in enumerate(p) if x and k<=i),F(0)) for i in range(M+1))
def sigdegree(s):return 2+s//2
def edge_allowed(a,b):return a==2 or b==2
REV=tuple((t,a,b,c) for t in (0,1) for a in range(6) for b in range(6) for c in range(a+1,6) if edge_allowed(sigdegree(a),sigdegree(b)) and edge_allowed(sigdegree(b),sigdegree(c)) and (not t or edge_allowed(sigdegree(a),sigdegree(c))))
STARS=tuple((d,x,T) for d in DEGREES for x in (0,1) for T in combinations_with_replacement(range(6) if d==2 else range(2),d))
CLASSES=tuple((d,h,k) for d in DEGREES for h in DEGREES for k in range(1,min(d,h)+1))
PATTERNS=tuple((2,T) for T in combinations_with_replacement(DEGREES,2))+((3,(2,2,2)),(4,(2,2,2,2)))
DESCRIPTORS=(('c',),)+tuple(('m',h) for h in DEGREES)+tuple(('q',h) for h in DEGREES)+tuple(('J',2,h) for h in (3,4))+tuple(('R',)+r for r in REV)+tuple(('S',)+s for s in STARS)+tuple(('BR',d,h,k,kind) for d,h,k in CLASSES for kind in ('X','B'))
INDEX={r:i for i,r in enumerate(DESCRIPTORS)}
NEQ=171;NFEAT=225
require((len(REV),len(STARS),len(CLASSES),len(PATTERNS),len(DESCRIPTORS))==(56,60,23,8,171),'Independent schema dimensions')
def subsets(elements):
    seq=tuple(elements)
    for size in range(len(seq)+1):
        for x in combinations(seq,size):yield frozenset(x)

def association(hs):
    # Retain the root-neighbor bits. For a proper selected subset S,
    # probability = z^|S|/(1+z)^d product_{i in S}(1-(z/(1+z))^(h_i-1)).
    # Clear each term to the largest denominator needed by any proper S.
    d=len(hs); a=[h-1 for h in hs]; M=d+max(sum(a[i] for i in S) for S in subsets(range(d)) if len(S)<d)
    R=(0,)
    for S in subsets(range(d)):
        if len(S)==d:continue
        term=oneplus(M-d-sum(a[i] for i in S))
        for i in S:term=mul(term,subtract(oneplus(a[i]),(0,)*a[i]+(1,)))
        R=add(R,shift(term,len(S)))
    require(M==d+sum(a)-min(a),'Association denominator exponent')
    require(R[0]==1 and all(x>=0 and isinstance(x,int) for x in R),'Positive association numerator')
    N=shift(R,1);D=add(oneplus(M),N)
    require(D[0]==1 and all(x>=0 for x in D),'Positive association denominator')
    # Independent common-denominator product-minus-full-subset expansion.
    # T=(prod[(1+z)^h-z^h]-z^d prod[(1+z)^(h-1)-z^(h-1)])/(1+z)^sum(h).
    fullM=sum(hs); numerator=(1,);full=(0,)*d+(1,)
    for h in hs:
        numerator=mul(numerator,subtract(oneplus(h),(0,)*h+(1,)))
        full=mul(full,subtract(oneplus(h-1),(0,)*(h-1)+(1,)))
    require(mul(R,oneplus(fullM-M))==subtract(numerator,full),'Association proper-subset/product identity')
    return dict(M=M,R=R,N=N,D=D)
ASSOCIATION=tuple(association(hs) for d,hs in PATTERNS)

def validate_structure(V):
    d=V['d'];deg=V['degrees'];nb=set(range(1,d+1));edges=V['edges'];parts=V['parts'];blocks=V['blocks'];choices=V['choices'];stubs=V['stubs']
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
    validate_structure(V);d=V['d'];edges=list(map(set,V['edges']));blocked=[set(range(d+1))]
    for u in range(1,d+1):
        if all(x for T,(h,x,b) in zip(V['blocks'],V['choices']) if u in T):
            clause={0,u}
            for e in edges:
                if u in e:clause|=e
            blocked.append(clause)
    blocked.extend(set(T) for T,(h,x,b) in zip(V['blocks'],V['choices']) if b)
    return tuple(S for S in subsets(range(d+1)) if not any(Q<=S for Q in blocked))
def independent_views():
    def partitions(stubs):
        def go(rem):
            if not rem:yield ();return
            first=min(rem);tail=sorted(rem-{first})
            for count in range(min(3,len(tail))+1):
                for other in combinations(tail,count):
                    part=(first,)+other
                    if len({stubs[i] for i in part})<len(part):continue
                    for rest in go(rem-set(part)):yield (part,)+rest
        yield from go(set(range(len(stubs))))
    for d in DEGREES:
        for deg in product(DEGREES,repeat=d):
            if not all(edge_allowed(d,h) for h in deg):continue
            for eset in subsets(combinations(range(1,d+1),2)):
                edges=sorted(eset)
                if any(not edge_allowed(deg[a-1],deg[b-1]) for a,b in edges):continue
                rem=[deg[u-1]-1-sum(u in e for e in edges) for u in range(1,d+1)]
                if min(rem)<0:continue
                stubs=[u for u,n in enumerate(rem,1) for _ in range(n)]
                for parts in partitions(stubs):
                    blocks=[tuple(stubs[j] for j in P) for P in parts]
                    options=[[(h,x,b) for h in DEGREES for x,b in ((0,0),(1,0),(1,1)) if h>=len(T) and all(edge_allowed(h,deg[u-1]) for u in T) and (h!=len(T) or x==b)] for T in blocks]
                    for choices in product(*options):
                        V=dict(d=d,degrees=list(deg),edges=[list(e) for e in edges],stubs=stubs,parts=[list(p) for p in parts],blocks=[list(t) for t in blocks],choices=[list(c) for c in choices])
                        V['masks']=sorted(sum(2**j for j in S) for S in inside_sets(V));yield V

def local_configuration(V,S):
    d=V['d'];deg={0:d,**{u:h for u,h in enumerate(V['degrees'],1)}};local={u:set() for u in range(d+1)}
    for u in range(1,d+1):local[0].add(u);local[u].add(0)
    for a,b in V['edges']:local[a].add(b);local[b].add(a)
    boundary=[]
    for j,(T,(h,x,b)) in enumerate(zip(V['blocks'],V['choices']),d+1):
        deg[j]=h;boundary.append((j,set(T),h,x,b))
        for u in T:local[u].add(j)
    occupied=set(S)|{j for j,T,h,x,b in boundary if x};symbols={u:2*(h-2)+int(u in occupied) for u,h in deg.items()};a=defaultdict(F)
    a[('c',)]=int(0 in S)
    for h in DEGREES:
        a[('m',h)]=int(d==h)-sum((F(1,deg[u]) for u in local[0] if deg[u]==h),F(0))
        a[('q',h)]=int(d==h and 0 in S)-sum((F(int(u in S),deg[u]) for u in local[0] if deg[u]==h),F(0))
    for h in (3,4):a[('J',2,h)]=sum(int(0 in S and u in S)*((d==2 and deg[u]==h)-(d==h and deg[u]==2)) for u in local[0])
    for u in local[0]:
        for w in local[u]-{0}:
            x,y,z=symbols[0],symbols[u],symbols[w]
            if x!=z:a[('R',int(w in local[0]),min(x,z),y,max(x,z))]+=1 if x<z else -1
    a[('S',d,int(0 in S),tuple(sorted(symbols[u] for u in local[0])))]+=1
    for u in local[0]:a[('S',deg[u],int(u in S),tuple(sorted(symbols[w] for w in local[u])))]-=F(1,deg[u])
    for j,T,h,x,b in boundary:
        k=len(T);a[('BR',d,h,k,'X')]+=x;a[('BR',d,h,k,'B')]+=b
        a[('BR',h,d,k,'X')]-=int(0 in S);a[('BR',h,d,k,'B')]-=int(({0}|(local[0]-T))<=S)
    require(all(k in INDEX for k in a),'Unknown independent feature')
    return a

def evaluate(V):
    allowed=inside_sets(V);require(sorted(sum(2**j for j in S) for S in allowed)==V['masks'],'Recorded local family mismatch')
    Z=[0]*5;polys=[[0]*5 for _ in DESCRIPTORS]
    for S in allowed:
        k=len(S);require(k<=4,'Local partition degree');Z[k]+=1
        for desc,val in local_configuration(V,S).items():
            scaled=12*val;require(scaled.denominator==1 if isinstance(scaled,F) else isinstance(scaled,int),'Scale12 clearing')
            polys[INDEX[desc]][k]+=int(scaled)
    ext=[]
    for d,h,k in CLASSES:
        choices=[(x,b) for T,(hh,x,b) in zip(V['blocks'],V['choices']) if V['d']==d and hh==h and len(T)==k]
        ext.append((len(choices),sum(x for x,b in choices),sum(b for x,b in choices)))
    p=PATTERNS.index((V['d'],tuple(sorted(V['degrees']))))
    Z=trim(Z);require(Z[0]==1 and all(x>=0 for x in Z),'Positive local partition')
    return V['d'],Z,tuple(map(trim,polys)),tuple(ext),p

def cleared(key):
    d,Z,F12,ex,p=key;O=tuple(F(x,12) for x in F12[0]);A=(1,4,6)
    target=scale(subtract(mul(A,O),mul((0,1,3),Z)),12)
    factors=list(F12[1:])
    for (d,h,k),(N,X,B) in zip(CLASSES,ex):
        factors.append(scale(mul(Z,mul(oneplus(2),(-X,N-X))),12))
        a=h-k;f=subtract((0,)*a+(X,),scale(oneplus(a),B))
        factors.append(scale(mul(Z,mul(oneplus(3-a),f)),12))
    for j,ass in enumerate(ASSOCIATION):
        factors.append(scale(subtract(mul(ass['D'],O),mul(ass['N'],Z)),12) if j==p else (0,))
    require(len(factors)==224,'Cleared support length')
    require(all(x.denominator==1 if isinstance(x,F) else isinstance(x,int) for P in [target]+factors for x in P),'Integer cleared factors')
    return tuple(map(int,target)),tuple(tuple(map(int,P)) for P in factors)

def encode_key(k):
    d,Z,F12,ex,p=k;return dict(d=d,Z=Z,features12=F12,exterior=ex,pattern=p)
def key_from_row(r):return r['d'],tuple(r['Z']),tuple(map(tuple,r['features12'])),tuple(map(tuple,r['exterior'])),r['pattern']

def source_path():
    # Read only the local copy of the views.
    candidates=[ROOT/'inputs/views.jsonl.gz',ROOT/'inputs/point/inputs/views.jsonl.gz',ROOT/'inputs/structural/views.jsonl.gz']
    for p in candidates:
        if p.exists():return p
    raise RuntimeError('View input not found')

def enumerate_all():
    start=time.monotonic();rawbytes=gzip.open(source_path(),'rb').read();require(hashlib.sha256(rawbytes).hexdigest()==DIGEST,'Fixed grammar digest')
    source=[json.loads(line) for line in rawbytes.splitlines()];lookup={clean(V):i for i,V in enumerate(source)}
    require(len(source)==len(lookup)==17553,'Fixed grammar cardinality')
    regenerated={clean(V) for V in independent_views()};require(regenerated==set(lookup),'Independent grammar differs')
    columns=defaultdict(list);keys=[]
    for i,V in enumerate(source):key=evaluate(V);columns[key].append(i);keys.append(key)
    ordered=sorted(columns);cid={k:i for i,k in enumerate(ordered)}
    rows=[dict(id=i,**encode_key(k),source_indices=columns[k],full_views=len(columns[k])) for i,k in enumerate(ordered)]
    bundle=dict(generation_sha256=DIGEST,rows=rows,full_to_column=[cid[k] for k in keys],association=ASSOCIATION,patterns=PATTERNS,classes=CLASSES)
    with gzip.open(HERE/'independent_columns.json.gz','wt') as f:json.dump(bundle,f,sort_keys=True)
    summary=dict(all_passed=True,full_views=len(source),columns=len(rows),generation_sha256=DIGEST,by_root=dict(Counter(str(V['d']) for V in source)),transport_rows=170,inequality_rows=54,zero_transport_indices=[j for j in range(170) if all(k[2][j+1]==(0,) for k in ordered)],elapsed_seconds=time.monotonic()-start)
    dump(HERE/'enumeration.json',summary);dump(HERE/'association.json',dict(all_passed=True,patterns=PATTERNS,polynomials=ASSOCIATION));print(json.dumps(summary,indent=2),flush=True)
    return bundle,lookup

def fixtures():
    out={f'C{n}':(n,[(u,(u+1)%n) for u in range(n)]) for n in (3,4,5)}
    out.update({f'K2,{k}':(k+2,[(u,w) for u in (0,1) for w in range(2,k+2)]) for k in (3,4)})
    out['two_triangles_one_root']=(5,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)])
    out['triangle_and_square_one_root']=(6,[(0,1),(1,2),(2,0),(0,3),(4,3),(4,5),(5,0)])
    out['once_subdivision_K4']=(10,[(a,w) for w,e in enumerate(combinations(range(4),2),4) for a in e])
    out['three_once_one_twice']=(7,[(u,w) for u in (0,1) for w in (2,3,4)]+[(0,5),(5,6),(6,1)])
    out['two_once_two_returning_triangles']=(8,[(0,2),(2,1),(0,3),(3,1),(0,4),(4,5),(5,0),(1,6),(6,7),(7,1)])
    return out

def graph_configuration(adj,I,v):
    d=len(adj[v]);a=defaultdict(F);sym=lambda w:2*(len(adj[w])-2)+int(w in I)
    a[('c',)]=int(v in I)
    for h in DEGREES:
        a[('m',h)]=int(d==h)-sum((F(1,len(adj[w])) for w in adj[v] if len(adj[w])==h),F(0))
        a[('q',h)]=int(d==h and v in I)-sum((F(int(w in I),len(adj[w])) for w in adj[v] if len(adj[w])==h),F(0))
    for h in (3,4):a[('J',2,h)]=int(v in I)*sum(int(w in I)*((d==2 and len(adj[w])==h)-(d==h and len(adj[w])==2)) for w in adj[v])
    for u in adj[v]:
        for w in adj[u]-{v}:
            if sym(v)!=sym(w):a[('R',int(w in adj[v]),min(sym(v),sym(w)),sym(u),max(sym(v),sym(w)))]+=1 if sym(v)<sym(w) else -1
    sigma=lambda u:(len(adj[u]),int(u in I),tuple(sorted(sym(w) for w in adj[u])))
    a[('S',)+sigma(v)]+=1
    for u in adj[v]:a[('S',)+sigma(u)]-=F(1,len(adj[u]))
    R=adj[v]|{v};outer=set().union(*(adj[u] for u in adj[v]))-R;ex=defaultdict(lambda:[0,0,0])
    for w in outer:
        h=len(adj[w]);T=adj[v]&adj[w];k=len(T);x=int(w in I);b=int(({w}|(adj[w]-R))<=I)
        a[('BR',d,h,k,'X')]+=x;a[('BR',d,h,k,'B')]+=b
        a[('BR',h,d,k,'X')]-=int(v in I);a[('BR',h,d,k,'B')]-=int(({v}|(adj[v]-T))<=I)
        ex[d,h,k][0]+=1;ex[d,h,k][1]+=x;ex[d,h,k][2]+=b
    require(all(k in INDEX for k in a),'Unknown graph descriptor')
    return a,tuple(tuple(ex[c]) for c in CLASSES),PATTERNS.index((d,tuple(sorted(len(adj[w]) for w in adj[v]))))

def graph_view(adj,v,E):
    ordered=[v]+sorted(adj[v]);index={w:j for j,w in enumerate(ordered)};R=set(ordered);d=len(adj[v])
    endpoints=[w for u in ordered[1:] for w in sorted(adj[u]-R)];stubs=[index[u] for u in ordered[1:] for w in sorted(adj[u]-R)]
    boundary=list(dict.fromkeys(endpoints));parts=[[i for i,w in enumerate(endpoints) if w==b] for b in boundary]
    V=dict(d=d,degrees=[len(adj[u]) for u in ordered[1:]],edges=[[i,j] for i,j in combinations(range(1,d+1),2) if ordered[j] in adj[ordered[i]]],stubs=stubs,parts=parts,blocks=[[stubs[i] for i in P] for P in parts],choices=[[len(adj[w]),int(w in E),int(({w}|(adj[w]-R))<=E)] for w in boundary])
    V['masks']=sorted(sum(2**j for j in S) for S in inside_sets(V));return V,ordered

def exact_graph_numerators(adj,I,v):
    a,ext,p=graph_configuration(adj,I,v);nums=[(12*a[k],) for k in DESCRIPTORS]
    for (d,h,k),(N,X,B) in zip(CLASSES,ext):
        nums.append((-X,N-X));nums.append(subtract((0,)*(h-k)+(X,),scale(oneplus(h-k),B)))
    for j,ass in enumerate(ASSOCIATION):nums.append(subtract(scale(ass['D'],int(v in I)),ass['N']) if j==p else (0,))
    return tuple(nums),ext,p

def exact_view_numerators(key):
    d,Z,F12,ext,p=key;O=tuple(F(x,12) for x in F12[0]);nums=list(F12)
    for (d,h,k),(N,X,B) in zip(CLASSES,ext):
        nums.append(mul(Z,(-X,N-X)));nums.append(mul(Z,subtract((0,)*(h-k)+(X,),scale(oneplus(h-k),B))))
    for j,ass in enumerate(ASSOCIATION):nums.append(subtract(mul(ass['D'],O),mul(ass['N'],Z)) if j==p else (0,))
    return tuple(nums)

def fixture_audit():
    start=time.monotonic();source=[json.loads(x) for x in gzip.open(source_path(),'rt')];raw_lookup={clean(V):i for i,V in enumerate(source)};reports={};total=0;allrecords=[]
    for name,(n,edges) in fixtures().items():
        adj=[set() for _ in range(n)]
        for a,b in edges:adj[a].add(b);adj[b].add(a)
        require(all(len(a) in DEGREES for a in adj) and all(edge_allowed(len(adj[a]),len(adj[b])) for a,b in edges),'Invalid fixed fixture')
        closed=[a|{v} for v,a in enumerate(adj)];family=[I for I in subsets(range(n)) if not any(N<=I for N in closed)]
        ZG=(0,)
        for I in family:ZG=add(ZG,(0,)*len(I)+(1,))
        direct=[(0,)]*NFEAT;graph_by={};weighted_size=(0,)
        for I in family:
            sums=[F(0)]*171;weighted_size=add(weighted_size,(0,)*len(I)+(len(I),))
            for v in range(n):
                nums,ext,p=exact_graph_numerators(adj,I,v);graph_by[v,I]=nums
                for j,num in enumerate(nums):direct[j]=add(direct[j],shift(num,len(I)))
                for j,num in enumerate(nums[:171]):sums[j]+=num[0]
            require(all(x==0 for x in sums[1:]),f'{name}: configurationwise transport identity')
        require(direct[0]==scale(weighted_size,12),f'{name}: direct occupancy polynomial')
        from_exposure=[(0,)]*225;exposures=0;repeat=0;c4={}
        for v in range(n):
            R=adj[v]|{v};groups=defaultdict(list)
            for I in family:groups[I-R].append(I)
            masssum=(0,)
            for E,completions in groups.items():
                V,ordered=graph_view(adj,v,E);raw=clean(V);require(raw in raw_lookup,f'{name}: uncovered actual view')
                actual=sorted(sum(2**j for j,w in enumerate(ordered) if w in I) for I in completions);require(actual==V['masks'],f'{name}: conditional masks')
                key=evaluate(V);Z=key[1];massnum=shift(Z,len(E));actualmass=(0,)
                for I in completions:actualmass=add(actualmass,(0,)*len(I)+(1,))
                require(massnum==actualmass,f'{name}: exact exposure weight polynomial');masssum=add(masssum,massnum)
                nums=exact_view_numerators(key);exact_exposure=[]
                for j,num in enumerate(nums):
                    predicted=shift(num,len(E));actualsum=(0,)
                    for I in completions:actualsum=add(actualsum,shift(graph_by[v,I][j],len(I)))
                    require(predicted==actualsum,f'{name}: complete feature polynomial {j} at root {v}')
                    from_exposure[j]=add(from_exposure[j],predicted)
                repeat+=len({tuple(T) for T in V['blocks']})<len(V['blocks'])
                if name=='C4':c4['present' if E else 'absent']=(Z,tuple(F(x,12) for x in key[2][0]))
                allrecords.append(dict(fixture=name,root=v,exterior=sorted(E),full_view_id=raw_lookup[raw],mass_numerator=massnum,mass_denominator_scale=n,global_partition=ZG));exposures+=1
            require(masssum==ZG,f'{name}: root exposure normalization polynomial')
        require(direct==from_exposure,f'{name}: averaged complete polynomial features')
        require(all(x==(0,) for x in direct[1:171]),f'{name}: weighted mean transports')
        # Independent exact sign diagnostics on the cleared averaged inequalities.
        # Nonnegative power coefficients suffice on the whole nonnegative axis;
        # use no sampled activities or fitting here.
        power_nonnegative=[all(x>=0 for x in p) for p in direct[171:]]
        if name=='C4':
            require(c4=={'absent':((1,3,3),(F(0),F(1),F(2))),'present':((1,3),(F(0),F(1)))},'C4 complete local polynomials')
            require(ZG==(1,4,6) and direct[0]==(0,48,144),'C4 exact weighted target')
        reports[name]=dict(vertices=n,feasible_sets=len(family),exposures=exposures,repeated_boundary_exposures=repeat,partition=ZG,all_feature_identity_checks=225*exposures,transport_polynomial_identities=170,inequality_power_coefficient_nonnegative=power_nonnegative,averaged_numerators=[[str(x) for x in p] for p in direct])
        total+=exposures
    require(total==1221,'Fixed exposure total')
    report=dict(all_passed=True,fixtures=reports,total_exposures=total,complete_feature_polynomial_identities=225*total,elapsed_seconds=time.monotonic()-start)
    dump(HERE/'fixtures.json',report)
    with gzip.open(HERE/'fixture_exposures.json.gz','wt') as f:json.dump(allrecords,f,sort_keys=True)
    print(json.dumps({k:v for k,v in report.items() if k!='fixtures'},indent=2),flush=True)
    return report

SUPPORT=([f'm{k}' for k in DEGREES]+[f'q{k}' for k in DEGREES]+['J23','J24']+[f'R_{t}_{a}_{b}_{c}' for t,a,b,c in REV]+[f'S_{d}_{x}_'+','.join(map(str,T)) for d,x,T in STARS]+[f'BR_{kind}_{d}_{h}_{k}' for d,h,k in CLASSES for kind in ('X','B')]+[f'E{kind}_{d}_{h}_{k}' for d,h,k in CLASSES for kind in ('X','B')]+[f'L_{d}_'+','.join(map(str,T)) for d,T in PATTERNS])
require(len(SUPPORT)==224,'Full independent multiplier support')

def primary_key(r):
    return r['d'],tuple(r['Z']),(scale(tuple(r['O']),12),)+tuple(map(tuple,r['F'])),tuple(map(tuple,r['exterior'])),r['pattern']

def compare():
    start=time.monotonic()
    with gzip.open(HERE/'independent_columns.json.gz','rt') as f:ind=json.load(f)
    with gzip.open(ROOT/'symbolic_columns.json.gz','rt') as f:pri=json.load(f)
    require(pri['generation_sha256']==ind['generation_sha256']==DIGEST,'Comparison scope digest')
    require(pri['equality_names']+pri['inequality_names']==SUPPORT,'Complete row order')
    require(clean(pri['association'])==clean(ind['association']),'Independent association coefficient comparison')
    imap=ind['full_to_column'];pmap=pri['full_to_column'];require(len(imap)==len(pmap)==17553,'Raw mapping cardinality')
    groups=defaultdict(list);cleared_hash=hashlib.sha256();counts=0
    for i,(a,b) in enumerate(zip(imap,pmap)):
        ik=key_from_row(ind['rows'][a]);pk=primary_key(pri['rows'][b]);require(ik==pk,f'Full complete polynomial key mismatch {i}');groups[b].append(i);counts+=225
    for row in pri['rows']:
        require(groups[row['id']]==row['source_indices'] and row['full_views']==len(groups[row['id']]),'Primary raw inverse/multiplicity')
        target,factors=cleared(primary_key(row));cleared_hash.update((clean([target,factors])+'\n').encode())
    require(len(set(primary_key(r) for r in pri['rows']))==len(pri['rows'])==len(ind['rows']),'Complete symbolic compression')
    c4=[]
    for r in pri['rows']:
        V=r['representative']
        if V['d']==2 and V['degrees']==[2,2] and not V['edges'] and V['blocks']==[[1,2]] and V['choices'][0][0]==2:c4.append(r['id'])
    require(sorted(c4)==sorted(pri['c4_columns']) and len(c4)==2,'C4 calibration columns')
    out=dict(all_passed=True,full_views=len(pmap),symbolic_columns=len(pri['rows']),complete_feature_comparisons=counts,all_raw_mappings_and_multiplicities=True,independent_cleared_polynomial_sha256=cleared_hash.hexdigest(),c4_columns=c4,elapsed_seconds=time.monotonic()-start)
    dump(HERE/'comparison.json',out);print(json.dumps(out,indent=2),flush=True);return out

def basis_to_power(values):
    d=len(values)-1;out=[F(0)]*(d+1)
    for i,value in enumerate(values):
        for k in range(d-i+1):out[i+k]+=value*comb(d,i)*comb(d-i,k)*((-1)**k)
    return trim(out)

def verify_certificate(certpath,output=None,full=True):
    start=time.monotonic();cert=json.loads(Path(certpath).read_text())
    require(cert['schema_version']=='degree_four_weighted_interval_v1','Certificate schema version')
    require(cert['interval']==['0','1'],'Certificate interval')
    require(cert['clearing_factor']=='12*A*Z_C','Certificate clearing factor')
    stage=cert['stage'];require(stage in ('global_d2','global_d4'),'Certificate stage')
    d,M=(2,16) if stage=='global_d2' else (4,20)
    require(cert['degree']==d and cert['bernstein_degree']==M,'Fixed degrees')
    require(cert['support']==SUPPORT and cert['generation_sha256']==DIGEST,'Certificate support/scope')
    vals=cert['bernstein_polynomials'];require(len(vals)==224 and all(len(p)==d+1 for p in vals),'Complete certificate dimensions')
    vals=[tuple(map(F,p)) for p in vals];require(all(v>=0 for p in vals[170:] for v in p),'Negative multiplier Bernstein coefficient')
    polynomials=list(map(basis_to_power,vals));scale_den=lcm(*(x.denominator for p in polynomials for x in p));ips=[tuple(int(x*scale_den) for x in p) for p in polynomials]
    with gzip.open(HERE/'independent_columns.json.gz','rt') as f:data=json.load(f)
    minimum=None;tight=[];slack_digest=hashlib.sha256();checked=0;c4=[]
    for row in data['rows']:
        key=key_from_row(row);target,factors=cleared(key);slack=scale(target,scale_den)
        for p,q in zip(ips,factors):slack=subtract(slack,mul(p,q))
        beta=bern(slack,M)
        require(all(x>=0 for x in beta),f'Negative exact Bernstein slack independent column {row["id"]}: '+str([(i,str(x)) for i,x in enumerate(beta) if x<0]))
        if minimum is None or min(beta)<minimum:minimum=min(beta)
        if all(x==0 for x in beta):tight.append(row['id'])
        if key[0]==2 and key[1] in ((1,3,3),(1,3)) and key[2][0] in ((0,12,24),(0,12)):
            # Actual C4 identification is validated below from its full view, not only these polynomials.
            pass
        slack_digest.update(('|'.join(str(x/F(scale_den)) for x in beta)+'\n').encode());checked+=len(beta)
    if full:
        raw=gzip.open(source_path(),'rb').read();require(hashlib.sha256(raw).hexdigest()==DIGEST,'Replay raw grammar digest')
        source=[json.loads(x) for x in raw.splitlines()];require(len(source)==len(data['full_to_column'])==17553,'Replay full scope')
        inverse=defaultdict(list)
        for i,V in enumerate(source):
            key=evaluate(V);cid=data['full_to_column'][i];require(key==key_from_row(data['rows'][cid]),f'Regenerated full weighted feature mismatch {i}');inverse[cid].append(i)
            if V['d']==2 and V['degrees']==[2,2] and V['edges']==[] and V['blocks']==[[1,2]] and V['choices'][0][0]==2:
                require(cid in tight,'C4 polynomial slack is not identically zero');c4.append(i)
        require(len(c4)==2,'Actual C4 calibration full views')
        require({clean(V) for V in independent_views()}==set(map(clean,source)),'Full independent grammar regeneration')
        for r in data['rows']:require(inverse[r['id']]==r['source_indices'] and len(inverse[r['id']])==r['full_views'],'Replay full inverse multiplicity')
    out=dict(all_passed=True,stage=stage,degree=d,bernstein_degree=M,full_views=17553 if full else None,symbolic_columns=len(data['rows']),checked_slack_coefficients=checked,checked_multiplier_coefficients=54*(d+1),minimum_scaled_bernstein=str(minimum),common_coefficient_denominator=str(scale_den),tight_columns=len(tight),tight_full_views=sum(data['rows'][i]['full_views'] for i in tight),c4_calibration_full_views=c4,slack_bernstein_sha256=slack_digest.hexdigest(),elapsed_seconds=time.monotonic()-start,optimizer_free=True)
    if output:dump(output,out)
    print(json.dumps(out,indent=2),flush=True);return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['enumerate','fixtures','compare','verify']);p.add_argument('--certificate');p.add_argument('--output');a=p.parse_args()
    if a.command=='enumerate':enumerate_all()
    elif a.command=='fixtures':fixture_audit()
    elif a.command=='compare':compare()
    elif a.command=='verify':verify_certificate(a.certificate,a.output)
