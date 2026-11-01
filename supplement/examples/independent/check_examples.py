#!/usr/bin/env python3
"""Independent exact explanatory-example checks.

No project generators, symbolic engines, or optimizers are imported.
K2,3 admissibility is checked by direct closed-neighborhood containment
for every subset. All polynomial arithmetic uses integers/Fraction.
Finite examples validate the formulas; the subdivision-family statement
is also proved by the elementary argument in RESULTS.md.
"""
from __future__ import annotations
if not __debug__:
    raise RuntimeError('Exact example audit rejects optimized Python')
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import Counter
import hashlib,json

HERE=Path(__file__).resolve().parent

def require(ok,message):
    if not ok:raise RuntimeError(message)

def powerset(vertices):
    v=tuple(vertices)
    for k in range(len(v)+1):
        for subset in combinations(v,k):yield frozenset(subset)

def make_graph(vertices,edges):
    graph={v:set() for v in vertices}
    for u,v in edges:
        require(u!=v and u in graph and v in graph,'Invalid edge')
        require(v not in graph[u],'Duplicate edge')
        graph[u].add(v);graph[v].add(u)
    return graph

def feasible(graph,chosen):
    return all(not (set(neighbors)|{v})<=chosen for v,neighbors in graph.items())

def poly_from_sets(sets,stat=None):
    stat=stat or (lambda subset:1)
    values=Counter()
    for subset in sets:values[len(subset)]+=stat(subset)
    return trim([values[k] for k in range(max(values,default=0)+1)])

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p

def add(*polys):
    out=[0]*max(map(len,polys))
    for p in polys:
        for i,x in enumerate(p):out[i]+=x
    return trim(out)

def scale(p,k):return trim([k*x for x in p])
def shift(p,k):return [0]*k+list(p)
def mul(p,q):
    out=[0]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):out[i+j]+=x*y
    return trim(out)

def evaluate(p,x):return sum(c*x**i for i,c in enumerate(p))
def derivative(p):return trim([i*p[i] for i in range(1,len(p))] or [0])
def serial_poly(p):return [str(x) if isinstance(x,Q) and x.denominator!=1 else int(x) for x in p]
def serial_sets(sets):return [sorted(S) for S in sorted(sets,key=lambda S:(len(S),tuple(sorted(S))))]

def main():
    vertices=('a','a_prime','b','c','d')
    edges=tuple((a,b) for a in vertices[:2] for b in vertices[2:])
    G=make_graph(vertices,edges)
    allsets=list(powerset(vertices));F=[S for S in allsets if feasible(G,S)]
    Z=poly_from_sets(F)
    require(Z==[1,5,10,7],'Global K2,3 polynomial mismatch')
    marginals={v:poly_from_sets(F,lambda S,v=v:int(v in S)) for v in vertices}
    require(all(marginals[v]==[0,1,4,3] for v in vertices[:2]),'Degree-three marginal mismatch')
    require(all(marginals[v]==[0,1,4,5] for v in vertices[2:]),'Degree-two marginal mismatch')
    # Direct dominating-set enumeration is a second predicate.
    dominating=[D for D in allsets if all(v in D or bool(set(G[v])&D) for v in vertices)]
    Dpoly=poly_from_sets(dominating)
    require(Dpoly==[0,0,7,10,5,1],'Domination polynomial mismatch')
    require({frozenset(set(vertices)-D) for D in dominating}==set(F),'Complement bijection mismatch')
    q2=lambda v,S:Q(int(len(G[v])==2 and v in S))-sum((Q(int(u in S),2) for u in G[v] if len(G[u])==2),Q(0))
    exposures={}
    for root in ('b','a'):
        inside=tuple(v for v in vertices if v==root or v in G[root])
        exterior=tuple(v for v in vertices if v not in inside)
        rows=[];recovered_Z=[0];recovered_O=[0];recovered_q2=[0]
        for E in powerset(exterior):
            local=[S for S in powerset(inside) if feasible(G,S|E)]
            require(all(S|E in F for S in local),'Conditional completion failure')
            Zc=poly_from_sets(local);Oc=poly_from_sets(local,lambda S:int(root in S))
            q2num=poly_from_sets(local,lambda S:q2(root,S|E))
            mass_numerator=shift(Zc,len(E));root_numerator=shift(Oc,len(E))
            recovered_Z=add(recovered_Z,mass_numerator)
            recovered_O=add(recovered_O,root_numerator)
            recovered_q2=add(recovered_q2,shift(q2num,len(E)))
            rows.append(dict(exterior=sorted(E),inside=list(inside),feasible_inside_sets=serial_sets(local),
                Z=Zc,O=Oc,q2_numerator=serial_poly(q2num),q2_denominator=Zc,
                fixed_root_mass_numerator=mass_numerator,fixed_root_mass_denominator=Z,
                uniform_root_pair_mass_denominator=scale(Z,5),
                fixed_root_mass_at_one=str(Q(sum(mass_numerator),sum(Z))),
                uniform_root_pair_mass_at_one=str(Q(sum(mass_numerator),5*sum(Z))),
                root_conditional_occupancy_at_one=str(Q(sum(Oc),sum(Zc))),
                conditional_q2_at_one=str(Q(sum(q2num),sum(Zc)))))
        require(recovered_Z==Z,'True exposure masses do not normalize')
        require(recovered_O==marginals[root],'Exposure marginal reconstruction failed')
        global_q2num=poly_from_sets(F,lambda S:q2(root,S))
        require(recovered_q2==global_q2num,'Transport reconstruction failed')
        exposures[root]=dict(degree=len(G[root]),rows=rows,
                global_q2_numerator=serial_poly(global_q2num),global_q2_denominator=Z,
                global_q2_at_one=str(Q(sum(global_q2num),sum(Z))))
    # Extra selected transport row for E={c}, using every ordered walk.
    E=frozenset({'c'});root='b';inside=('a','a_prime','b')
    local=[S for S in powerset(inside) if feasible(G,S|E)]
    def reversal(S):
        selected=S|E
        symbol=lambda v:2*(len(G[v])-2)+int(v in selected)
        total=0
        for u in G[root]:
            for w in G[u]:
                if w==root or w in G[root]:continue
                triple=(symbol(root),symbol(u),symbol(w))
                total+=int(triple==(0,2,1))-int(triple==(1,2,0))
        return total
    feature_nums={
      'm2':poly_from_sets(local,lambda S:Q(1)),
      'm3':poly_from_sets(local,lambda S:-Q(2,3)),
      'q2':poly_from_sets(local,lambda S:Q(int(root in S))),
      'q3':poly_from_sets(local,lambda S:-sum((Q(int(u in S),3) for u in G[root]),Q(0))),
      'j':poly_from_sets(local,lambda S:sum(int(root in S and u in S) for u in G[root])),
      'R_0_0_2_1':poly_from_sets(local,reversal)}
    require(feature_nums['m2']==[1,3,2] and feature_nums['m3']==[Q(-2,3),-2,Q(-4,3)],'Degree transport formulas mismatch')
    require(feature_nums['q2']==[0,1,2] and feature_nums['q3']==[0,Q(-2,3),Q(-2,3)],'Occupancy transport formulas mismatch')
    require(feature_nums['j']==[0,0,2] and feature_nums['R_0_0_2_1']==[2,0,-2],'Edge/reversal formulas mismatch')
    worked_feature_data=dict(root=root,exterior=sorted(E),Z=[1,3,2],
        feature_numerators={k:serial_poly(v) for k,v in feature_nums.items()},
        per_inside_reversal=[dict(inside=sorted(S),statistic=reversal(S)) for S in sorted(local,key=lambda S:(len(S),tuple(sorted(S))))])
    expected_b=[([1,3,3],[0,1,2]),([1,3,2],[0,1,2]),([1,3,2],[0,1,2]),([1,3],[0,1])]
    require([(r['Z'],r['O']) for r in exposures['b']['rows']]==expected_b,'Degree-two exposure formulas mismatch')
    require([r['fixed_root_mass_at_one'] for r in exposures['b']['rows']]==['7/23','6/23','6/23','4/23'],'Degree-two exposure masses mismatch')
    require([(r['Z'],r['O']) for r in exposures['a']['rows']]==[([1,4,6,4],[0,1,3,3]),([1,4,3],[0,1])],'Degree-three exposures mismatch')
    # Transport is an identity for each feasible set, not only at z=1.
    require(all(sum((q2(v,S) for v in vertices),Q(0))==0 for S in F),'q2 pointwise root-sum identity failed')
    q2rootnum=add(*(poly_from_sets(F,lambda S,v=v:q2(v,S)) for v in vertices))
    require(q2rootnum==[0],'Polynomial average transport failure')
    occupancy_numerator=shift(derivative(Z),1)
    require(occupancy_numerator==add(*(marginals.values())),'Global mean identity mismatch')
    A=[1,4,6];c4_num=[0,1,3]
    diff_degree3=add(mul(marginals['a'],A),scale(mul(c4_num,Z),-1))
    diff_degree2=add(mul(marginals['b'],A),scale(mul(c4_num,Z),-1))
    diff_mean=add(mul(occupancy_numerator,A),scale(mul(c4_num,Z),-5))
    require(diff_degree3==[0,0,0,0,-1,-3],'Degree-three comparison mismatch')
    require(diff_degree2==[0,0,0,2,7,9],'Degree-two comparison mismatch')
    require(diff_mean==[0,0,0,6,19,21],'Mean comparison mismatch')
    # Independent finite graph constructions for the stated general family.
    originals={
      'K4':make_graph(tuple('0123'),tuple(combinations(tuple('0123'),2))),
      'K3,3':make_graph(tuple('abcdef'),tuple((u,v) for u in 'abc' for v in 'def')),
      'triangular_prism':make_graph(tuple('012345'),(('0','1'),('1','2'),('2','0'),('3','4'),('4','5'),('5','3'),('0','3'),('1','4'),('2','5'))),
      'cube':make_graph(tuple(f'{i:03b}' for i in range(8)),tuple((f'{i:03b}',f'{j:03b}') for i in range(8) for j in range(i+1,8) if (i^j).bit_count()==1))}
    family=[]
    for name,H in originals.items():
        require(all(len(N)==3 for N in H.values()),'Original is not cubic')
        es=sorted({tuple(sorted((u,v))) for u,N in H.items() for v in N})
        old=list(H);new=[f'edge_{i}' for i in range(len(es))]
        subdiv_edges=tuple(edge for i,(u,v) in enumerate(es) for edge in ((u,new[i]),(v,new[i])))
        J=make_graph(old+new,subdiv_edges);counts=Counter(map(len,J.values()))
        N=len(H)
        require(len(es)==3*N//2 and len(J)==5*N//2,'Subdivision size formula failed')
        require(counts=={2:3*N//2,3:N},'Subdivision degree counts failed')
        require(all((u in old)!=(v in old) for u,v in subdiv_edges),'Bipartition failure')
        require(Q(counts[2],len(J))==Q(3,5)>Q(7,187),'Profile comparison failed')
        # Every degree-two inserted vertex forces both incident edges in a 2-factor.
        forced={tuple(sorted((w,v))) for w in new for v in J[w]}
        forced_degrees=Counter(v for edge in forced for v in edge)
        require(all(forced_degrees[v]==3 for v in old),'Forced original degree is not three')
        row=dict(original=name,N=N,original_edges=len(es),subdivided_order=len(J),
                 n2=counts[2],n3=counts[3],n2_fraction='3/5',bipartition_sizes=[len(old),len(new)],
                 no_two_factor_forced_original_degrees=sorted(set(forced_degrees[v] for v in old)))
        if name=='K4':
            checked=0;factors=[]
            for selected_edges in powerset(subdiv_edges):
                checked+=1;degrees=Counter(v for e in selected_edges for v in e)
                if all(degrees[v]==2 for v in J):factors.append(selected_edges)
            require(checked==4096 and not factors,'Subdivision K4 exhaustive 2-factor control failed')
            row.update(edge_subsets_exhausted=checked,two_factors_found=len(factors))
        family.append(row)
    result=dict(status='PASS',method='Direct ordinary-set exhaustive closed-neighborhood predicates; integer/Fraction arithmetic only',
       K23=dict(vertices=list(vertices),edges=edges,subsets_exhausted=len(allsets),admissible_subsets=len(F),
          feasible_sets=serial_sets(F),Z=Z,D=Dpoly,marginal_numerators=marginals,worked_E_c_transports=worked_feature_data,
          marginals_at_one={v:str(Q(sum(p),sum(Z))) for v,p in marginals.items()},
          alpha_numerator=occupancy_numerator,alpha_denominator=scale(Z,5),alpha_at_one='2/5',
          domination_mean_density_at_one='3/5',exposures=exposures,
          q2_root_sum_zero_for_each_feasible_set=True,q2_root_sum_polynomial=serial_poly(q2rootnum),
          generic_difference_numerators=dict(degree3_minus_c4=diff_degree3,degree2_minus_c4=diff_degree2,alpha_minus_c4=diff_mean),
          generic_difference_denominators=dict(degree3_minus_c4=mul(A,Z),degree2_minus_c4=mul(A,Z),alpha_minus_c4=scale(mul(A,Z),5))),
       subdivision_controls=family,
       general_subdivision_argument='Every original cubic graph has N vertices and 3N/2 edges. Full once-subdivision has N degree-three originals and 3N/2 degree-two inserted vertices. A spanning 2-factor would have to use both incident edges at every inserted vertex, forcing all three edges at every original vertex, impossible. Equivalently the bipartition has unequal sizes N and 3N/2, while every cycle uses equally many of both parts.',
       limits='Explanatory checks only; not a new theorem, antecedent clearance, or extension of the delivered manuscript.')
    (HERE/'RESULTS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('PASS: all 32 K2,3 subsets, both complement/domination predicates, per-degree marginals, all six fixed-root exposure cases, true masses, and q2 transport.')
    print('PASS: exact generic-z differences and once-subdivided cubic counts/forced-degree contradiction; all 4096 edge subsets for subdivided K4.')

if __name__=='__main__':main()
