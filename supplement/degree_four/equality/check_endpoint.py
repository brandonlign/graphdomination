#!/usr/bin/env python3
"""Endpoint-only evidence from fixed stored views and the sealed candidate."""
if not __debug__:raise RuntimeError('Endpoint audit refuses Python -O/-OO')
from fractions import Fraction as Q
import gzip,json,time
import check_equality as e

def main():
    t=time.monotonic();d=e.read('primary_columns.json.gz');cert=e.read('certificate.json')
    views=[json.loads(line) for line in gzip.open(e.INP/'views.jsonl.gz','rt')]
    p=[Q(a[-1]) for a in cert['bernstein_polynomials']];assoc=e.associations(d['patterns'])
    endpoint=json.loads((e.BASE/'ENDPOINT_SUPPORT.json').read_text());end=set(endpoint['columns'])
    records=[]
    for i,v in enumerate(views):
        if v['d']!=4 or any(x for h,x,b in v['choices']):continue
        e.need(all(b==0 for h,x,b in v['choices']),'Empty exterior flags')
        r=e.reconstruct(v,d['support'],d['classes'],d['patterns']);col=d['full_to_column'][i]
        e.need(e.same_feature(r,d['rows'][col]),'Empty degree4 independent feature mismatch')
        target,ff=e.factors(r,d['classes'],assoc)
        H=Q(sum(target))-sum((a*sum(f) for a,f in zip(p,ff)),Q())
        e.need(H>0 and col not in end,'Degree4 empty view not strict at endpoint')
        records.append(dict(raw_index=i,column=col,view_sha256=e.sh(v),slack_at_one=str(H)))
    e.need(len(records)==370,'Degree-four empty-view fixed count')
    j=cert['support'].index('L_3_2,2,2');e.need(p[j]>0,'Degree-three endpoint multiplier not positive')
    a=assoc[d['patterns'].index([3,[2,2,2]])]
    mu=p[j]*Q(sum(a['D']),11)
    report=dict(all_passed=True,degree4_empty_views_reconstructed=len(records),
                independent_local_polynomials_checked=len(records)*172,
                minimum_degree4_empty_slack_at_one=str(min(Q(r['slack_at_one']) for r in records)),
                degree3_multiplier_polynomial_at_one=str(p[j]),degree3_actual_multiplier_at_one=str(mu),
                degree3_association_lower_bound_at_one=str(Q(sum(a['N']),sum(a['D']))),
                degree4_multiplier_polynomial_at_one=str(p[cert['support'].index('L_4_2,2,2,2')]),
                no_optimization=True,no_new_census=True,elapsed_seconds=time.monotonic()-t)
    e.dump('ENDPOINT_EMPTY_DEGREE4.json',dict(summary=report,views=records));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
