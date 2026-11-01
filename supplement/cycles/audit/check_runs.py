#!/usr/bin/env python3
"""Independent exact check on the fixed cycle orders 3..40.

No primary checker is imported. Coefficients come from rooted-zero runs:
m*[z^k]Z_n = n*[x^k](1+x+x^2)^m, m=n-k>=1.
"""
import argparse
import json
from pathlib import Path

if not __debug__:
    raise SystemExit('Refusing optimized Python: exact-check assertions must remain enabled.')

LOW, HIGH = 3, 40
HERE = Path(__file__).resolve().parent


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(*polys):
    c = [0] * max(map(len, polys))
    for p in polys:
        for j, v in enumerate(p):
            c[j] += v
    return trim(c)


def shift(p, degree, factor=1):
    return trim([0]*degree + [factor*v for v in p])


def get(a, j):
    return a[j] if 0 <= j < len(a) else 0


def runs(n):
    # powrows[m] is the coefficient vector of (1+x+x^2)^m.
    powrows = [[1]]
    for m in range(1, n+1):
        old = powrows[-1]
        row = [0] * (len(old)+2)
        for j, v in enumerate(old):
            row[j] += v
            row[j+1] += v
            row[j+2] += v
        powrows.append(row)
    out = []
    for k in range(n):
        m = n-k
        numerator = n*get(powrows[m], k)
        q, rem = divmod(numerator, m)
        assert rem == 0, ('rooted-zero integrality', n, k, numerator, m)
        out.append(q)
    return trim(out)


def gap(n, a):
    # Coefficient-level differentiation, independent of polynomial products.
    return trim([(j+1)*get(a, j+1)+(4*j-n)*get(a,j)
                 +(6*(j-1)-3*n)*get(a,j-1)
                 for j in range(len(a)+1)])


def generate():
    z = {1:[1], 2:[1,2]}
    z.update({n:runs(n) for n in range(LOW,HIGH+1)})
    h = {n:gap(n,z[n]) for n in range(LOW,HIGH+1)}
    for n in range(4,HIGH+1):
        assert z[n] == add(z[n-1],shift(z[n-2],1),shift(z[n-3],2)), ('Z recurrence',n)
    for n in range(LOW,HIGH+1):
        assert all(v>=0 for v in h[n]), ('negative coefficient',n,h[n])
        assert (h[n]==[0]) == (n==4), ('equality classification',n,h[n])
        if n>=5:
            assert all(get(h[n],j)==0 for j in range(3)), ('lower powers',n)
            assert get(h[n],3)==n, ('cubic coefficient',n)
    for n in range(7,HIGH+1):
        forced = add(h[n-1],shift(h[n-2],1),shift(h[n-3],2),
                     shift(z[n-5],3),shift(z[n-5],4,3),
                     shift(z[n-6],4,2),shift(z[n-6],5,3))
        assert h[n] == forced, ('positive-forcing recurrence',n,h[n],forced)
    data = {
        'method':'Rooted-zero run count; integer arithmetic only',
        'fixed_orders':[LOW,HIGH],
        'formal_initial_traces':{'1':z[1],'2':z[2]},
        'records':[{'n':n,'Z':z[n],'H':h[n]} for n in range(LOW,HIGH+1)],
        'checks':{
            'run_integrality':sum(range(LOW,HIGH+1)),
            'Z_recurrence':HIGH-4+1,
            'H_nonnegativity_and_equality':HIGH-LOW+1,
            'positive_forcing':HIGH-7+1,
            'leading_cubic_coefficient':HIGH-5+1
        }
    }
    return data


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=HERE/'RUN_CHECKS.json')
    p.add_argument('--compare',type=Path,help='Read primary JSON only after independent generation and checks')
    args = p.parse_args()
    data = generate()
    args.output.write_text(json.dumps(data,indent=2)+'\n')
    print('PASS: independent rooted-zero exact counting for every n=3..40.')
    print('PASS: every divisibility, Z recurrence, coefficient sign, equality, leading-power, and forcing identity check.')
    print(json.dumps(data['checks'],sort_keys=True))
    print('Saved independently generated records:',args.output)
    if args.compare is not None:
        primary = json.loads(args.compare.read_text())
        assert primary['coefficient_orders']==list(range(LOW,HIGH+1)), 'Primary range differs from the fixed range'
        for record in data['records']:
            n = record['n']
            supplied = primary['all_exact_polynomials'][str(n)]
            assert record['Z']==supplied['Z'], ('primary Z disagreement',n)
            assert record['H']==supplied['H'], ('primary H disagreement',n)
        for n,formal in data['formal_initial_traces'].items():
            assert formal==primary['formal_initial_traces'][n], ('formal initial trace disagreement',n)
        print('PASS: all 38 independently generated Z and H vectors agree with primary EXACT_CHECKS.json.')

if __name__=='__main__':
    main()
