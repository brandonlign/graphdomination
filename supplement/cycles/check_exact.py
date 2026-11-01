#!/usr/bin/env python3
"""Exact checks on a fixed range of cycle orders, using the standard library.

Run with python3 check_exact.py. Polynomial coefficients are ascending.
The range checks validate implementation, not the unbounded theorem.
"""
import hashlib
import json
from pathlib import Path

if not __debug__:
    raise RuntimeError('Exact checks require assertions; Python -O is rejected.')

ROOT = Path(__file__).resolve().parent
RANGE_ENUMERATION = (3, 12)
RANGE_COEFFICIENTS = (3, 40)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(*polys):
    result = [0] * max(map(len, polys))
    for p in polys:
        for i, c in enumerate(p):
            result[i] += c
    return trim(result)


def scale(p, c):
    return trim([c * x for x in p])


def shift(p, k):
    return trim([0] * k + list(p))


def mul(p, q):
    result = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            result[i+j] += a*b
    return trim(result)


def deriv(p):
    return trim([i * p[i] for i in range(1, len(p))] or [0])


def matmul(m, n):
    return [[add(*(mul(m[i][k], n[k][j]) for k in range(3)))
             for j in range(3)] for i in range(3)]


def matrix_power(m, n):
    out = [[(int(i == j),) for j in range(3)] for i in range(3)]
    for _ in range(n):
        out = matmul(out, m)
    return out


def trace(m):
    return add(*(m[i][i] for i in range(3)))


def independent_enumeration(n):
    """Neighborhood-mask feasibility versus independent vertex domination."""
    full = (1 << n) - 1
    neighborhoods = [sum(1 << j for j in {(v-1) % n, v, (v+1) % n})
                     for v in range(n)]
    allowed_counts = [0] * (n+1)
    dominating_counts = [0] * (n+1)
    allowed_sets = set()
    dominating_sets = set()
    for mask in range(1 << n):
        if all((mask & nb) != nb for nb in neighborhoods):
            allowed_counts[mask.bit_count()] += 1
            allowed_sets.add(mask)
        # Independent predicate: every vertex is selected or has a selected neighbor.
        chosen = {v for v in range(n) if (mask >> v) & 1}
        if all(v in chosen or (v-1) % n in chosen or (v+1) % n in chosen
               for v in range(n)):
            dominating_counts[len(chosen)] += 1
            dominating_sets.add(mask)
    assert {full ^ mask for mask in allowed_sets} == dominating_sets
    assert allowed_counts == list(reversed(dominating_counts))
    return trim(allowed_counts), trim(dominating_counts)


def symbolic_forcing_check():
    """Formal coefficients in independent symbols U_j=Z_(n-j).

    Reduce the original forcing by the three recurrence substitutions,
    without generating any particular Z_n or using numerical activities.
    """
    original = {
        1: (-1, -3),
        2: (1, 2),
        3: (0, 2, 5, 3),
    }
    state = original.copy()
    history = [{str(k): list(v) for k, v in sorted(state.items())}]
    for j in (1, 2, 3):
        coeff = state.pop(j, (0,))
        # U_j = U_(j+1) + z U_(j+2) + z² U_(j+3)
        for offset in (1, 2, 3):
            target = j + offset
            state[target] = add(state.get(target, (0,)), shift(coeff, offset-1))
        state = {k: v for k, v in state.items() if v != (0,)}
        history.append({str(k): list(v) for k, v in sorted(state.items())})
    expected = {5: (0, 0, 0, 1, 3), 6: (0, 0, 0, 0, 2, 3)}
    assert state == expected, (state, expected)
    return history


def main():
    a, b = (1, 4, 6), (1, 3)
    z = {1: (1,), 2: (1, 2), 3: (1, 3, 3)}
    for n in range(4, RANGE_COEFFICIENTS[1] + 1):
        z[n] = add(z[n-1], shift(z[n-2], 1), shift(z[n-3], 2))
    # Matrix computation independently validates the initial trace convention.
    s = [[(1,), (0, 1), (0,)],
         [(1,), (0,), (0, 1)],
         [(1,), (0,), (0,)]]
    for n in (1, 2, 3):
        assert trace(matrix_power(s, n)) == z[n]
    # Direct matrix Cayley-Hamilton identity: S³-S²-zS-z²I=0.
    s2, s3 = matrix_power(s, 2), matrix_power(s, 3)
    for i in range(3):
        for j in range(3):
            assert add(s3[i][j], scale(s2[i][j], -1),
                       scale(shift(s[i][j], 1), -1),
                       (0, 0, -int(i == j))) == (0,)

    enumerations = []
    for n in range(RANGE_ENUMERATION[0], RANGE_ENUMERATION[1] + 1):
        direct_z, direct_d = independent_enumeration(n)
        assert z[n] == direct_z
        assert trace(matrix_power(s, n)) == direct_z
        enumerations.append({'n': n, 'Z': list(direct_z), 'D': list(direct_d),
                             'subsets_checked_in_each_predicate': 1 << n})

    h = {n: add(mul(a, deriv(z[n])), scale(mul(b, z[n]), -n))
         for n in range(1, RANGE_COEFFICIENTS[1] + 1)}
    bases = {3: (0, 0, 6, 9), 4: (0,), 5: (0, 0, 0, 5, 15),
             6: (0, 0, 0, 6, 30, 18)}
    for n, expected in bases.items():
        assert h[n] == expected
    # Validate exact agreement with prior fixed fixtures; do not change those files.
    prior = ROOT / 'inputs' / 'FIXTURE_RESULTS.json'
    prior_hash = hashlib.sha256(prior.read_bytes()).hexdigest()
    prior_results = {row['name']: row for row in json.loads(prior.read_text())['results']}
    for n in range(3, 8):
        assert tuple(prior_results[f'C{n}']['Z_coefficients']) == z[n]
        assert tuple(prior_results[f'C{n}']['H_coefficients']) == h[n]

    for n in range(RANGE_COEFFICIENTS[0], RANGE_COEFFICIENTS[1] + 1):
        assert all(c >= 0 for c in h[n])
        assert (h[n] == (0,)) == (n == 4)
        if n >= 5:
            assert h[n][:3] == (0, 0, 0) and h[n][3] == n
        if n >= 7:
            forcing = add(mul((0, 0, 0, 1, 3), z[n-5]),
                          mul((0, 0, 0, 0, 2, 3), z[n-6]))
            rhs = add(h[n-1], shift(h[n-2], 1), shift(h[n-3], 2), forcing)
            assert h[n] == rhs
        if n >= 4:
            raw_forcing = add(mul((-1, -3), z[n-1]), mul((1, 2), z[n-2]),
                              mul((0, 2, 5, 3), z[n-3]))
            assert h[n] == add(h[n-1], shift(h[n-2], 1), shift(h[n-3], 2), raw_forcing)

    report = {
        'status': 'all exact assertions passed',
        'proof_not_finite_evidence': 'The all-order theorem is proved by the identity and induction in the article.',
        'enumeration_orders': list(range(3, 13)),
        'coefficient_orders': list(range(3, 41)),
        'forcing_recurrence_orders': list(range(7, 41)),
        'subset_masks_per_predicate': sum(1 << n for n in range(3, 13)),
        'formal_initial_traces': {str(n): list(z[n]) for n in (1, 2, 3)},
        'bases': {str(n): list(h[n]) for n in bases},
        'symbolic_forcing_reduction': symbolic_forcing_check(),
        'enumerations': enumerations,
        'all_exact_polynomials': {str(n): {'Z': list(z[n]), 'H': list(h[n])}
                                  for n in range(3, 41)},
        'prior_fixture_results_sha256_read_only': prior_hash,
    }
    (ROOT / 'EXACT_CHECKS.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS: formal traces Z1,Z2,Z3 and polynomial Cayley-Hamilton matrix identity')
    print('PASS: exact subset/domination enumeration C3-C12; 8184 masks in each predicate')
    print('PASS: setwise complementation and coefficient reversal for C3-C12')
    print('PASS: formal-symbol forcing reduction (three recurrence substitutions)')
    print('PASS: C3-C7 agree exactly with read-only prior fixture coefficient records')
    print('PASS: direct H, coefficient nonnegativity/equality and z^3 coefficient for n=3..40')
    print('PASS: raw H recurrence n=4..40 and positive forcing recurrence n=7..40')
    print('Finite checks validate implementation; the article gives the proof for all orders.')


if __name__ == '__main__':
    main()
