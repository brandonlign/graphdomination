#!/usr/bin/env python3
"""Exact sparse-polynomial check of the rational-bound certificate.

No symbolic-algebra package and no import of primary implementation. The
numerator is independently expanded from its five cleared rational summands.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

if not __debug__:
    raise RuntimeError('Optimized Python is not accepted.')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent / 'middle_interval'


def add(*polynomials):
    result = {}
    for poly in polynomials:
        for key, value in poly.items():
            result[key] = result.get(key, 0) + value
    return {key: value for key, value in result.items() if value}


def scale(poly, value):
    return {key: value * coeff for key, coeff in poly.items() if value * coeff}


def mul(*polynomials):
    result = {(0, 0): 1}
    for poly in polynomials:
        product = {}
        for (i, j), value in result.items():
            for (k, ell), coeff in poly.items():
                key = (i + k, j + ell)
                product[key] = product.get(key, 0) + value * coeff
        result = {key: value for key, value in product.items() if value}
    return result


def power(poly, exponent):
    return mul(*([poly] * exponent))


def compare(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    # Initial bivariate variables are q,x.
    one = {(0, 0): 1}
    q = {(1, 0): 1}
    x = {(0, 1): 1}
    r = add(one, x)
    a = add(scale(one, 4), scale(q, -3))
    xp2 = add(x, scale(one, 2))
    rmq = add(r, scale(q, -1))
    omq = add(one, scale(q, -1))
    omq3 = add(one, scale(power(q, 3), -1))
    # D=48*r^2*(x+2)^3*(r-q)*(1-q)*(1-q^3)>0.
    # N=D*(rational lower bound)/q.
    numerator = add(
        scale(mul(a, add(scale(one, 2), mul(power(q, 3), a)),
                  add(x, scale(power(q, 2), 2)), power(r, 2), power(xp2, 3),
                  rmq, omq, omq3), 3),
        scale(mul(power(x, 2), add(scale(power(x, 2), 3), scale(x, 10), scale(one, 12)),
                  r, rmq, omq, omq3), -16),
        scale(mul(q, power(x, 2), add(scale(r, 3), scale(mul(r, q), -1), scale(q, -2)),
                  power(xp2, 3), omq3), -8),
        scale(mul(power(q, 3), add(scale(one, 2), scale(power(q, 3), -1)),
                  r, power(xp2, 3), rmq, omq), -24),
        scale(mul(x, add(scale(one, 6), scale(q, -3), scale(power(q, 2), -1)),
                  power(xp2, 3), rmq, omq3), -4),
    )
    # Substitute x=3*q^2+s, yielding powers q,s.
    shifted = {}
    for (q_degree, x_degree), coeff in numerator.items():
        for s_degree in range(x_degree + 1):
            key = (q_degree + 2 * (x_degree - s_degree), s_degree)
            value = coeff * math.comb(x_degree, s_degree) * 3 ** (x_degree - s_degree)
            shifted[key] = shifted.get(key, 0) + value
    shifted = {key: value for key, value in shifted.items() if value}
    clearing_power = max(i for i, _ in shifted)
    # q=(1+2u)/(4*(1+u)); clear by [4*(1+u)]^clearing_power.
    transformed = {}
    for (q_degree, s_degree), coeff in shifted.items():
        for first_degree in range(q_degree + 1):
            first_coeff = math.comb(q_degree, first_degree) * 2 ** first_degree
            for second_degree in range(clearing_power - q_degree + 1):
                second_coeff = math.comb(clearing_power - q_degree, second_degree)
                key = (first_degree + second_degree, s_degree)
                value = coeff * 4 ** (clearing_power - q_degree) * first_coeff * second_coeff
                transformed[key] = transformed.get(key, 0) + value
    transformed = {key: value for key, value in transformed.items() if value}
    endpoint = {}
    for (q_degree, x_degree), coeff in numerator.items():
        for s_degree in range(x_degree + 1):
            value = (Fraction(coeff, 2 ** q_degree) * math.comb(x_degree, s_degree)
                     * Fraction(3, 4) ** (x_degree - s_degree))
            endpoint[s_degree] = endpoint.get(s_degree, Fraction(0)) + value
    endpoint = {key: value for key, value in endpoint.items() if value}
    primary = json.loads((ROOT / 'CERTIFICATE.json').read_text())
    primary_transformed = {(item['u_degree'], item['s_degree']): Fraction(item['coefficient'])
                           for item in primary['terms']}
    primary_endpoint = {item['s_degree']: Fraction(item['coefficient'])
                        for item in primary['endpoint_terms']}
    compare(primary_transformed == transformed, 'Transformed coefficient mismatch')
    compare(primary_endpoint == endpoint, 'Endpoint coefficient mismatch')
    compare(primary['clearing_power'] == clearing_power, 'Clearing-power mismatch')
    coefficient_pass = (all(value >= 0 for value in transformed.values())
                        and transformed.get((0, 0), 0) > 0)
    endpoint_pass = all(value >= 0 for value in endpoint.values()) and endpoint.get(0, 0) > 0
    compare(primary['coefficient_certificate_passed'] == coefficient_pass, 'Pass-status mismatch')
    compare(primary['endpoint_passed'] == endpoint_pass, 'Endpoint-status mismatch')
    output = {
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'method': 'Independent standard-library sparse integer polynomial reconstruction',
        'numerator_degrees_q_x': [max(i for i, _ in numerator), max(j for _, j in numerator)],
        'clearing_power': clearing_power,
        'nonzero_terms': len(transformed),
        'negative_terms': [dict(u_degree=i, s_degree=j, coefficient=str(value))
                           for (i, j), value in sorted(transformed.items()) if value < 0],
        'strict_constant': str(transformed.get((0, 0), 0)),
        'all_transformed_coefficients_match': True,
        'all_endpoint_coefficients_match': True,
        'coefficient_certificate_passed': coefficient_pass,
        'endpoint_passed': endpoint_pass,
        'full_compact_interval_passed': coefficient_pass and endpoint_pass,
        'endpoint_terms': [dict(s_degree=i, coefficient=str(value))
                           for i, value in sorted(endpoint.items())],
        'status': 'PASS' if coefficient_pass and endpoint_pass else 'MATCHED_BUT_CERTIFICATE_NOT_POSITIVE',
        'scope': 'The rational lower bound on the middle interval.'
    }
    (HERE / 'rational_independent_audit.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
