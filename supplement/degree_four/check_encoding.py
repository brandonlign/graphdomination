#!/usr/bin/env python3
"""Exact descriptive metrics for the selected degree-four witness; no search."""
if not __debug__:raise RuntimeError('Encoding metrics refuse optimized Python')
from pathlib import Path
from fractions import Fraction as Q
from math import lcm
import json,hashlib
R=Path(__file__).resolve().parent;p=R/'interval/certificate.json';c=json.loads(p.read_text())
if hashlib.sha256(p.read_bytes()).hexdigest()!='d110d283ea339b45bdf03879b6bff26fb80bfa1d3d5b73c6b579d67a386dae9b':raise RuntimeError('Selected certificate differs')
polys=[[Q(x)for x in row]for row in c['bernstein_polynomials']];v=[x for p in polys for x in p]
nums=sum(abs(x.numerator).bit_length()for x in v);dens=sum(x.denominator.bit_length()for x in v);L=lcm(*(x.denominator for x in v))
o=dict(all_passed=True,certificate_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),nonzero_polynomials=sum(any(p)for p in polys),nonzero_coefficients=sum(bool(x)for x in v),total_numerator_bits=nums,total_denominator_bits=dens,total_fraction_bits=nums+dens,largest_integer_bits=max(max(abs(x.numerator).bit_length(),x.denominator.bit_length())for x in v),common_denominator_bits=L.bit_length(),canonical_json_bytes=len(json.dumps(c,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()),maximum_denominator_decimal_digits=max(len(str(x.denominator))for x in v),coefficient_slots=len(v))
expected=dict(nonzero_polynomials=196,nonzero_coefficients=538,total_numerator_bits=26124,total_denominator_bits=27268,total_fraction_bits=53392,largest_integer_bits=140,common_denominator_bits=143,canonical_json_bytes=22875,coefficient_slots=672)
if any(o[k]!=x for k,x in expected.items()):raise RuntimeError('Encoding metrics differ')
(R/'CERTIFICATE_ENCODING.json').write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2))
