#!/usr/bin/env python3
"""Portable exact replay; the analytic premises are proved in the manuscript.

This wrapper uses the unchanged primary enumeration/checking modules, without
using an optimizer.
It certifies finite arithmetic, not the mathematical coverage proof.
"""
if not __debug__:
    raise RuntimeError('Exact replay requires ordinary Python, not -O or -OO')
from pathlib import Path
import argparse, json, time
import local_views_v3 as primary
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('kind',choices=('primary','secondary'))
a=p.parse_args(); start=time.time()
cert=json.loads((HERE/f'certificate_{a.kind}_v3.json').read_text())
result=primary.evaluate_certificate(cert,primary.load(),regenerate=True)
result['elapsed_seconds']=time.time()-start
(HERE/f'replay_{a.kind}.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('equality_sums','inequality_sums')},indent=2))
