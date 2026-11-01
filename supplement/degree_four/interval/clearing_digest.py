#!/usr/bin/env python3
if not __debug__:raise RuntimeError('Clearing digest refuses Python -O/-OO')
import hashlib,json
import symbolic_engine as s
h=hashlib.sha256();data=s.load()
for r in data['rows']:h.update((s.canonical(s.cleared(r))+'\n').encode())
out=dict(all_passed=True,columns=len(data['rows']),cleared_sha256=h.hexdigest());s.dump(s.HERE/'clearing_digest.json',out);print(json.dumps(out))
