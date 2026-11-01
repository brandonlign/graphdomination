#!/usr/bin/env python3
if not __debug__:raise RuntimeError('Primary controls refuse Python -O/-OO')
import copy,json,sys,subprocess
from pathlib import Path
import symbolic_engine as s
from solver_math import Problem
import verify_primary as v
path=Path(sys.argv[1]);cert=json.loads(path.read_text());data=s.load();p=Problem(cert['stage'],data);out=[]
def reject(name,fn):
 try:fn()
 except Exception as err:out.append(dict(control=name,rejected=True,type=type(err).__name__,message=str(err)));return
 raise RuntimeError('Corruption accepted: '+name)
zero=copy.deepcopy(cert);zero['bernstein_polynomials']=[['0']*(p.d+1)for _ in s.SUPPORT]
reject('zero_candidate',lambda:s.require(p.verify(zero)['all_passed'],'Zero candidate rejected'))
negative=copy.deepcopy(cert);negative['bernstein_polynomials'][170][0]='-1'
reject('negative_multiplier',lambda:s.require(p.verify(negative)['all_passed'],'Negative multiplier rejected'))
wrong=copy.deepcopy(cert);j=s.SUPPORT.index('BR_X_2_2_2');wrong['bernstein_polynomials'][j][0]=str(__import__('fractions').Fraction(wrong['bernstein_polynomials'][j][0])+1)
reject('changed_transport',lambda:s.require(p.verify(wrong)['all_passed'],'Changed transport rejected'))
wrong=copy.deepcopy(cert);wrong['support']=wrong['support'][:-1]
reject('missing_support',lambda:p.verify(wrong))
wrong=copy.deepcopy(data);wrong['rows'][wrong['full_to_column'][0]]['Z'][0]+=1
reject('corrupt_polynomial',lambda:v.full_views(wrong))
wrong=copy.deepcopy(data);wrong['full_to_column'][0]=(wrong['full_to_column'][0]+1)%len(wrong['rows'])
reject('corrupt_mapping',lambda:v.full_views(wrong))
wrong=copy.deepcopy(data);wrong['rows'][0]['full_views']+=1
reject('corrupt_multiplicity',lambda:v.full_views(wrong))
for script in('symbolic_engine.py','solver_math.py','verify_primary.py','primary_controls.py'):
 a=subprocess.run([sys.executable,'-O',str(s.HERE/script)],capture_output=True,text=True)
 s.require(a.returncode!=0 and 'refuse' in a.stderr.lower(),'-O control '+script)
 out.append(dict(control='optimized_'+script,rejected=True,returncode=a.returncode,stderr=a.stderr))
s.dump(s.HERE/'primary_controls.json',dict(all_passed=True,controls=out));print('All primary corruption and -O controls passed')
