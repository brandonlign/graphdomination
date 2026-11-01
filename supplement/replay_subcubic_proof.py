#!/usr/bin/env python3
"""Run the subcubic proof checks on a fresh copy of the files listed in the manifest."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Proof replay refuses optimized Python')
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parent
CERT_SHA='c6f20f9d56dd89013c8e1b186e1f0614368b7585c3243fa12ddc8b64fbfa0205'
def check_manifest(base,name):
 entries=[]
 for ln in (base/name).read_text().splitlines():
  want,rel=ln.split(maxsplit=1);p=base/rel
  if hashlib.sha256(p.read_bytes()).hexdigest()!=want:raise RuntimeError('Manifest mismatch: '+rel)
  entries.append(rel)
 return entries
def main(destination=None):
 start=time.time();entries=check_manifest(ROOT,'SUBCUBIC_PROOF_MANIFEST.sha256')
 dest=Path(destination)if destination else Path(tempfile.mkdtemp(prefix='weighted_audit_proof_'))
 if destination:
  if dest.exists():raise RuntimeError('Refusing existing replay destination')
  dest.mkdir(parents=True)
 W=dest/'proof';W.mkdir()
 for rel in entries+['SUBCUBIC_PROOF_MANIFEST.sha256']:
  p=W/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,p)
 C=W/'interval_certificate';U=C/'inputs/upstream';A=C/'audit/final_average_1000';A.mkdir(parents=True,exist_ok=True)
 cert=C/'subcubic_certificate.json'
 if hashlib.sha256(cert.read_bytes()).hexdigest()!=CERT_SHA:raise RuntimeError('Wrong certificate')
 # Invalid certificates for the negative controls are generated from the certificate.
 original=json.loads(cert.read_text());controls=C/'tests';controls.mkdir(exist_ok=True)
 for label in ['all_zero','negative_multiplier']:
  bad=json.loads(json.dumps(original))
  if label=='all_zero':bad['bernstein_polynomials']=[['0','0','0']for _ in bad['support']]
  else:
   j=next(i for i,n in enumerate(bad['support'])if n.startswith(('X_','B_','L_')));bad['bernstein_polynomials'][j][0]='-1'
  (controls/(label+'.json')).write_text(json.dumps(bad,indent=2)+'\n')
 env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONPATH=str(U/'tests/no_optimizer'))
 tests=[]
 def add(name,cwd,args,ok=True):tests.append((name,cwd,args,ok))
 def modcall(module,folder,body):return ['-c',f"import sys;sys.path.insert(0,{str(folder)!r});import {module} as m;{body}"]
 add('optimizer_import_guard',W,['-c','import scipy.optimize'],False)
 add('selected_support_dictionary',W/'examples',['check_selected_dictionary.py',str(W)])
 add('primary_symbolic_regeneration',U,['symbolic_engine.py'])
 add('independent_symbolic_regeneration',U,['audit/symbolic_features_independent.py'])
 add('complete_key_comparison',U,['audit/compare_complete_keys.py'])
 add('complete_cleared_polynomial_comparison',U,['audit/compare_clearing.py'])
 add('integer_bernstein_unit_tests',U,['audit/test_integer_bernstein.py'])
 add('selected_primary_full_certificate',U,['verify_primary.py',str(cert)])
 add('selected_independent_certificate',U,['audit/verify_interval_independent.py',str(cert),'--out',str(A/'independent_replay.json')])
 add('selected_independent_zero_support',C,['independent_zero_support.py',str(cert),str(A/'independent_replay.json'),str(A/'zero_support.json')])
 add('selected_independent_complete_zero_structure',C,['verify_low_bit_zero_structure.py'])
 add('direct_cycle_exposure_equality',W/'equality',['check_cycle_exposures.py'])
 add('worked_examples_primary',W/'examples',['check_examples.py'])
 add('worked_examples_independent',W/'examples/independent',['check_examples.py'])
 add('independent_view_grammar',W/'examples/independent',['count_view_grammar.py'])
 for kind in ['primary','secondary']:add('endpoint_'+kind,W/'endpoint',['replay_primary.py',kind])
 for name,args in [('validate',['validate_v3.py']),('enumerate',['enumerate_v3_independent.py']),('primary',['verify_v3.py','primary']),('secondary',['verify_v3.py','secondary'])]:add('endpoint_independent_'+name,W/'endpoint/independent',args)
 add('endpoint_third_checker',W/'endpoint',['challenge/check_exact.py'])
 add('cycle_primary',W/'cycles',['check_exact.py'])
 add('cycle_independent',W/'cycles',['audit/check_runs.py','--compare','EXACT_CHECKS.json'])
 for label in ['all_zero','negative_multiplier']:
  p=controls/(label+'.json')
  add(label+'_primary_rejected',U,['verify_primary.py',str(p)],False)
  add(label+'_independent_rejected',U,['audit/verify_interval_independent.py',str(p),'--out',str(A/(label+'.json'))],False)
 for name,folder,script in [('primary',U,'verify_primary.py'),('independent',U,'audit/verify_interval_independent.py'),('zero_support',C,'independent_zero_support.py'),('zero_structure',C,'verify_low_bit_zero_structure.py'),('cycle_exposure',W/'equality','check_cycle_exposures.py'),('endpoint',W/'endpoint','replay_primary.py'),('cycle',W/'cycles','check_exact.py'),('complete_keys',U,'audit/compare_complete_keys.py'),('harness',W,'replay_subcubic_proof.py')]:add('optimized_refusal_'+name,folder,['-O',script],False)
 logs=dest/'logs';logs.mkdir();results=[]
 for name,cwd,args,success in tests:
  t=time.time();p=subprocess.run([sys.executable,*args],cwd=cwd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  output=p.stdout.replace(str(W).encode(),b'<proof-copy>').replace(str(dest).encode(),b'<replay-work>')
  (logs/(name+'.log')).write_bytes(output);passed=(p.returncode==0)==success
  if name=='optimizer_import_guard':passed=passed and b'Optimizer import deliberately disabled' in output
  if name.startswith('optimized_refusal_'):passed=passed and (b'optimized' in output.lower()or b'-o' in output.lower())
  results.append(dict(name=name,expected_success=success,exit_code=p.returncode,passed=passed,elapsed_seconds=round(time.time()-t,3),log_sha256=hashlib.sha256(output).hexdigest()))
  print(name,'PASS'if passed else'FAIL',flush=True)
  if not passed:break
 report=dict(all_passed=len(results)==len(tests)and all(t['passed']for t in results),selected_certificate_sha256=CERT_SHA,proof_files_checked=len(entries),history_copied_or_required=False,old_interval_certificate_required=False,optimizer_imports_deliberately_disabled=True,source_files_mutated=False,planned_tests=len(tests),tests=results,elapsed_seconds=round(time.time()-start,3))
 (dest/'REPLAY_REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print('Report:',dest/'REPLAY_REPORT.json')
 if not report['all_passed']:raise RuntimeError('Proof replay failed; all outcomes retained')
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--destination');a=p.parse_args();main(a.destination)
