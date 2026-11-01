#!/usr/bin/env python3
"""Run all proof checks on a fresh copy, with numerical optimizers disabled."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Combined proof replay refuses optimized Python')
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check_manifest(base,name):
 entries=[]
 for line in (base/name).read_text().splitlines():
  wanted,relative=line.split(maxsplit=1);p=Path(relative)
  if p.is_absolute()or '..'in p.parts:raise RuntimeError('Unsafe manifest path')
  if digest(base/p)!=wanted:raise RuntimeError('Manifest mismatch: '+relative)
  entries.append(relative)
 return entries

def main(destination=None):
 started=time.monotonic();entries=check_manifest(ROOT,'PROOF_MANIFEST.sha256')
 if destination:
  dest=Path(destination).resolve()
  if dest.exists():raise RuntimeError('Refusing existing destination')
  dest.mkdir(parents=True)
 else:dest=Path(tempfile.mkdtemp(prefix='weighted_domination_replay_'))
 work=dest/'proof';work.mkdir();logs=dest/'logs';logs.mkdir()
 for rel in entries+['PROOF_MANIFEST.sha256']:
  p=work/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,p)
 guard=work/'interval_certificate/inputs/upstream/tests/no_optimizer'
 env=os.environ.copy();env.update(PYTHONPATH=str(guard),PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
 tests=[]
 def add(name,folder,args,success=True):tests.append((name,folder,args,success))
 def command(body):return ['-c',body]
 add('optimizer_guard',work,['-c','import scipy.optimize'],False)
 add('subcubic_complete_proof',work,['replay_subcubic_proof.py','--destination',str(dest/'subcubic')])
 D=work/'degree_four';S=D/'structural';P=D/'point';I=D/'interval';E=D/'equality'
 add('reduced_grammar_and_exposures',S,['enumerate_validate.py'])
 add('independent_reduced_grammar',S,['audit/grammar_count_independent.py'])
 add('count_transfer_controls',S,['audit/transfer_fixture_checks.py'])
 add('point_primary_full_replay',P,['engine.py','replay'])
 for cmd in ['replay','compare','validate','witness','controls']:
  add('point_independent_'+cmd,P,['audit/independent_check.py',cmd])
 add('point_supplemental_controls',P,['audit/supplemental_checks.py'])
 add('weighted_primary_symbolic_regeneration',I,['symbolic_engine.py','enumerate'])
 add('weighted_independent_symbolic_regeneration',I,['audit/independent_weighted.py','enumerate'])
 add('weighted_independent_fixture_identities',I,['audit/independent_weighted.py','fixtures'])
 add('weighted_complete_feature_comparison',I,['audit/independent_weighted.py','compare'])
 add('weighted_cleared_factor_digest',I,['clearing_digest.py'])
 add('weighted_independent_gate_controls',I,['audit/gate_controls.py'])
 add('weighted_exact_unit_tests',I,['solver_unit_tests.py'])
 add('degree_four_selected_encoding',D,['check_encoding.py'])
 add('weighted_primary_full_certificate',I,['verify_primary.py','certificate.json'])
 add('weighted_independent_full_certificate',I,['audit/independent_weighted.py','verify','--certificate','certificate.json','--output','audit/certificate_verification.json'])
 add('weighted_primary_corruption_controls',I,['primary_controls.py','certificate.json'])
 add('weighted_independent_corruption_controls',I,['audit/certificate_controls.py','certificate.json'])
 add('equality_bind_fresh_replay_evidence',D,command("from pathlib import Path;import shutil;I=Path('interval');E=Path('equality/inputs');[(shutil.copy2(I/a,E/b))for a,b in [('certificate.json','certificate.json'),('primary_certificate_replay.json','primary_replay.json'),('audit/certificate_verification.json','independent_certificate_verification.json'),('symbolic_columns.json.gz','primary_columns.json.gz'),('audit/independent_columns.json.gz','independent_columns.json.gz')]]"))
 add('equality_exact_slacks_and_zero_views',E,['check_equality.py'])
 add('endpoint_empty_degree_four_reconstruction',E,['check_endpoint.py'])
 add('equality_structural_predicate',E,['check_structure.py'])
 add('endpoint_independent_direct_arithmetic',E,['audit/check_endpoint.py'])
 add('cycle_analytic_identity',E,['audit/check_analytic.py'])
 add('source_implication_exact_comparisons',D/'comparisons',['checks/check_implications.py'])
 add('combined_cycle_singleton_exposures',D,['check_cycle_exposures.py'])
 for name,folder,script in [('grammar',S,'enumerate_validate.py'),('point',P,'engine.py'),('weighted_primary',I,'verify_primary.py'),('weighted_independent',I,'audit/independent_weighted.py'),('equality',E,'check_equality.py'),('endpoint',E,'check_endpoint.py'),('structure',E,'check_structure.py'),('comparisons',D/'comparisons','checks/check_implications.py'),('cycles',D,'check_cycle_exposures.py'),('harness',work,'replay_proof.py')]:
  add('optimized_refusal_'+name,folder,['-O',script],False)
 results=[]
 for name,cwd,args,expected in tests:
  begin=time.monotonic();r=subprocess.run([sys.executable,*args],cwd=cwd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  out=r.stdout.replace(str(work).encode(),b'<proof-copy>').replace(str(dest).encode(),b'<replay-work>').replace(str(ROOT).encode(),b'<release>')
  passed=(r.returncode==0)==expected
  if name=='optimizer_guard':passed=passed and b'Optimizer import deliberately disabled'in out
  if name.startswith('optimized_refusal_'):passed=passed and (b'refuse'in out.lower()or b'optimized'in out.lower()or b'without -o'in out.lower())
  (logs/(name+'.log')).write_bytes(out)
  results.append(dict(name=name,expected_success=expected,exit_code=r.returncode,passed=passed,seconds=round(time.monotonic()-begin,3),log_sha256=hashlib.sha256(out).hexdigest()))
  print(name,'PASS'if passed else'FAIL',flush=True)
  if not passed:break
 report=dict(all_passed=len(results)==len(tests)and all(t['passed']for t in results),planned_tests=len(tests),proof_files_checked=len(entries),proof_manifest_sha256=digest(ROOT/'PROOF_MANIFEST.sha256'),history_copied_or_required=False,optimizer_imports_deliberately_disabled=True,numerical_optimization_or_recovery_executed=False,source_files_mutated=False,subcubic_certificate_sha256=digest(work/'interval_certificate/subcubic_certificate.json'),degree_four_certificate_sha256=digest(I/'certificate.json'),tests=results,elapsed_seconds=round(time.monotonic()-started,3))
 if (dest/'subcubic/REPLAY_REPORT.json').exists():report['subcubic_report']=json.loads((dest/'subcubic/REPLAY_REPORT.json').read_text())
 for label,path in [('weighted_primary',I/'primary_certificate_replay.json'),('weighted_independent',I/'audit/certificate_verification.json'),('equality',E/'EXACT_CHECKS.json'),('endpoint',E/'ENDPOINT_EMPTY_DEGREE4.json'),('structure',E/'STRUCTURAL_CHECKS.json')]:
  if path.exists():
   obj=json.loads(path.read_text());report[label]=obj.get('summary',obj)
 (dest/'REPLAY_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k not in ['tests','subcubic_report','weighted_primary','weighted_independent','equality','endpoint','structure']},indent=2))
 if not report['all_passed']:raise RuntimeError('Proof replay failed; all outcomes retained')
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--destination');a=p.parse_args();main(a.destination)
