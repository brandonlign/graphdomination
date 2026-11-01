#!/usr/bin/env python3
"""Independent optimizer-free rejection controls for an already accepted candidate."""
from __future__ import annotations
if not __debug__:raise RuntimeError('Independent certificate controls refuse optimized Python')
from pathlib import Path
from fractions import Fraction as F
import argparse,contextlib,copy,gzip,io,json,subprocess,sys,time
import independent_weighted as a
HERE=Path(__file__).resolve().parent

def main(path):
    start=time.monotonic();base=json.loads(Path(path).read_text());dest=HERE/'corruptions';dest.mkdir(exist_ok=True);out=[]
    def reject(name,certificate,bundle=None,full=False):
        cp=dest/(name+'.json');a.dump(cp,certificate);oldhere=a.HERE
        if bundle is not None:
            directory=dest/name;directory.mkdir(exist_ok=True)
            with gzip.open(directory/'independent_columns.json.gz','wt') as f:json.dump(bundle,f,sort_keys=True)
            a.HERE=directory
        captured=io.StringIO()
        try:
            try:
                with contextlib.redirect_stdout(captured):a.verify_certificate(cp,full=full)
            except (RuntimeError,ValueError,KeyError,TypeError,IndexError) as e:
                (dest/(name+'.log')).write_text(captured.getvalue()+type(e).__name__+': '+str(e)+'\n')
                out.append(dict(name=name,rejected=True,exception=type(e).__name__,message=str(e)));return
            raise RuntimeError('Corruption was accepted: '+name)
        finally:a.HERE=oldhere
    for field,value in [('schema_version','unsupported'),('interval',['0','2']),('clearing_factor','6*A*Z_C')]:
        bad=copy.deepcopy(base);bad[field]=value;reject('wrong_'+field,bad)
    bad=copy.deepcopy(base);bad['generation_sha256']='0'*64;reject('wrong_scope_digest',bad)
    bad=copy.deepcopy(base);bad['support']=bad['support'][:-1];reject('missing_transport_support',bad)
    bad=copy.deepcopy(base);bad['bernstein_polynomials'][170][0]='-1';reject('negative_inequality_multiplier',bad)
    bad=copy.deepcopy(base);bad['bernstein_polynomials']=[['0']*len(p) for p in bad['bernstein_polynomials']];reject('all_zero_candidate',bad)
    bad=copy.deepcopy(base);j=a.SUPPORT.index('BR_X_2_2_2');bad['bernstein_polynomials'][j]=[str(F(x)+1) for x in bad['bernstein_polynomials'][j]];reject('boundary_reversal_coefficient_changed',bad)
    with gzip.open(HERE/'independent_columns.json.gz','rt') as f:bundle=json.load(f)
    bad=copy.deepcopy(bundle);bad['rows'][0]['Z'][0]+=1;reject('local_partition_corruption',base,bad,True)
    bad=copy.deepcopy(bundle);bad['full_to_column'][0]=(bad['full_to_column'][0]+1)%len(bad['rows']);reject('full_view_mapping_corruption',base,bad,True)
    bad=copy.deepcopy(bundle);bad['rows'][0]['full_views']+=1;reject('multiplicity_corruption',base,bad,True)
    r=subprocess.run([sys.executable,'-O',str(HERE/'certificate_controls.py')],capture_output=True,text=True)
    (dest/'optimized_rejection.log').write_text(r.stdout+r.stderr);a.require(r.returncode!=0 and 'refuse' in r.stderr.lower(),'-O controls refusal')
    out.append(dict(name='optimized_certificate_controls',rejected=True,returncode=r.returncode))
    result=dict(all_passed=True,controls=out,elapsed_seconds=time.monotonic()-start,optimizer_free=True)
    a.dump(HERE/'certificate_controls.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('certificate');args=p.parse_args();main(args.certificate)
