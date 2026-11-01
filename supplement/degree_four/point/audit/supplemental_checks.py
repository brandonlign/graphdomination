#!/usr/bin/env python3
"""Additional optimizer-free graph-witness and checker-sensitivity tests."""
if not __debug__:raise RuntimeError('Supplemental exact audit refuses Python -O/-OO')
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import copy,json
import independent_check as checker
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
bundle=checker.load_independent()
fixture=json.loads((HERE/'fixture_checks.json').read_text())['fixtures']['C4']
expected=defaultdict(F)
for exposure in fixture['exposure_records']:
    expected[bundle['full_to_column'][exposure['full_view_id']]]+=F(exposure['mass'])
c4=json.loads((ROOT/'c4_primal.json').read_text());actual=defaultdict(F)
for item in c4['support']:actual[item['id']]+=F(item['weight'])
checker.require(dict(actual)==dict(expected),'C4 primal does not equal true graph-exposure pushforward')
cert=json.loads((ROOT/'certificate.json').read_text());perturbed=copy.deepcopy(cert)
perturbed['equality_dual'][1]=str(F(cert['equality_dual'][1])+1000000)
try: checker.verify_witness(perturbed,bundle)
except RuntimeError as error:
    checker.require(str(error).startswith('Negative dual slack'),'Transport corruption rejected for unexpected reason')
    corruption=dict(rejected=True,error=str(error))
else:raise RuntimeError('Bad transport multiplier accepted')
primal=json.loads((ROOT/'primal.json').read_text());primal_report=checker.verify_primal(primal,bundle)
report=dict(all_passed=True,c4_graph_support_exact=True,c4_support={str(k):str(v) for k,v in sorted(expected.items())},transport_slack_corruption=corruption,recovered_primal_objective=primal_report['objective'],recovered_primal_support_size=primal_report['support_size'],nonzero_dual_coordinates=sum(F(x)!=0 for x in cert['equality_dual']+cert['inequality_dual']))
checker.dump(HERE/'supplemental_checks.json',report);print(json.dumps(report,indent=2))
