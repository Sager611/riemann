# dumpplan.py -- run cert_ch.py's (untrusted, float) planner for the two record thetas and dump the plan; no certify, no writes in the repo.
import sys, os, json, math, time
HP = '/Users/adriansager/playground/riemann/Riemann/Deep/HP141/MertensCH'
sys.path.insert(0, HP)
import cert_ch as CC
from flint import arb
t0 = time.time()
pairs = [(cf['name'], CC.load(cf)) for cf in CC.CFG]
import bprime
assert bprime.BPs == "0.87334585610671858", bprime.__file__
pairs += [(cf['name'] + '_Bp', bprime.make(bprime.translate(cf))) for cf in CC.CFG]
certs = []
for nm, c in pairs:
    lC = c.tail(CC.TP)
    certs.append(dict(name=nm, c=c, a=c.a, s=c.s1, lC=lC, CS=CC.CLstar(c.s1), af=float(c.a.mid()), sf=float(c.s1.mid()),
                      lCf=float(lC.upper()), lCSf=float(CC.CLstar(c.s1).log().upper()), thstar=float(((c.s1 + c.a)/(1 + c.a)).upper())))
JL0 = arb(json.load(open(os.path.join(HP, 'logs', 'JL0.json')))['JL0_upper'])
RHP = JL0 + 2*CC.CL*CC.TP**CC.AL/CC.AL
consts = dict(RHPf=float(RHP.upper()), lTP=math.log(int(CC.TPS)))
out = dict(names=[c['name'] for c in certs], sigma1=[c['c'].d['sigma1'] for c in certs], JL0=str(JL0), plans={})
for th in [0.9409, 0.9382]:
    p = CC.make_plan(th, certs, consts)
    out['plans'][str(th)] = dict(L0=p['L0'], L1=p['L1'], log10c=p['log10c'], ranges=p['ranges'], fin=list(p['fin']))
    used = sorted(set(r[2] for r in p['ranges']))
    print(th, 'plan log10c', round(p['log10c'], 4), 'L0', p['L0'], 'L1', p['L1'], 'ranges', len(p['ranges']), 'fin', p['fin'][0],
          'used', [(certs[i]['name'], sum(1 for r in p['ranges'] if r[2] == i)) for i in used], round(time.time() - t0, 1), 's', flush=True)
json.dump(out, open('/Users/adriansager/playground/riemann/Riemann/Deep/MertensAudit/logs/plans.json', 'w'))
