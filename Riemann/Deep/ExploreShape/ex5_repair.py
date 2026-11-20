import sys, json; sys.path.insert(0, '.')
from fractions import Fraction as Fr
from exact5 import check5
l = json.load(open('logs/scanJ5_r2_d64_k3.355.json'))[0]
a0 = [Fr(x).limit_denominator(10**12) for x in l['a5']]
ph0 = {k.strip('()').replace(', ', ''): Fr(v).limit_denominator(10**12) for k, v in l['phi5'].items()}
def show(d0, e12, e32, d5):
    a = a0[:]; a[0] += d0; a[5] += d5
    ph = dict(ph0); ph['12'] += e12; ph['32'] += e32; ph['14'] -= e32
    r = check5(a, ph, ng=300, KM=4000)
    L = [a[2] + ph['12'], (ph['32'] + ph['14']) / 2, a[4] + ph['34']]
    print(float(d0), float(e12), float(e32), float(d5), 'pmin %.4g' % r['pmin'], 'lead', {k: '%.3g' % v for k, v in r['lead'].items()},
          'Ceig %.3g' % r['C_eigs'][0], 'Linf det %.3g' % float(L[0] * L[2] - L[1] ** 2), 'worst %.3g @%.3f %s' % r['worst'],
          'DETroots', ['%.4g' % x for x in r['roots']['DET']], 'N11roots', ['%.4g' % x for x in r['roots']['N11']], flush=True)
for args in [(Fr(1,10**4), Fr(2,10**6), Fr(2,10**7), 0), (Fr(1,10**4), Fr(1,10**6), Fr(2,10**7), 0), (Fr(2,10**4), Fr(5,10**6), Fr(5,10**7), 0)]:
    show(*args)
