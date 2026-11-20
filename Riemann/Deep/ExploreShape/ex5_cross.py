import sys, json; sys.path.insert(0, '.')
from fractions import Fraction as Fr
from validate import minform
from exact5 import check5
c = json.load(open('logs/cert5_m3.355.json'))
a = [Fr(x) for x in c['a']]; ph = {k: Fr(v) for k, v in c['phi'].items()}
phv = {(1, 2): float(ph['12']), (1, 4): float(ph['14']), (3, 2): float(ph['32']), (3, 4): float(ph['34'])}
print('certified form: validate N30 (rel. eig, g in [3e-3, pi]):', minform([float(x) for x in a], phv, 5, 30, 300, smin=1e-6))
for da in (Fr(-1, 1000), Fr(-1, 10000)):
    b = a[:]; b[0] += da
    r = check5(b, ph, ng=300, KM=3000)
    print('a0 %+g: Lemma 3 worst %.3g at g=%.3f (%s); validate N30: %s' % (float(da), r['worst'][0], r['worst'][1], r['worst'][2], minform([float(x) for x in b], phv, 5, 30, 300, smin=1e-3)))
