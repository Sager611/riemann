import sys, json; sys.path.insert(0, '.')
from fractions import Fraction as Fr
from exact6 import check6
import os; l = json.load(open(os.environ.get('SRC', 'logs/scanJ6m_r2_d80_k3.365.json')))[0]
a = [Fr(x).limit_denominator(10**12) for x in l['a6']]
ph = {k.strip('()').replace(', ', ''): Fr(v).limit_denominator(10**12) for k, v in l['phi6'].items()}
r = check6(a, ph, ng=400, KM=3000)
print(r)
for d0 in (Fr(1, 10**4), Fr(3, 10**4)):
    b = a[:]; b[0] += d0
    print('a0 +', float(d0), check6(b, ph, ng=400, KM=3000)['worst'])
