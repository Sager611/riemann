import sys, json; sys.path.insert(0, '.')
from fractions import Fraction as Fr
from exact5 import check5
l = json.load(open('logs/scanJ5_r2_d64_k3.355.json'))[0]
a = [Fr(x).limit_denominator(10**12) for x in l['a5']]
ph = {k.strip('()').replace(', ', ''): Fr(v).limit_denominator(10**12) for k, v in l['phi5'].items()}
print('raw SDP form:', check5(a, ph, ng=200, KM=3000))
