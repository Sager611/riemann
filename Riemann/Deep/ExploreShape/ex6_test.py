import sys, json; sys.path.insert(0, '.')
from fractions import Fraction as Fr
from exact5 import numer5
from exact6 import numer6
# consistency: J6 formulas with a6 = phi52 = phi54 = 0 must reproduce Lemma 3 numerators exactly
c = json.load(open('logs/cert5_m80_3.358.json'))
a = [Fr(x) for x in c['a']]; ph = {k: Fr(v) for k, v in c['phi'].items()}
r5 = numer5(a, ph); r6 = numer6(a + [Fr(0)], ph)
print('J6->J5 reduction identical:', all(list(x) == list(y) for x, y in zip(r5[:5], r6[:5])) and r5[5] == r6[5])
