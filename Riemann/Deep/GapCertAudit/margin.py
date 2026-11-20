"""gapcertaudit: margin of H for a certJ/certJ5 output: terms a_j kappa^(2j) c_j, their radii, rad(H)/H, and (optional second
file) the comparison of two truncations.  usage: margin.py certJ.json [certJ_other.json]"""
import sys, json
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 40


def ball(s):
    s = s.strip().strip('[]')
    if '+/-' in s:
        m, r = s.split('+/-'); return mp.mpf(m), mp.mpf(r)
    return mp.mpf(s), mp.mpf(0)


A = json.load(open(sys.argv[1]))
k2 = mp.mpf(Fr(A['kappa']).numerator) / Fr(A['kappa']).denominator; k2 = k2 * k2
a = [mp.mpf(Fr(x).numerator) / Fr(x).denominator for x in A['a']]
c = [ball(x) for x in A['c']]; H = ball(A['H']); s = [ball(x) for x in A['s']]
print('D', A['D'], 'prec', A['prec'], 'certified', A['certified'], 'H', A['H'], 'H/c0', A['H_over_c0'])
tot = mp.mpf(0); radsum = mp.mpf(0)
for j in range(len(a)):
    t = -a[j] * k2 ** j * c[j][0]; rt = abs(a[j]) * k2 ** j * c[j][1]
    tot += t; radsum += rt
    print('  j=%d  -a_j k^2j c_j = %s  (rad %s, rel %s)' % (j, mp.nstr(t, 12), mp.nstr(rt, 3), mp.nstr(c[j][1] / abs(c[j][0]), 3)))
print('  sum of terms %s = H mid %s; |H|/max|term| = %s' % (mp.nstr(tot, 10), mp.nstr(H[0], 10),
      mp.nstr(abs(H[0]) / max(abs(a[j] * k2 ** j * c[j][0]) for j in range(len(a))), 4)))
print('  rad(H) %s <= sum_j |a_j| k^2j rad(c_j) = %s; H/rad(H) = %s' % (mp.nstr(H[1], 3), mp.nstr(radsum, 3), mp.nstr(H[0] / H[1], 4)))
print('  s radii (rel):', [mp.nstr(x[1] / x[0], 2) for x in s])
print('  raw tails:', A['tails'])
if len(sys.argv) > 2:
    B = json.load(open(sys.argv[2])); H2 = ball(B['H'])
    print('second: D', B['D'], 'prec', B['prec'], 'certified', B['certified'], 'H', B['H'])
    print('  |H1 - H2| = %s, rad1 + rad2 = %s, overlap: %s' % (mp.nstr(abs(H[0] - H2[0]), 3), mp.nstr(H[1] + H2[1], 3),
          abs(H[0] - H2[0]) <= H[1] + H2[1]))
    for m in range(len(s)):
        s2 = ball(B['s'][m]); print('  s_%d: |diff|/s = %s, rad1/s = %s, rad2/s = %s' % (m, mp.nstr(abs(s[m][0] - s2[0]) / s[m][0], 2),
              mp.nstr(s[m][1] / s[m][0], 2), mp.nstr(s2[1] / s[m][0], 2)))
