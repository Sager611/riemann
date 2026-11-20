"""explore53: summarise a certJ output (ball midpoints -> standardized central moments; relative radii)."""
import sys, json
from math import comb
d = json.load(open(sys.argv[1]))
def mid_rad(s):
    s = s.strip().lstrip('[').rstrip(']')
    if '+/-' in s:
        a, b = s.split('+/-'); return float(a), float(b)
    return float(s), 0.0
s = [mid_rad(x)[0] for x in d['s']]; rr = [mid_rad(x)[1] for x in d['s']]
m = [x / s[0] for x in s]; mu = m[1]; K = len(m) - 1
cen = [sum(comb(k, i) * m[i] * (-mu) ** (k - i) for i in range(k + 1)) for k in range(K + 1)]
sd = cen[2] ** 0.5
print('mean %.7f var %.7f' % (mu, cen[2]), 'standardized', ' '.join('%.5f' % (cen[k] / sd ** k) for k in range(K + 1)))
print('relative radii of s_m:', ' '.join('%.1e' % (rr[k] / abs(s[k])) for k in range(K + 1)))
print('H =', d['H'], ' H/c0 =', d['H_over_c0'], ' certified:', d['certified'], ' imag:', d['imag_check'], ' times', d['times_s'], ' peak', d['peak_MiB'])
