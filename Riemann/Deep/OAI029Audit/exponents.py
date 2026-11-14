# Exact rational recomputation of every exponent in 029's zero-exclusion core (02a l.69-76, 06-zero-free l.165-683).
from fractions import Fraction as Fr
a, b = Fr(9, 10), Fr(133, 1000); M0 = a + b; h0 = 1 - a; d1 = Fr(1, 10**7); s1 = 1 + d1; z1 = Fr(1, 6) + d1
E = lambda sx, sw, sz: a/2 + sx - 1 + h0*sz + b*(sw - 1)      # (8.x) contour exponent
C1 = a/2 + 1 - 1 + h0/6
out = {}
out['M0'] = M0; out['C*(1)'] = C1; out['C*(1)==7/15'] = (C1 == Fr(7, 15))
out['reflection (2M0-1-b)/2'] = (2*M0 - 1 - b)/2; out['margin 7/15-933/2000'] = C1 - (2*M0 - 1 - b)/2
out['margin == (4-4a-3b)/6'] = (C1 - (2*M0 - 1 - b)/2 == (4 - 4*a - 3*b)/6)
out['R_w'] = E(s1, Fr(95, 100), z1) - C1; out['R_z'] = E(s1, 1, Fr(16, 100)) - C1
base = E(s1, Fr(3, 1000), Fr(4, 10)) - C1; out['nonprincipal contour'] = base
out['U<=Z^.09 (alpha=1)'] = base + Fr(9, 100)*Fr(11, 10)
out['exceptional (alpha=1/2)'] = base + Fr(10001, 100000)*Fr(6, 10)
out['good central (xi=.998)'] = E(Fr(998, 1000), Fr(3, 1000), Fr(4, 10)) - C1 + Fr(10001, 100000)*Fr(11, 10)
out['good central == -377/300000'] = (out['good central (xi=.998)'] == Fr(-377, 300000))
out['good central without xi-move'] = base + Fr(10001, 100000)*Fr(11, 10)
R0 = 350000; far = E(4, 4, R0) - h0*R0; out['far: E*(4,4,R0)-.1R0'] = far
out['far at U=Z^.10001'] = far + h0*R0 + Fr(10001, 100000)*(1 - R0) - C1
out['detector 10(1/2-.996)+3/2'] = 10*(Fr(1, 2) - Fr(996, 1000)) + Fr(3, 2)
out['density 11*(2-2*.996)'] = 11*(2 - 2*Fr(996, 1000))
out['final saving 1/6000-1e-5 > 1/20000'] = (Fr(1, 6000) - Fr(1, 10**5) > Fr(1, 20000))
out['nonprincipal slack after 1e-4 loss vs .0004'] = -out['good central (xi=.998)'] - Fr(1, 10**4) - Fr(4, 10**4)
for k, v in out.items():
    print('%-45s %s' % (k, v if isinstance(v, bool) else '%s = %.9f' % (v, float(v))))
