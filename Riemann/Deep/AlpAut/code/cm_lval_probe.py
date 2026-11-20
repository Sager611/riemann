#!/usr/bin/env python3
# alpaut probe P4 (CM example): odd-side L-values L(0, chi_{-D}) = 2h(-D)/w(-D) of imaginary quadratic (CM) fields can be divisible by
# a large power of l, but only up to log_l of their archimedean size (h < sqrt(D) log D): polynomial in log D. Scan fundamental -D < DMAX.
import sys, math, cypari2
pari = cypari2.Pari(); pari.allocatemem(64 * 10**6)
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
best = {3: (0, 0, 0), 5: (0, 0, 0), 7: (0, 0, 0)}
ratio = {3: 0.0, 5: 0.0, 7: 0.0}
n = 0
for D in range(3, DMAX):
    if not int(pari('isfundamental(%d)' % (-D))): continue
    h = int(pari('qfbclassno(%d)' % (-D))); n += 1
    for l in best:
        v = 0; x = h
        while x % l == 0: x //= l; v += 1
        if v > best[l][0]: best[l] = (v, D, h)
        bound = math.log(math.sqrt(D) * math.log(D) + 1) / math.log(l)
        ratio[l] = max(ratio[l], v / bound)
print('fundamental -D, D <', DMAX, ':', n, 'fields')
for l in best:
    v, D, h = best[l]
    print('l=%d: max v_l(h(-D)) = %d at D=%d (h=%d); max v_l(h)/log_l(sqrt(D) log D) = %.3f' % (l, v, D, h, ratio[l]))
