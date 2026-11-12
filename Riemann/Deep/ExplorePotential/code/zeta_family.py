# zeta_family.py (explore35, adversarial): move a phase-locked family of zeta zeros off the line (ordinates kept) and measure what the
# ordinate law and the on-line share see: Delta(all-ordinate law) = -2 E_off; Delta(on-line share) = -sum_F 2 x^{i gamma} cosh(delta log x).
import sys, math, numpy as np
g = np.load(sys.argv[1]); T = float(g.max()); N = len(g)
def Lam(n): return {2: math.log(2), 3: math.log(3), 4: math.log(2), 5: math.log(5), 7: math.log(7), 8: math.log(2), 9: math.log(3)}.get(n, 0.0)
fam = g[np.cos(g*math.log(2)) < -0.9]
print('N = %d, family {cos(gamma log 2) < -0.9}: %d zeros (%.1f%%)' % (N, len(fam), 100*len(fam)/N))
print('| delta | x | Landau P | -2 E_off (change of the all-ordinate law) | on-line share change | Landau error scale max abs(S - P) = 6.5 |')
for d in (0.01, 0.1, 0.3):
    for x in (2, 3, 6, 9):
        L = math.log(x); P = -(T/(2*math.pi))*Lam(x)/math.sqrt(x)
        E2 = -(2*np.exp(1j*fam*L)*(math.cosh(d*L) - 1)).sum().real
        sh = -(2*np.exp(1j*fam*L)*math.cosh(d*L)).sum().real
        print('| %.2f | %d | %+.1f | %+.2f | %s |  |' % (d, x, P, E2, ('%+.3f' % (sh/P)) if P else ('%+.1f (abs)' % sh)))
