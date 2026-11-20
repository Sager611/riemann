"""C1: the brief's s(d,n): largest backward Hermite time keeping J*^{d,n} real-rooted, in units of its variance,
sigma* = sup{sigma in [0,1]: exp(-(1-sigma)D^2/2 + sum_k rho_k D^k) x^d real-rooted}  (rho_k from ffcum, 30-digit rationals).
Compared for xi (t = 0 and de Bruijn-Newman time t = 0.2) and DH: a Gamma-factor (Taylor) quantity, blind to zeros and to t."""
from ffcum import *
from c2_basin import coeffs_from_logsym, is_real_rooted
from fractions import Fraction as Fr
def sigma_star(rho, d, steps=24):
    lo, hi = Fr(0), Fr(1)
    for _ in range(steps):
        mid = (lo+hi)/2
        l = {2: -(1-mid)/2}
        for k in range(3, d+1): l[k] = Fr(mp.nstr(rho[k], 30))
        if is_real_rooted(coeffs_from_logsym(l, d), d): lo = mid
        else: hi = mid
    return float(lo)
print(' d     n   sigma*(xi,t=0)  sigma*(xi,t=0.2)  sigma*(DH)')
for d in (3, 4, 5, 8):
    for n in (0, 10, 100, 1000):
        row = []
        for name in ('xi_t0', 'xi_t0.2', 'dh_t0'):
            rho, _ = rho_rev(seq[name], n, d)
            row.append('%.6f' % sigma_star(rho, d))
        print('%2d %5d   %s' % (d, n, '        '.join(row)), flush=True)
