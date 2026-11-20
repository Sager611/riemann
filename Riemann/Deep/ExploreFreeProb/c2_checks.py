"""C2 checks: (i) rho_3 of the unreversed J^{d,n} vs 1/(3 sqrt(n+d)); (ii) the cubic basin is an interval (spot checks)."""
import sys
from ffcum import *
from c2_basin import coeffs_from_logsym, is_real_rooted, r3
from fractions import Fraction as Fr
g = seq['xi_t0']
for d, n in ((4, 0), (4, 100), (4, 1000), (20, 100), (20, 1000), (60, 1000)):
    r, l = rho_orig(g, n, d, 4)
    rr, _ = rho_rev(g, n, 4)
    print('d=%3d n=%5d  rho3_orig = %s  *3sqrt(n+d) = %s   rho3_rev = %s' % (d, n, mp.nstr(r[3], 6), mp.nstr(3*r[3]*mp.sqrt(n+d), 6), mp.nstr(rr[3], 6)))
for d in (10, 40):
    lo, hi = r3(d, 30)
    pts = [lo*Fr(k, 8) for k in range(1, 8)] + [-lo*Fr(k, 8) for k in range(1, 8)]
    ok = all(is_real_rooted(coeffs_from_logsym({2: Fr(-1, 2), 3: e}, d), d) for e in pts)
    bad = [float(e) for e in (hi*Fr(9, 8), hi*Fr(3, 2), hi*2, hi*4) if is_real_rooted(coeffs_from_logsym({2: Fr(-1, 2), 3: e}, d), d)]
    print('d=%d: real-rooted at 14 points in (-r3, r3): %s; real-rooted points among 1.125,1.5,2,4 x r3: %s' % (d, ok, bad))
