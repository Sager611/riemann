# B1 cross-check: finite zero-pad chain (Euler stages, wb = R line through (R, z)) with the twisted bookkeeping, vs the ODE of b1_twist.py;
# and the rho needed for theta = 0.873.
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from b1_twist import make, contact
import numpy as np
def chain(rho, STEP=1 / 40000, LEFT=0.725):
    c = rho / 2; k0 = 1 / (2 - c); R = (5 / 6 - 0.5 * k0) / (1.25 - k0); z = 0.0
    rs = [R]; zs = [0.0]
    while R - STEP > LEFT:
        k = (1.25 * R - 5 / 6) / (R - 0.5 - z)          # k = thr (gate 1 tight), line through (R, z) => wb = R
        r = R - STEP; H = r / R; u = (1 - H) / (6 - rho)
        z = H * z + u * (1 / k - 2 + rho / 2); R = r; rs.append(R); zs.append(z)
    rs = np.array(rs[::-1]); zs = np.array(zs[::-1])
    return lambda w: np.interp(w, rs, zs)
if __name__ == '__main__':
    for rho in (0.0, 24 / 289, 0.14516):
        Gc = chain(rho); Go, r0 = make(rho, 'correct')
        print("rho=%.5f chain z(0.7523)=%.6e ode z(0.7523)=%.6e  chain B*=%.12f ode B*=%.12f" % (rho, Gc(0.7523), Go(0.7523), contact(Gc)[0], contact(Go)[0])); sys.stdout.flush()
    lo, hi = 0.5, 1.5
    for _ in range(30):
        m = (lo + hi) / 2; B = contact(make(m, 'correct')[0])[0]
        if B > 0.873: lo = m
        else: hi = m
    print("rho with theta = 0.873 (consistent bookkeeping): %.5f" % ((lo + hi) / 2))
