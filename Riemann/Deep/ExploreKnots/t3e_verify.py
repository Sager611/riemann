# t3e_verify.py — independent mp.dirichlet check of the first Poincare-sphere H_- zero (from t3d_poincare.py)
import numpy as np, mpmath as mp
exec(open('t3d_poincare.py').read().split("fam = {")[0])
mp.mp.dps = 40
co = [mp.mpf(float(H[n % N])) for n in range(N)]
g = lambda s: mp.dirichlet(s, co)
z = mp.findroot(g, mp.mpc(0.604987667, 19.74619692), tol=1e-30)
print('Poincare H_minus zero', mp.nstr(z, 15), '|value|', mp.nstr(abs(g(z)), 3))
