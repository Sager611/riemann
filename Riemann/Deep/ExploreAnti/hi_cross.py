# explore55: bilinear form of Conjecture H at the 10^12 window: B_ij - u_i u_j = -P_i(1)P_j(1)/12 for P in {x, x^2, x^3, P_BHB, P_2D}.
import numpy as np, sys, json, re, mpmath as mp
from moll import mobius, bound_complex, bound_real
mp.mp.dps = 40; TWOPI = 2*np.pi; t0 = 267653395647; path = sys.argv[1]; th = float(sys.argv[2])
d = np.array([float(l) for l in open(path) if re.fullmatch(r'\s*[0-9]+\.[0-9]+\s*', l)]); t = t0 + d; T = t.max(); L = np.log(T/TWOPI)
zp, _ = np.load(path.replace('.txt', '_zp10000.npy'))
thm = np.array([float(mp.fmod(mp.siegeltheta(mp.mpf(t0)+mp.mpf(str(x))), 2*mp.pi)) for x in d]); D = -1j*np.exp(-1j*thm)*zp
y = T**th; mu = mobius(int(y)); n = np.nonzero(mu[1:int(y)+1])[0] + 1
r = np.array([float(mp.fmod(t0*mp.log(int(k)), 2*mp.pi)) for k in n]); ln = np.log(n.astype(float)); x = np.log(y/n)/np.log(y)
polys = {'x': x, 'x2': x**2, 'x3': x**3, 'BHB': -th*x**2 + (1+th)*x, 'P2D': ((0.5+th)*x - th*x**2/2)/(0.5+th/2)}
Wt = np.stack([mu[n]*n**-0.5*v for v in polys.values()], 1); M = np.empty((len(t), len(polys)), complex)
for a in range(0, len(t), 8):
    M[a:a+8] = np.exp(-1j*(r[None, :] + d[a:a+8, None]*ln[None, :])) @ Wt
Y = D[:, None]*M; ph = Y[:, 0].mean()/abs(Y[:, 0].mean()); Z = Y/ph; N = len(Y); names = list(polys)
U = Z.mean(0)/L; Bm = (Z.T @ Z)/N/L**2; out = dict(L=L, theta=th, y=y, N=N, Imu=[float(v.imag) for v in U])
out['B-uu'] = {f'{names[i]},{names[j]}': round(float((Bm[i, j] - U[i]*U[j]).real), 4) for i in range(5) for j in range(i, 5)}
out['CGG_BHB'] = bound_complex(Y[:, [3]]); out['2D_P2D'] = bound_real(Y[:, [4]]); out['CGG_J3'] = bound_complex(Y[:, :3]); out['2D_J3'] = bound_real(Y[:, :3])
out['pred_CGG'] = 1-(1+th)**-3; out['pred_2D'] = 1-(1+2*th)**-3
print(json.dumps(out), flush=True)
