# nb2fals e0: first look at per-N data (N <= NMAX): noise structure and spectrum in u = log N.
import sys, math, numpy as np, mpmath as mp
O = '/Users/adriansager/playground/riemann/.scratch/nb2fals/'
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
P = np.load(O + f'perN{NMAX}.npz')['perN']; m = P[:, 0] >= 100; P = P[m]
N, u, C2, kap, JmS, R = P[:, 0], P[:, 1], P[:, 2], P[:, 3], P[:, 4], P[:, 5]; Q = C2 - JmS
try: Z = np.load(O + 'zeros100.npz'); gam = Z['gam']; zp = Z['zp']
except Exception:
    mp.mp.dps = 20; gam = []; zp = []
    for k in range(1, 101):
        r = mp.zetazero(k); gam.append(float(r.imag)); zp.append(complex(mp.zeta(r, derivative=1)))
    gam = np.array(gam); zp = np.array(zp); np.savez(O + 'zeros100.npz', gam=gam, zp=zp)
print('gamma_1..5', np.round(gam[:5], 4), '|zeta\'(rho)| 1..5', np.round(np.abs(zp[:5]), 4))
for name, y in (('C2', C2), ('Q=C2-(J-Sig)', Q), ('kappa', kap), ('R', R), ('J-Sig', JmS)):
    dy = np.diff(y); print(f'{name}: mean(1/N-weighted) {np.sum(y/N)/np.sum(1/N):+.5f} sd {y.std():.5f}; |dy| rms {np.sqrt(np.mean(dy**2)):.2e} (N in [100,500]: {np.sqrt(np.mean(dy[:400]**2)):.2e}, [3000,4000]: {np.sqrt(np.mean(dy[-1000:]**2)):.2e})')
# arithmetic noise: correlation of C2(N+1)-C2(N) with mu(N+1)^2 (new squarefree column)
# periodogram on a uniform u grid (linear interpolation), Hann taper
ug = np.linspace(u[0], u[-1], 8192); T = ug[-1] - ug[0]; win = np.hanning(len(ug))
om = np.linspace(0.5, 80, 1600)
for name, y in (('C2', C2), ('Q', Q), ('kappa', kap)):
    yg = np.interp(ug, u, y); yg = yg - np.sum(yg * win) / np.sum(win)
    F = np.array([abs(np.sum(win * yg * np.exp(-1j * w * ug))) for w in om]) * 2 / np.sum(win)
    pk = [i for i in range(1, len(om) - 1) if F[i] > F[i - 1] and F[i] > F[i + 1]]
    pk = sorted(pk, key=lambda i: -F[i])[:12]
    print(f'{name} periodogram peaks (omega, amp): ' + ', '.join(f'({om[i]:.2f}, {F[i]:.5f})' for i in sorted(pk)))
pairs = sorted([(gam[j] - gam[i], i + 1, j + 1) for i in range(10) for j in range(i + 1, 10)])
print('pair freqs (zeros <= 10) < 20:', ', '.join(f'{w:.2f}({i},{j})' for w, i, j in pairs if w < 20))
print(f'T = {T:.3f} (u in [{u[0]:.3f}, {u[-1]:.3f}]), resolution 2pi/T = {2*math.pi/T:.2f}')
