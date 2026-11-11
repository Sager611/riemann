# cue_visc.py -- jerbyconf: random-matrix (CUE) model of Gram-point viscosity at confined isolated bad Gram points.
# Eigenangles theta_j of a Haar unitary N x N; Gram lattice x_m = 2 pi m / N (smooth count = m, as theta(g_n)/pi = n); S_m = #{theta_j in (0, x_m]} - m;
# good iff S_m even (zeta: (-1)^n Z(g_n) > 0 iff S(g_n) even); isolated bad: S_m odd, S_{m-1}, S_{m+1} even; confined: exactly two angles in (x_{m-1}, x_{m+1}].
# Viscosity: Z_U'/Z_U(x) = sum_j (1/2) cot((x - theta_j)/2); normalised s*mu = (2 pi/N) |Z_U'/Z_U(x_m)| (zeta: s = 2 pi/log(t/2pi); the zeta drift
# pi/4 is cancelled to +-0.03 by the far zeros, logs/decomp_low.txt).  usage: cue_visc.py N MATRICES SEED
import sys, numpy as np, time
N = int(sys.argv[1]); M = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
x = 2*np.pi*np.arange(N)/N
conf, zf, iso_tot, gp = [], [], 0, 0
t0 = time.time()
for it in range(M):
    A = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2)
    Q, R = np.linalg.qr(A); Q = Q*(np.diag(R)/np.abs(np.diag(R)))
    th = np.sort(np.mod(np.angle(np.linalg.eigvals(Q)), 2*np.pi))
    cnt = np.searchsorted(th, x, side='right')          # #{theta <= x_m}
    S = cnt - np.arange(N)
    good = (S % 2 == 0)
    gl = np.roll(good, 1); gr = np.roll(good, -1)
    isob = (~good) & gl & gr
    isob[0] = isob[-1] = False                            # avoid wrap-around at the cut
    idx = np.nonzero(isob)[0]
    if len(idx) == 0: continue
    nblk = cnt[idx + 1] - cnt[idx - 1]                    # zeros in (x_{m-1}, x_{m+1}]
    v = 0.5/np.tan((x[idx][:, None] - th[None, :])/2)
    smu = (2*np.pi/N)*np.abs(v.sum(axis=1))
    conf.append(smu[nblk == 2]); zf.append(smu[nblk == 0]); iso_tot += len(idx); gp += N - 2
conf = np.concatenate(conf); zf = np.concatenate(zf)
print(f'N={N} matrices={M} lattice points={gp} isolated bad={iso_tot} ({iso_tot/gp:.4f} per point) confined={len(conf)} zero-free={len(zf)} ({len(zf)/iso_tot:.2e}) secs={time.time()-t0:.0f}')
q = np.quantile(conf, [1e-4, 1e-3, 1e-2, 0.1, 0.5])
print('confined s*mu quantiles 1e-4,1e-3,1e-2,0.1,0.5:', ' '.join(f'{v:.3f}' for v in q), ' min', f'{conf.min():.3f}')
for L in (2.5, 2.2, 2.0, 1.948, 1.8, 1.7, 1.5, 1.2, 1.0):
    print(f'  P(confined and s*mu < {L}) per isolated bad = {np.sum(conf < L)/iso_tot:.3e}  (count {np.sum(conf < L)})')
if len(zf): print('zero-free s*mu min', f'{zf.min():.3f}', 'median', f'{np.median(zf):.3f}')
np.save(sys.argv[4] if len(sys.argv) > 4 else '/tmp/jc_cue.npy', conf)
