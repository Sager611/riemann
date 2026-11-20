# nb2fals f5: (i) head constant: fits of Sigma_K on K in [40, NMAX] (Sigma_inf - A/log K; Sigma_inf - A K^-theta), J - Sigma_NMAX;
# (ii) the one-sided combination C2 <= (J - Sigma_NMAX) + Qbar with Q = C2 - (J - Sigma_N) from the GP fits (fitted shape and B), trend and none.
import os, sys, math, numpy as np
os.environ.setdefault('NB2F_NMAX', '4000'); sys.argv = ['f3', 'none']
exec(open('/Users/adriansager/playground/riemann/Riemann/Deep/NB2Fals/code/f3_robust.py').read())
J = 3.5626634565363677; Sig = np.load(O + f'perN{NMX}.npz')['Sig']
KK = np.arange(40, NMX + 1); S = Sig[KK - 1]
def lsq(X): coef, *_ = np.linalg.lstsq(X, S, rcond=None); return coef, float(np.sqrt(np.mean((X @ coef - S) ** 2)))
c1, e1 = lsq(np.column_stack([np.ones(len(KK)), 1 / np.log(KK)]))
best = min(((th,) + lsq(np.column_stack([np.ones(len(KK)), KK ** (-th)])) for th in np.linspace(0.05, 1.5, 59)), key=lambda t: t[2])
c3, e3 = lsq(np.column_stack([np.ones(len(KK)), np.log(np.log(KK))]))
c4, e4 = lsq(np.column_stack([np.ones(len(KK)), 1 / np.log(KK), 1 / np.log(KK) ** 2]))
print(f'(i) Sigma_K, K in [40,{NMX}]: Sigma_{NMX} = {Sig[-1]:.6f}, J - Sigma_{NMX} = {J - Sig[-1]:.5f}; Sigma_1280 {Sig[1279]:.6f} Sigma_2560 {Sig[2559]:.6f}')
print(f'    S_inf - A/log K: S_inf {c1[0]:.5f} (J - S_inf {J - c1[0]:.5f}), A {-c1[1]:.4f}, rms {e1:.1e}')
print(f'    S_inf - A/log K - B/log^2 K: S_inf {c4[0]:.5f} (J - S_inf {J - c4[0]:.5f}), rms {e4:.1e}')
print(f'    S_inf - A K^-th: th {best[0]:.3f}, S_inf {best[1][0]:.5f} (J - S_inf {J - best[1][0]:.5f}), rms {best[2]:.1e}')
print(f'    a + b log log K (divergent): b {c3[1]:.4f}, rms {e3:.1e}')
for trend in (0, 1):
    for fs in (False, True):
        b, s, ll, x = fit('Q', 50, trend, fs)
        print(f'(ii) Q{" + C3/L" if trend else ""} K=50 shape {"fitted" if fs else "B"}: Qbar = {b[0]:+.5f} +- {s[0]:.5f}' + (f' (alpha {x[3]:.2f} beta {x[4]:.2f})' if fs else '') + f' => C2 <= {J - Sig[-1] + b[0]:+.5f} +- {s[0]:.5f} ({(J - Sig[-1] + b[0]) / s[0]:+.1f} sigma); with J - S_inf from 1/log K fit: {J - c1[0] + b[0]:+.5f}', flush=True)
