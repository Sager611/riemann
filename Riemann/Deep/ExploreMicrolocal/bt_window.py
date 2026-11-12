# explore47 scale test: CCM/Bérard bound pieces on Odlyzko windows. Usage: python3 bt_window.py <file> [k_lo k_hi]
import sys, os, math, time, json
import numpy as np
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from btlib import prime_powers, load, TB
t0 = time.time()
name = sys.argv[1]; T0, k0, vals = load(name)
lo, hi = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, len(vals))
vals = vals[lo:hi]; K = len(vals); kidx = np.arange(k0 + lo + 1, k0 + lo + K + 1)
mp.mp.dps = 50
tau_d = np.array([float(mp.mpf(v)) for v in vals])            # γ − T0 (double, |·| ≤ 8e4)
th = np.array([float(mp.siegeltheta(T0 + mp.mpf(v)) / mp.pi - (k0 + lo + 1 + i)) for i, v in enumerate(vals)])
Splus = -th - 1.0                                             # S(γ_k^+) = k − ϑ/π − 1
Sminus = Splus - 1.0                                          # S(γ_k^−)
tmid = float(T0) + float(np.mean(tau_d)); L = math.log(tmid / (2 * math.pi))
DL = list(range(2, 15)); NMAX = int(math.exp(14)) + 1
n, lam, base = prime_powers(NMAX); logn = np.log(n.astype(float)); N = len(n)
A = np.zeros((N, len(DL))); B = np.zeros_like(A)
for j, Dt in enumerate(DL):
    m = logn < Dt; e = logn[m] / Dt
    A[m, j] = (1 / math.pi + (1 - e) / np.tan(math.pi * e)) * lam[m] / np.sqrt(n[m]) / Dt
    B[m, j] = (1 - e) * lam[m] / np.sqrt(n[m]) / Dt
if T0 != 0:
    tp = 2 * mp.pi
    ph0 = np.array([float(mp.fmod(T0 * mp.log(int(q)), tp)) for q in n])
else:
    ph0 = np.zeros(N)
print(f'{name} K={K} t≈{tmid:.6g} log(t/2π)={L:.4f} N_pp={N} prep {time.time()-t0:.1f}s', flush=True)
eps = 1e-7
def orbit(tvec, phases_fn):
    Op = np.zeros((len(tvec), len(DL))); Om = np.zeros_like(Op)
    for i in range(0, len(tvec), 32):
        ph = phases_fn(tvec[i:i + 32])
        s = np.sin(ph); c = np.cos(ph); sa = s @ A; cb = c @ B
        Op[i:i + 32] = -(sa + cb); Om[i:i + 32] = -(sa - cb)
    return Op, Om
step = max(1, K // 2500); sel = np.arange(0, K, step)
zp = lambda tv: ph0[None, :] + np.outer(tv, logn)
Opp, Omp = orbit(tau_d[sel] + eps, zp)      # at γ_k^+
Opm, Omm = orbit(tau_d[sel] - eps, zp)      # at γ_k^−
tauD = np.array([math.log((tmid) / (2 * math.pi)) / (2 * Dt) for Dt in DL])
Lk = np.log((float(T0) + tau_d[sel]) / (2 * math.pi)); tk = Lk[:, None] / (2 * np.array(DL, dtype=float))[None, :]
Bup = tk + Opp; Blo = -tk + Omm
marg_up = (Bup - Splus[sel, None]).min(axis=0); marg_lo = (Sminus[sel, None] - Blo).min(axis=0)
# random-phase (Bohr–Jessen) model: fresh iid phases θ_p per prime for each sample, phase of p^k = kθ_p
rng = np.random.default_rng(47); pb = np.unique(base)
kk = np.round(logn / np.log(base.astype(float))).astype(int); bi = np.searchsorted(pb, base)
def rphases(tv):
    th_p = rng.uniform(0, 2 * np.pi, (len(tv), len(pb)))
    return th_p[:, bi] * kk[None, :]
Rp, Rm = orbit(np.zeros(len(sel)), rphases)
Dfine = np.arange(1.0, 17.01, 0.25)
n2, lam2, _ = prime_powers(int(math.exp(17)) + 1); logn2 = np.log(n2.astype(float))
TBf = np.array([TB(Dt, logn2, lam2, n2) for Dt in Dfine])
tot = L / (2 * Dfine) + TBf; jm = int(np.argmin(tot))
res = dict(name=name, K=K, lo=lo, t=tmid, L=L, Dtri=float(Dfine[jm]), Btri=float(tot[jm]),
           ccm=0.25 * math.log(tmid) / math.log(math.log(tmid)), maxS=float(Splus.max()), minS=float(Sminus.min()),
           DL=DL, TB=[TB(Dt, logn2, lam2, n2) for Dt in DL], tau=tauD.tolist(),
           supO_zeta=np.maximum(np.abs(Opp).max(axis=0), np.abs(Omm).max(axis=0)).tolist(),
           supO_rand=np.maximum(np.abs(Rp).max(axis=0), np.abs(Rm).max(axis=0)).tolist(),
           rmsO_zeta=np.sqrt((Opp ** 2).mean(axis=0)).tolist(), rmsO_rand=np.sqrt((Rp ** 2).mean(axis=0)).tolist(),
           maxBup=Bup.max(axis=0).tolist(), minBlo=Blo.min(axis=0).tolist(), marg_up=marg_up.tolist(), marg_lo=marg_lo.tolist(),
           nsel=int(len(sel)), secs=round(time.time() - t0, 1))
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs'), exist_ok=True)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs', f'bt_{name}_{lo}.json'), 'w') as f: json.dump(res, f)
print(f"t={tmid:.4g} Δ_tri={res['Dtri']} B_tri={res['Btri']:.3f} CCM¼={res['ccm']:.3f} maxS={res['maxS']:.3f} minS={res['minS']:.3f}")
for j, Dt in enumerate(DL):
    print(f"Δ={Dt:2d} τ={tauD[j]:.3f} TB={res['TB'][j]:9.3f} sup|O|ζ={res['supO_zeta'][j]:.3f} rand={res['supO_rand'][j]:.3f} "
          f"rms ζ={res['rmsO_zeta'][j]:.3f} rand={res['rmsO_rand'][j]:.3f} maxB+={res['maxBup'][j]:.3f} minB-={res['minBlo'][j]:.3f} "
          f"marg+={marg_up[j]:.3f} marg-={marg_lo[j]:.3f}")
print('secs', res['secs'], flush=True)
