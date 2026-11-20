# b3cert falsification of Lemma B (exact-profile dual term) and of the residue value M(s), theta near 0.49, T <= 1e8.
# Based on Riemann/Deep/CeilAudit/code/b3_profile.py (copied, read-only); changes: complex beta kept in the profile amplitude,
# Hankel branches T_+ and T_- reported separately, ratio scanned over r = n/N0 (T(n) depends on (n,pq) only via n/pq).
import sys, math, time, numpy as np
from scipy.special import roots_genlaguerre, gamma as cgamma
NL = 80
def apm(X, d):
    t, wt = roots_genlaguerre(NL, d.real - 0.5)
    base = wt * np.exp(1j * d.imag * np.log(t)) / cgamma(d + 0.5)
    tt = t[None, :] / (2 * X[:, None])
    return ((1 + 1j * tt) ** (d - 0.5)) @ base, ((1 - 1j * tt) ** (d - 0.5)) @ base
def Bpm(X, d):   # the two Hankel parts of B = sin(pi d/2)J_d + cos(pi d/2)Y_d (chunked: memory)
    bp = np.empty(len(X), complex); bm = np.empty(len(X), complex)
    for i in range(0, len(X), 4000):
        x = X[i:i+4000]; ap, am = apm(x, d); pre = np.sqrt(2 / (np.pi * x)) * (-0.5j)
        bp[i:i+4000] = pre * np.exp(1j * (x - np.pi / 4)) * ap; bm[i:i+4000] = -pre * np.exp(-1j * (x - np.pi / 4)) * am
    return bp, bm
def profile_exact(T, H, beta, du, N):   # P(w) = pi^{-1/2} int e^{-y^2+(1-2beta)y/H} e^{2iwy} e^{-i(T-Hw)g(y)} dy, beta complex
    m = np.arange(N); u0 = -N * du / 2; u = u0 + m * du
    ok = u > -H * 0.999; uu = np.where(ok, u, 0)
    y = np.where(ok, 0.5 * H * np.log1p(uu / H), -1e9)
    amp = np.where(np.abs(y) < 9, np.exp(-y ** 2 + (1 - 2 * beta) * np.where(np.abs(y) < 9, y, 0) / H), 0) * np.where(ok, 0.5 / (1 + uu / H), 0)
    G = amp * np.exp(-1j * T * np.where(ok, uu / H - np.log1p(uu / H), 0))
    k = np.fft.fftfreq(N) * N; w = 2 * np.pi * k / (N * du)
    Phi = du * np.exp(1j * w * u0) * N * np.fft.ifft(G) / math.sqrt(math.pi)
    o = np.argsort(w); return w[o], Phi[o]
def brute(T, H, beta, w, dy=2e-4):
    y = np.arange(-9, 9, dy); g = np.expm1(2 * y / H) - 2 * y / H
    return np.array([np.sum(np.exp(-y ** 2 + (1 - 2 * beta) * y / H + 2j * wi * y - 1j * (T - H * wi) * g)) * dy / math.sqrt(math.pi) for wi in w])
def setup(T, H, beta):
    tau = T / H ** 2; umax = H * (math.exp(18 / H) - 1)
    du = min(2e-3, 2 * np.pi / (8 * max(tau, 0.05) * umax))
    N = 1 << int(math.ceil(math.log2(2 * np.pi / (0.02 * du))))
    w, Phi = profile_exact(T, H, beta, du, N)
    keep = np.abs(w) <= min(8 * math.sqrt(1 + 4 * tau ** 2) + 12, 0.5 * T / H)
    return w[keep].copy(), Phi[keep].copy(), N
def Tpm(T, H, npq, alpha, beta, w, Phi):   # T_+(n), T_-(n) for n/(pq) = npq
    z, d = alpha - beta, alpha + beta; dw = w[1] - w[0]; xi = (T - H * w) / (2 * np.pi)
    A = xi ** (-beta) * Phi * xi ** (-z / 2); bp, bm = Bpm(4 * np.pi * np.sqrt(npq * xi), d); c = -np.pi * (H / (2 * np.pi)) * dw
    return c * np.sum(A * bp), c * np.sum(A * bm)
if __name__ == '__main__':
    t0 = time.time()
    shifts = [(0.01, 0.02), (0.05 + 0.1j, -0.03 + 0.07j), (0.1 + 0.2j, -0.1 - 0.15j), (0.2 - 0.1j, 0.03 + 0.2j)]
    rs = np.exp(np.linspace(math.log(0.003), math.log(8.0), 29))
    jobs = [tuple(map(float, j.split(':'))) for j in sys.argv[1].split(',')] if len(sys.argv) > 1 else [(T, th) for T in [1e6, 1e7, 1e8] for th in [0.49, 0.495, 0.47, 0.45, 0.40]]
    for T, th in jobs:
        if True:
            H = T ** th; K = T / H; tau = T / H ** 2
            for (al, be) in shifts:
                al, be = complex(al), complex(be)
                w, Phi, N = setup(T, H, be)
                if th == 0.49 and al == shifts[1][0]:
                    ws = np.array([0.0, 1.3, -2.7, 0.4 * w[-1]]); ex = np.interp(ws, w, Phi.real) + 1j * np.interp(ws, w, Phi.imag)
                    print(f"  validate T={T:.0e} th={th}: max|FFT-brute|/max|P| = {np.abs(ex - brute(T, H, be, ws)).max() / np.abs(Phi).max():.1e}", flush=True)
                # residue check: M(w') = int P(x) x^{w'-1} dx vs (H/(2 sqrt pi))(T/2pi)^{w'-1}, w' = 1-alpha-beta and 1
                xi = (T - H * w) / (2 * np.pi); dw = w[1] - w[0]
                res = [abs((H / (2 * np.pi)) * dw * np.sum(Phi * xi ** (ww - 1)) / ((H / (2 * math.sqrt(math.pi))) * (T / (2 * np.pi)) ** (ww - 1)) - 1) for ww in (1 - al - be, 1.0 + 0j)]
                rat = []; worst = (0, 0, 0)
                for r in rs:
                    c_r = 2 * math.sqrt(2 * math.pi) * math.sqrt(r); X0 = K * c_r; npq = r * T / H ** 2
                    tp, tm = Tpm(T, H, npq, al, be, w, Phi); s = H * X0 ** -0.5
                    rat.append((r, abs(tp + tm) / s, abs(tp) / s, abs(tm) / s))
                a = np.array(rat); i = int(np.argmax(a[:, 1]))
                big = a[a[:, 0] >= 4][:, 1].max()
                print(f"T={T:.0e} th={th} tau={tau:.3g} K={K:.4g} shifts=({al},{be}): max_r |T|/(H X0^-1/2) = {a[i,1]:.3f} at r={a[i,0]:.3g}; "
                      f"max|T+|={a[:,2].max():.3f} max|T-|={a[:,3].max():.3f}; min_r|T|ratio(r<=1)={a[a[:,0]<=1][:,1].min():.3f}; max ratio r>=4: {big:.1e}; "
                      f"residue |M/M_pred-1| = {res[0]:.1e}, {res[1]:.1e} (H/T={H/T:.1e}); N={N}", flush=True)
    print(f"done {time.time()-t0:.0f}s", flush=True)
