# b3nu: exact dual transform T(x) = T_+(x) + T_-(x), x = n/(pq), zero shifts (alpha = beta = 0, d = z = 0, B = Y_0).
# Profile by FFT copied from Riemann/Deep/B3Cert/code/b3c_fals.py (read-only); Hankel parts by the d = 0 asymptotic series (X >= 100).
# Demodulation: D_s(x) = T_s(x) exp(-i Phi_s(x)), Phi_s(x) = T[log(1+v) + v/(1+v)], v^2 = (2 pi x/T)(1+v), sign v = s
# (2D stationary value of the joint phase 2Ty/H - 2pi(e^{2y/H}-1)xi + s 4pi sqrt(x xi)); D_s is smooth -> spline in log x.
import math, numpy as np
from scipy.interpolate import CubicSpline
def profile_exact(T, H, du, N):
    m = np.arange(N); u0 = -N * du / 2; u = u0 + m * du
    ok = u > -H * 0.999; uu = np.where(ok, u, 0)
    y = np.where(ok, 0.5 * H * np.log1p(uu / H), -1e9)
    amp = np.where(np.abs(y) < 9, np.exp(-y ** 2 + np.where(np.abs(y) < 9, y, 0) / H), 0) * np.where(ok, 0.5 / (1 + uu / H), 0)
    G = amp * np.exp(-1j * T * np.where(ok, uu / H - np.log1p(uu / H), 0))
    k = np.fft.fftfreq(N) * N; w = 2 * np.pi * k / (N * du)
    Phi = du * np.exp(1j * w * u0) * N * np.fft.ifft(G) / math.sqrt(math.pi)
    o = np.argsort(w); return w[o], Phi[o]
def setup(T, H, sub_dw=0.04):
    tau = T / H ** 2; umax = H * (math.exp(18 / H) - 1)
    du = min(2e-3, 2 * np.pi / (8 * max(tau, 0.05) * umax))
    N = 1 << int(math.ceil(math.log2(2 * np.pi / (0.02 * du))))
    w, Phi = profile_exact(T, H, du, N)
    keep = np.abs(w) <= min(8 * math.sqrt(1 + 4 * tau ** 2) + 12, 0.5 * T / H)
    w, Phi = w[keep].copy(), Phi[keep].copy()
    step = max(1, int(sub_dw / (w[1] - w[0])))
    return w, Phi, w[::step].copy(), Phi[::step].copy()
def hank(X, s):   # s = +1: sqrt(2/(pi X))(-i/2)e^{i(X-pi/4)}a_+ ; s = -1: -sqrt(2/(pi X))(-i/2)e^{-i(X-pi/4)}a_-   (Y_0 = b_+ + b_-)
    a = 1 - 1j * s / (8 * X) - 9 / (128 * X ** 2) + 1j * s * 225 / (3072 * X ** 3) + 11025 / (98304 * X ** 4)
    return s * np.sqrt(2 / (np.pi * X)) * (-0.5j) * np.exp(1j * s * (X - np.pi / 4)) * a
def Tbr(T, H, xs, w, Phi, s, chunk=None):
    chunk = chunk or max(8, int(4e5 // len(w))); dw = w[1] - w[0]; xi = (T - H * w) / (2 * np.pi); c = -np.pi * (H / (2 * np.pi)) * dw; out = np.empty(len(xs), complex)
    for i in range(0, len(xs), chunk):
        x = xs[i:i + chunk]; X = 4 * np.pi * np.sqrt(np.outer(x, xi)); out[i:i + chunk] = c * (hank(X, s) @ Phi)
    return out
def Phis(T, x, s):
    b = 2 * np.pi * x / T; v = (b + s * np.sqrt(b * b + 4 * b)) / 2
    return T * np.log1p(v) + T * v / (1 + v)
class Tinterp:
    def __init__(self, T, H, xmin, xmax, npts=1500):
        self.T = T; w, Phi, ws, Ps = setup(T, H); self.full = (w, Phi); self.sub = (ws, Ps)
        self.lx = np.linspace(math.log(xmin), math.log(xmax), npts); x = np.exp(self.lx); self.sp = {}
        for s in (1, -1):
            D = Tbr(T, H, x, ws, Ps, s) * np.exp(-1j * Phis(T, x, s)) * x ** 0.25
            self.sp[s] = (CubicSpline(self.lx, D.real), CubicSpline(self.lx, D.imag))
    def __call__(self, x, s=None):
        lx = np.log(x); out = 0
        for t in ((1, -1) if s is None else (s,)):
            a, b = self.sp[t]; out = out + (a(lx) + 1j * b(lx)) * x ** -0.25 * np.exp(1j * Phis(self.T, x, t))
        return out
