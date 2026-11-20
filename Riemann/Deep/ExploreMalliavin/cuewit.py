# explore62 (d): Theorem 3's witness at bounded scale: Gaussian background (flat spectrum on [-pi,pi], FFT synthesis) + K independent CUE(N) windows
# (real char. poly, normalised <W^2>=1, amplitude lambda = eta^-1/2, length fraction p = eta^2, C^1 cosine junctions of 1 unit).
# Measures zeros per unit, second-order means vs sinc, G(x) vs the exact CUE G_N (Toeplitz), kurtosis vs M4_CUE(N).
import numpy as np, sys
rng = np.random.default_rng(62)
N, LOG2N, eta = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]); S = 16; h = 1.0/S
n = 2**LOG2N; K = int(round(eta**2*n*h/N)); p = K*N/(n*h); eta = np.sqrt(p)  # power-of-2 FFT (Bluestein buffers broke the cap)
lam = eta**-0.5; cg = (1-eta)/(1-p)
freq = 2*np.pi*np.fft.rfftfreq(n, d=h); amp = (np.abs(freq) <= np.pi).astype(float)
Zc = (rng.standard_normal(len(freq)) + 1j*rng.standard_normal(len(freq)))*amp
g = np.fft.irfft(Zc, n); g *= np.sqrt(cg/np.mean(g*g)); del Zc, amp, freq
f = g.copy(); del g
def cue(N):
    A = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2); Q, R = np.linalg.qr(A)
    Q = Q*(np.diag(R)/np.abs(np.diag(R))); return np.angle(np.linalg.eigvals(Q))
x = np.arange(N*S)*h; ph = 2*np.pi*x/N
ramp = np.ones(N*S); r1 = np.arange(S)*h; ramp[:S] = 0.5 - 0.5*np.cos(np.pi*r1); ramp[-S:] = ramp[:S][::-1]
slots = np.arange(0, n - N*S, N*S + 2*S); starts = np.sort(rng.choice(slots, K, replace=False))
w4 = []
for s0 in starts:
    th = cue(N); W = np.prod(2*np.sin((ph[:, None] - th[None, :])/2), axis=1)/np.sqrt(N+1); w4.append(np.mean(W**4))
    f[s0:s0+N*S] = ramp*lam*W + (1-ramp)*f[s0:s0+N*S]
w4 = np.array(w4); print(f"windows: empirical <W^4> = {w4.mean():.1f} (ensemble M4_CUE(N) = {(N+2)**2*(N+3)/(12*(N+1)):.1f}); median per window {np.median(w4):.1f}, max {w4.max():.0f}: heavy tail (moments of moments)")
zeros = int((np.sign(f[:-1])*np.sign(f[1:]) < 0).sum()); m2 = np.mean(f*f); m4 = np.mean(f**4)
print(f"N={N} K={K} eta={eta} p={p:.4f} lambda={lam:.3f}: length {n*h:.0f} units, zeros per unit {zeros/(n*h):.4f} "
      f"(predicted (1-p)/sqrt3 + p = {(1-p)/np.sqrt(3)+p:.4f}); <f^2> = {m2:.4f}; M4 = {m4/m2**2:.2f} (M4_CUE(N) = {(N+2)**2*(N+3)/(12*(N+1)):.2f})")
print("  C(u)/<f^2> vs sinc: " + " ".join(f"u={u}: {np.mean(f[:n-int(u*S)]*f[int(u*S):])/m2:+.4f}/{np.sinc(u):+.4f}" for u in (0.25, 0.5, 1, 1.5, 2.5)))
def logD(N, phi):  # = cue4.logD (Toeplitz determinant of |1-e^{it}|^2|1-e^{i(t-phi)}|^2)
    q = np.exp(-1j*phi); fc = {2: q, -2: np.conj(q), 1: -2*q-2, -1: -2*np.conj(q)-2, 0: 4+2*np.cos(phi)}
    T = np.zeros((N, N), complex); idx = np.arange(N)
    for m, v in fc.items():
        j = idx[(idx-m >= 0) & (idx-m < N)]; T[j, j-m] = v
    s, ld = np.linalg.slogdet(T); return ld, s
l0, _ = logD(N, 0.0)
for xx in (0.5, 1, 2, 4):
    k = int(xx*S); Gf = np.mean(f[:n-k]**2*f[k:]**2)/m4
    lx, sg = logD(N, 2*np.pi*xx/N); GN = (sg*np.exp(lx-l0)).real
    M4N = (N+2)**2*(N+3)/(12*(N+1)); pred = (M4N*GN + (1-p)*cg**2*(1+2*np.sinc(xx)**2))/(M4N + 3*(1-p)*cg**2)
    print(f"  G_f({xx}) = {Gf:.4f}; exact CUE G_N = {GN:.4f}; formula (windows + Wick background) = {pred:.4f}")
