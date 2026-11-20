# explore62 (a): exact CUE fourth-order data vs the squared-chaos (MNPP) and hypercontractivity bounds.
# K_N(phi) = E|Lambda(1)|^2|Lambda(e^{i phi})|^2 / (E|Lambda|^2)^2 = D_N(|1-e^{it}|^2 |1-e^{i(t-phi)}|^2)/(N+1)^2 (Heine), G_N = K_N/K_N(0).
import numpy as np
from fractions import Fraction
def logD(N, phi):
    q = np.exp(-1j*phi); f = {2: q, -2: np.conj(q), 1: -2*q-2, -1: -2*np.conj(q)-2, 0: 4+2*np.cos(phi)}
    T = np.zeros((N, N), complex); idx = np.arange(N)
    for m, v in f.items():
        j = idx[(idx-m >= 0) & (idx-m < N)]; T[j, j-m] = v
    s, ld = np.linalg.slogdet(T); return ld, s
def D_exact_pi(N):  # symbol 2-2cos(2t): Toeplitz splits into even/odd tridiagonal (2,-1) blocks, det(m)=m+1
    def tri(m):
        a, b = 1, 2  # det of size 0, 1
        for _ in range(m-1): a, b = b, 2*b-a
        return b if m >= 1 else 1
    return tri((N+1)//2)*tri(N//2)
for N in (1, 2, 3, 10, 51, 100):
    ex = D_exact_pi(N); fl = np.exp(logD(N, np.pi)[0]) if N > 1 else None
    print(f"N={N}: D_N(pi) exact={ex} = (ceil(N/2)+1)(floor(N/2)+1)? {ex == ((N+1)//2+1)*(N//2+1)}; float={fl}; K_N(pi)={Fraction(ex,(N+1)**2)} = {ex/(N+1)**2:.6f}")
for N in (100, 200, 400):
    l0, _ = logD(N, 0.0)
    G = []
    for x in (0.5, 1, 2, 4, 8, 16):
        lx, s = logD(N, 2*np.pi*x/N); G.append((x, (s*np.exp(lx-l0)).real))
    M4 = (N+2)**2*(N+3)/(12*(N+1))
    print(f"N={N} M4={M4:.1f} " + " ".join(f"G({x})={g:.4f}" for x, g in G) + "  x^2 G(x): " + " ".join(f"{x*x*g:.4f}" for x, g in G[3:]) + f"  (3/pi^2={3/np.pi**2:.4f})")
    ks = []
    for phi in np.linspace(0.05, np.pi, 64):
        lx, s = logD(N, phi); ks.append((phi, (s*np.exp(lx)).real/(N+1)**2))
    cross = min(p for p, k in ks if k < 1)
    print(f"   K_N(phi) < 1 from phi = {cross:.3f} (FH limit pi/3 = {np.pi/3:.3f}); K_N(pi) = {ks[-1][1]:.5f}; FH 1/(4 sin^2(phi/2)) at pi/2: {1/(4*np.sin(np.pi/4)**2):.4f} vs {[k for p,k in ks if abs(p-np.pi/2)<0.03]}")
# chaos bounds: q=1 (Gaussian) G>=1/3; q=2 G>=1/15 (max kurtosis 15); q-th chaos G>=9^-q (hypercontractivity + MNPP)
for q, lb in ((1, 1/3), (2, 1/15), (3, 9.0**-3), (4, 9.0**-4)):
    xq = np.sqrt(3/np.pi**2/lb)
    print(f"q={q}: fixed-chaos bound G>={lb:.5f}; zeta/CUE G(x)~3/(pi^2 x^2) falls below it for x > {xq:.2f}")
