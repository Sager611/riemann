# ceilaudit7 A: exact CUE off-diagonal shifted mean D(a,b) = <sum_gamma Z(g+a)Z(g+b)>/<Z^2> (normalised units)
import sys, numpy as np
def Dab(N, a, b):
    u = 2*np.pi*a/N; v = 2*np.pi*b/N
    c = (u+v)/2; dl = (u-v)/2; co = np.cos(dl)
    f = {2: np.exp(-1j*c), -2: np.exp(1j*c), 1: -2*np.exp(-1j*c)-2*co, -1: -2*np.exp(1j*c)-2*co, 0: 4*co+2*np.cos(c)}
    M = N-1; T = np.zeros((M, M), complex); idx = np.arange(M)
    for m, val in f.items():
        j = idx[(idx-m >= 0) & (idx-m < M)]
        T[j, j-m] = val
    s, ld = np.linalg.slogdet(T)
    det = s*np.exp(ld)
    S = 4*np.sin(u/2)*np.sin(v/2)*det
    return (S/(N*(N+1))).real, abs(det.imag)/abs(det)
mode = sys.argv[1]
if mode == 'check':
    for N in (50, 200, 800):
        for a in (0.25, 0.5, 1.0, 1.5):
            d, im = Dab(N, a, a)
            print(f"N={N} a={a} D(a,a)={d:.6f} 1-sinc^2={1-(np.sinc(a))**2:.6f} im={im:.1e}")
        d, _ = Dab(N, 1e-3, 1e-3); print(f"N={N} D(e,e)/e^2={d/1e-6:.5f} pi^2/3={np.pi**2/3:.5f}")
elif mode == 'anti':
    N = int(sys.argv[2])
    for a in np.arange(0.05, 1.551, 0.05):
        d, im = Dab(N, a, -a); e, _ = Dab(N, a, 1e-4)
        print(f"a={a:.2f} D(a,-a)={d:+.6f} E(a)~D(a,eps)/eps={e/1e-4:+.6f} im={im:.1e}")
elif mode == 'grid':
    N = int(sys.argv[2]); g = np.arange(0.05, 1.201, 0.05)
    best = None
    for a in g:
        row = []
        for b in g:
            d, _ = Dab(N, a, -b); row.append(d)
            if d > 0 and (best is None or max(a, b) < best[0]): best = (max(a, b), a, b, d)
        print(f"a={a:.2f} " + " ".join(f"{x:+.3f}" for x in row))
    print("smallest max(a,b) with D(a,-b)>0:", best)
