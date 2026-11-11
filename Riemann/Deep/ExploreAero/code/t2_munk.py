# T2: Munk / Glauert (lifting-line) structure of Weil's window form; ground state vs elliptic loading and Riemann's Phi;
# tip load (Kutta-type edge behaviour); FT zeros; moved-zero adversary; Gershgorin ("diagonal dominance") certificate windows.
import numpy as np, sys
sys.path.insert(0, '.')
from weilleg import WeilLeg
from scipy.special import iv, spherical_jn
GAM = [14.134725141734693, 21.022039638771555, 25.010857580145688, 30.424876125859513, 32.935061587739189,
       37.586178158825671, 40.918719012147495, 43.327073280914999, 48.005150881167159, 49.773832477672302]
def Phi(u):
    u = np.asarray(u, float); n = np.arange(1, 40)[:, None]; e = np.exp(4*u)[None, :]
    return (2*np.pi**2*n**4*np.exp(9*u)[None, :] - 3*np.pi*n**2*np.exp(5*u)[None, :])*np.exp(-np.pi*n**2*e) if False else \
        ((2*np.pi**2*n**4*np.exp(9*u)[None, :] - 3*np.pi*n**2*np.exp(5*u)[None, :])*np.exp(-np.pi*n**2*e)).sum(0)
def glauert(W, M):   # Legendre coefficients S (K x M) of phi_k = sin(k theta), x = a cos theta, k = 1,3,..; exact Gram G_G
    N = 2*W.K + 2*M + 20; j = np.arange(1, N + 1); th = j*np.pi/(N + 1); t = np.cos(th); w = np.pi/(N + 1)*np.sin(th)**2
    ks = 2*np.arange(M) + 1
    U = np.array([np.sin(k*th)/np.sin(th) for k in ks])                  # U_{k-1}(t)
    Pv = W.vals(W.a*t)                                                   # K x N, includes nrm
    S = W.a*(Pv*w) @ U.T
    kk, ll = np.meshgrid(ks, ks, indexing='ij')
    GG = W.a*(1/(1 - (kk - ll)**2) - 1/(1 - (kk + ll)**2))
    return ks, S, GG
def lowest(Q, G=None):
    if G is None: return np.linalg.eigh(Q)
    Lc = np.linalg.cholesky(G); Li = np.linalg.inv(Lc)
    e, V = np.linalg.eigh(Li @ Q @ Li.T); return e, Li.T @ V
def gersh(Q, G=None):  # min over rows of Q_ii - sum_{j!=i}|Q_ij| in the G-orthonormalised basis (Lowdin)
    if G is not None:
        e, V = np.linalg.eigh(G); Gi = V @ np.diag(e**-0.5) @ V.T; Q = Gi @ Q @ Gi
    return min(Q[i, i] - (np.abs(Q[i]).sum() - abs(Q[i, i])) for i in range(len(Q)))
K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
print(f"K = {K} even Legendre modes (degree <= {2*K-2}); Glauert projections use the same K")
print(" a      c      lam0        lam1       tip(K/2)  tip(K)   cos(f0,Phi) cos(f0,ell) GlauertRes(M=1,4,8)       Munk: Dmin   D_ell/Dmin  2+Dmin   lam(R|kerl)")
rows = []
for a in [0.1, 0.2, 0.3, np.log(2)/2, 0.4, 0.45, 0.5, np.log(3)/2, 0.6, 0.65]:
    W = WeilLeg(a, K).build(); e, V = np.linalg.eigh(W.Q); c0 = V[:, 0]*np.sign(V[0, 0])
    Wh = WeilLeg(a, K//2).build(); eh, Vh = np.linalg.eigh(Wh.Q); ch = Vh[:, 0]*np.sign(Vh[0, 0])
    tipK = (c0 @ W.nrm)/(c0 @ (W.nrm*W.vals(np.array([0.0]))[:, 0]/W.nrm)); tiph = (ch @ Wh.nrm)/(ch @ (Wh.vals(np.array([0.0]))[:, 0]))
    # Phi|_W and elliptic loading in Legendre coefficients (Gauss-Legendre projection)
    xg, wg = np.polynomial.legendre.leggauss(400); x = a*xg; wx = a*wg; Pv = W.vals(x)
    cPhi = (Pv*wx) @ Phi(x); cPhi /= np.linalg.norm(cPhi)
    ks, S, GG = glauert(W, 8); s1 = S[:, 0]
    cosPhi = abs(c0 @ cPhi); cosEll = abs(c0 @ s1)/np.linalg.norm(s1)
    res = []
    for M in [1, 4, 8]:
        b = np.linalg.solve(GG[:M, :M], S[:, :M].T @ c0); res.append(np.sqrt(max(0.0, 1 - b @ GG[:M, :M] @ b)))
    # Munk problem: R = A - P; D_min = min{c^T R c : l.c = 1}; elliptic loading's drag at unit lift
    R = W.A - W.P; l = W.l
    Pk = np.eye(K) - np.outer(l, l)/(l @ l); eR = np.linalg.eigvalsh(Pk @ R @ Pk)
    lamker = sorted(eR, key=abs)[1:]; lamker = min(v for v in eR if abs(v) > 1e-12) if True else 0
    # constrained min via KKT
    Mk = np.block([[2*R, l[:, None]], [l[None, :], np.zeros((1, 1))]]); sol = np.linalg.solve(Mk, np.r_[np.zeros(K), 1.0])
    cm = sol[:K]; Dmin = cm @ R @ cm
    Dell = (s1 @ R @ s1)/(l @ s1)**2
    print(f"{a:5.3f} {np.exp(2*a):6.3f} {e[0]:11.4e} {e[1]:10.3e} {tiph:8.4f} {tipK:8.4f}  {cosPhi:9.6f}  {cosEll:9.6f}  " +
          ' '.join(f"{r:7.4f}" for r in res) + f"   {Dmin:10.5f} {Dell/Dmin if Dmin!=0 else 0:9.4f} {2+Dmin:10.3e} {lamker:10.3e}")
    rows.append((a, W, c0, e))
# Glauert-basis structure of the archimedean operator at a = 0.5 (lifting line: (pi/2) k delta_kl exactly)
for a in [0.3, 0.6]:
    W = WeilLeg(a, 60).build(); ks, S, GG = glauert(W, 6)
    AG = S.T @ W.A @ S; SG = S.T @ S
    d = np.sqrt(np.abs(np.diag(AG)))
    print(f"a={a}: Glauert projection check max|S^T S - G_G|/max G_G = {np.abs(SG - GG).max()/np.abs(GG).max():.2e}")
    print("  Weil arch A in Glauert modes k=1,3,..,11, normalised |A_kl|/sqrt(|A_kk A_ll|):")
    for i in range(6): print('   ', ' '.join(f"{abs(AG[i, j])/(d[i]*d[j]):6.3f}" for j in range(6)), f"  A_kk={AG[i, i]:8.4f}  (lifting line (pi/2)k = {np.pi/2*ks[i]:6.3f})")
    eGG = np.linalg.eigvalsh(GG); print(f"  Glauert Gram G_G (also off-diagonal) eigenvalue range {eGG.min():.3e}..{eGG.max():.3e}")
# FT zeros of the ground state vs zeta zeros
print("Real zeros of F0(r) = FT of ground state on [0, 52] vs gamma_k:")
for a, W, c0, e in rows:
    if a < 0.29: continue
    r = np.linspace(0.01, 52, 52001); F = c0 @ W.ft(r); z = r[:-1][np.sign(F[:-1]) != np.sign(F[1:])]
    print(f"  a={a:.3f}: zeros {np.round(z, 3)[:7]}")
print("gamma:", np.round(GAM[:7], 3))
# adversary: move one zero pair off the line: rho = 1/2 +- delta + i gamma (quartet); Q_mv = Q - 2 v(g)v(g)^T + 4 Re v(g+id)v(g+id)^T
print("Moved-zero adversary: lam0(Q_moved) (negative = moved function fails positivity on this window)")
for gi in [0, 1, 4]:
    g = GAM[gi]
    for d in [0.02, 0.1, 0.25]:
        out = []
        for a, W, c0, e in rows:
            v = W.ft(np.array([g]))[:, 0].real; vz = W.ft(np.array([g + 1j*d]))[:, 0]
            Qm = W.Q - 2*np.outer(v, v) + 4*np.real(np.outer(vz, vz))
            out.append(f"{np.linalg.eigvalsh(Qm)[0]:10.2e}")
        print(f"  gamma_{gi+1}={g:7.3f} delta={d:4.2f}: " + ' '.join(out))
print("  windows a:", ' '.join(f"{a:10.3f}" for a, *_ in rows))
print("zeta lam0:  ", ' '.join(f"{e[0]:10.2e}" for a, W, c0, e in rows))
# Gershgorin (diagonal-dominance, 'Munk-type') certificates in Legendre and Glauert(+pole) bases
print("Gershgorin margin min_i(Q_ii - sum|Q_ij|) (positive = diagonal dominance certifies Q>0 on the span):")
for a, W, c0, e in rows:
    ks, S, GG = glauert(W, 8)
    QG = S.T @ W.Q @ S
    print(f"  a={a:.3f}: Legendre K={K}: {gersh(W.Q):10.3e}   Glauert M=8 (Lowdin): {gersh(QG, S.T @ S):10.3e}")
