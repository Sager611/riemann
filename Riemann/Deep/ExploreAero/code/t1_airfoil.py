# T1: thin-airfoil operator = i tanh(pi H), H = dilation (Berry-Keating xp) generator in u=(1+x)/(1-x).
# (a) generalized eigenfunctions g_tau = ((1+x)/(1-x))^{i tau}/sqrt(1-x^2): T g = i tanh(pi tau) g, T f = (1/pi) PV int f(t)/(t-x) dt
# (b) Legendre-Galerkin (closed form) truncations: spectrum in tau = artanh(mu)/pi, counting law vs Riemann-von Mangoldt
# (c) discrete Kutta mode (exact zero eigenvector, K odd) vs projection of 1/sqrt(1-x^2)
import numpy as np, mpmath as mp, sys
mp.mp.dps = 30
def T_apply(f, x):
    # PV int_{-1}^{1} f(t)/(t-x) dt for f with integrable endpoint singularities: subtract f(x)
    fx = f(x)
    g = lambda t: (f(t) - fx)/(t - x)
    I = mp.quad(g, [-1, x, 1])
    I += fx*mp.log((1 - x)/(1 + x))
    return I/mp.pi
print("(a) generalized eigenfunctions")
for tau in [0, 0.35, 0.7, 1.5]:
    g = lambda t, tau=tau: mp.power((1+t)/(1-t), 1j*tau)/mp.sqrt(1 - t*t)
    for x in [mp.mpf('-0.6'), mp.mpf('0.3')]:
        lhs = T_apply(g, x); rhs = 1j*mp.tanh(mp.pi*tau)*g(x)
        print(f"  tau={tau:4} x={float(x):5}: T g / g = {complex(lhs/g(x)):.12f}  i tanh(pi tau) = {complex(1j*mp.tanh(mp.pi*tau)):.12f}")
# (b) Legendre-Galerkin matrix of T (orthonormal p_n): <p_m, T p_n>, T P_n = -(2/pi) Q_n (Neumann), int P_m Q_n = (1-(-1)^{m+n})/((m-n)(m+n+1))
def Tmat(K):
    m = np.arange(K)[:, None]; n = np.arange(K)[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        PQ = np.where(m == n, 0.0, (1 - (-1.0)**(m + n))/((m - n)*(m + n + 1.0)))
    return -(2/np.pi)*np.sqrt((2*m + 1)*(2*n + 1))/2*PQ
print("(b) spectrum of K x K Legendre-Galerkin truncation, tau_k = artanh(mu_k)/pi (mu_k>0)")
for K in [41, 81, 161, 321, 641]:
    M = Tmat(K)
    assert np.allclose(M, -M.T)
    mu = np.sort(np.linalg.eigvals(M).imag); mup = mu[mu > 1e-12]
    tau = np.arctanh(np.minimum(mup, 1 - 1e-16))/np.pi
    tau = np.sort(tau)
    # counting function N(tau) for tau <= 1.0, fit slope
    sel = tau[tau < 1.0]
    k = np.arange(1, len(sel) + 1)
    slope = np.polyfit(sel, k, 1)[0] if len(sel) > 3 else float('nan')
    sp = np.diff(tau[:8])
    print(f"  K={K:4}: max mu={mup.max():.6f}, #tau<1: {len(sel):3}, slope dN/dtau={slope:7.3f}, (4/pi)log K={4/np.pi*np.log(K):7.3f}, first spacings {np.round(sp[:5],4)}")
# compare shapes: box-xp N ~ c*tau (linear) vs RvM N(T) = (T/2pi)log(T/2pi e)+7/8 (superlinear)
M = Tmat(641); mu = np.sort(np.linalg.eigvals(M).imag); tau = np.sort(np.arctanh(np.minimum(mu[mu>1e-12],1-1e-16))/np.pi)
q = [np.searchsorted(tau, t) for t in [0.25, 0.5, 1.0, 2.0, 4.0]]
print("  K=641 counts N(tau) at tau=0.25,0.5,1,2,4:", q, " (linear => ratios 1:2:4:8:16 until saturation at tau~(1/pi)artanh(1-eps))")
# (c) discrete Kutta mode
print("(c) discrete Kutta mode (zero eigenvector, K odd) vs projection of 1/sqrt(1-x^2)")
for K in [41, 161, 641]:
    M = Tmat(K)
    w, V = np.linalg.eig(M)
    i0 = np.argmin(abs(w)); v = np.real_if_close(V[:, i0]); v = np.real(v)/np.linalg.norm(np.real(v))
    # projection coefficients c_n = int p_n(x)/sqrt(1-x^2) dx via Gauss-Chebyshev (exact)
    N = 2*K + 10; th = (np.arange(N) + 0.5)*np.pi/N; x = np.cos(th)
    P = np.polynomial.legendre.legvander(x, K - 1)*np.sqrt((2*np.arange(K) + 1)/2)
    c = (np.pi/N)*P.sum(axis=0); nc = np.linalg.norm(c); c /= nc
    print(f"  K={K:4}: |mu_0|={abs(w[i0]):.2e}, overlap |<v,c>|={abs(v@c):.12f}, ||proj||={nc:.4f} (grows like sqrt(log K): {np.sqrt(np.log(K)):.4f})")
