
# Quick test: complex Q (Bernard's class Q = 1 + odd real + i*even real in z=1-2y) with complex P. Untrusted numerics.
import numpy as np, sys
from numpy.polynomial import chebyshev as C
from numpy.polynomial import legendre as Lg
NQ = 600
xg, wg = Lg.leggauss(NQ); yg = (xg+1)/2; wg = wg/2
def Qbasis(R, M):
    z = 1-2*yg; E = np.exp(R*yg)
    Phi=[E]; dPhi=[R*E]; kinds=['one']
    for m in range(1, 2*M+1):
        c = np.zeros(m+1); c[m] = 1.0
        T = C.chebval(z, c); dT = C.chebval(z, C.chebder(c))*(-2.0)
        Phi.append(E*(T-1)); dPhi.append(R*E*(T-1)+E*dT); kinds.append('odd' if m%2 else 'even')
    return np.array(Phi), np.array(dPhi), kinds
def Pbasis(K):
    # P = x + sum_k p_k (x^k - x), k=2..K, on nodes
    X = yg
    B = [X] + [X**k - X for k in range(2, K+1)]
    dB = [np.ones_like(X)] + [k*X**(k-1) - 1 for k in range(2, K+1)]
    return np.array(B), np.array(dB)
def cval(v, p, Phi, dPhi, PB, dPB, xi):
    w = v@Phi; dw = v@dPhi; P = p@PB; dP = p@dPB
    A = np.sum(wg*abs(w)**2); Cc = np.sum(wg*abs(dw)**2); B = np.sum(wg*w*np.conj(dw))
    aP = np.sum(wg*abs(dP)**2); bP = np.sum(wg*abs(P)**2); D = np.sum(wg*dP*np.conj(P))
    return 1 + (A*aP + xi**2*Cc*bP + 2*xi*np.real(B*D))/xi, A, Cc, B
def minimize_quadratic(f, n):
    # f: real vector u -> real quadratic; recover by finite differences (exact for quadratics)
    f0 = f(np.zeros(n)); g = np.zeros(n); H = np.zeros((n,n)); I = np.eye(n)
    fp = np.array([f(I[i]) for i in range(n)]); fm = np.array([f(-I[i]) for i in range(n)])
    g = (fp - fm)/2; Hd = fp + fm - 2*f0
    for i in range(n):
        H[i,i] = Hd[i]
        for j in range(i):
            H[i,j] = H[j,i] = (f(I[i]+I[j]) - fp[i] - fp[j] + f0)/2*2/2
    # f(u) = f0 + g.u + u^T (H/2) u ; note H[i,j] computed as f(ei+ej)-f(ei)-f(ej)+f0 = H_ij (off-diag of Hessian)
    for i in range(n):
        for j in range(i):
            H[i,j] = H[j,i] = f(I[i]+I[j]) - fp[i] - fp[j] + f0
    u = np.linalg.solve(H, -g)
    return u, f(u)
def run(R, xi, M, K, cplx, iters=60):
    Phi, dPhi, kinds = Qbasis(R, M); PB, dPB = Pbasis(K)
    nq = len(kinds)-1
    qidx = [i for i in range(1, nq+1) if kinds[i]=='odd' or cplx]
    def v_of(u):
        v = np.zeros(nq+1, complex); v[0] = 1
        for t, i in enumerate(qidx): v[i] = u[t] if kinds[i]=='odd' else 1j*u[t]
        return v
    np_ = K-1
    def p_of(u):
        p = np.zeros(K, complex); p[0] = 1
        if cplx: p[1:] = u[:np_] + 1j*u[np_:]
        else: p[1:] = u[:np_]
        return p
    # init: P = sinh(0.77x)/sinh(0.77) approx via fit
    p = p_of(np.zeros(2*np_ if cplx else np_)); uq = np.zeros(len(qidx))
    lp = np.linalg.lstsq(PB[1:].T, np.sinh(0.77*yg)/np.sinh(0.77) - yg, rcond=None)[0]
    p[1:] = lp
    if cplx: p[1:] = p[1:] + 1j*np.random.default_rng(int(sys.argv[5]) if len(sys.argv)>5 else 1).normal(0, 0.5, np_)
    for it in range(iters):
        uq, cq = minimize_quadratic(lambda u: cval(v_of(u), p, Phi, dPhi, PB, dPB, xi)[0], len(qidx)); v = v_of(uq)
        nu = 2*np_ if cplx else np_
        up, cp = minimize_quadratic(lambda u: cval(v, p_of(u), Phi, dPhi, PB, dPB, xi)[0], nu); p = p_of(up)
    c, A, Cc, B = cval(v, p, Phi, dPhi, PB, dPB, xi)
    return 1 - np.log(c)/R, c, B
xi = eval(sys.argv[1]); R = float(sys.argv[2]); M = int(sys.argv[3]); K = int(sys.argv[4])
for cplx in (False, True):
    k, c, B = run(R, xi, M, K, cplx)
    print(f"xi={xi:.6f} R={R} M={M} K={K} complex={cplx}: kappa={k:.9f} c={c:.6f} ImB={B.imag:.3e}", flush=True)
