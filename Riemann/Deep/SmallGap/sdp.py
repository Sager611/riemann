# SDP for Lamzouri 2611.14486 Thm Main2 extremal problem (float; not a certificate)
# maximize J(f) = f(0) + 2 int_0^theta f(a) a da over f(x) = g(x/s), g = sum_{k<=d} a_k psi_{2k}
# psi_n: L2-normalized Hermite functions, psi_n(y) = (2pi)^{1/4} h_n(sqrt(2pi) y), psi_n^ = (-i)^n psi_n.
# constraints (exact polynomial identities, matched in the orthonormal psi basis):
#   g(y)      = sum_{i,j<=d, i=j mod 2} X_ij phi_i phi_j(y)                  (f >= 0)
#   -g^(t)/1  = sum X'_ij phi_i phi_j(t) + (t^2 - c) sum Y'_ij phi_i phi_j(t), c = (s lam)^2   (fhat <= 0 on |xi| >= lam)
#   phi_i(y) = h_i(sqrt(pi) y)  (so phi_i phi_j = poly * e^{-pi y^2});  fhat(0) = s * g^(0) = 1.
# usage: OMP_NUM_THREADS=1 uv run --no-project --with cvxpy --with numpy --with scipy python3 sdp.py d lam [s]
import numpy as np, cvxpy as cp, sys, time
from scipy.special import roots_hermite
from numpy.polynomial.legendre import leggauss
THETA = 4/7 - 1e-5

def hfun(u, n):
    """[len(u), n+1]: orthonormal Hermite functions h_m(u) = (2^m m! sqrt pi)^{-1/2} H_m(u) e^{-u^2/2}."""
    u = np.atleast_1d(np.asarray(u, float)); out = np.empty((len(u), n+1))
    out[:, 0] = np.pi**-0.25*np.exp(-u*u/2)
    if n >= 1: out[:, 1] = np.sqrt(2.0)*u*out[:, 0]
    for m in range(1, n):
        out[:, m+1] = np.sqrt(2.0/(m+1))*u*out[:, m] - np.sqrt(m/(m+1))*out[:, m-1]
    return out

def psi_even(y, d):
    return (2*np.pi)**0.25*hfun(np.sqrt(2*np.pi)*np.atleast_1d(np.asarray(y, float)), 2*d)[:, 0::2]

_cache = {}
def tensors(d):
    """T[k,i,j] = <phi_i phi_j, psi_2k>, T2[k,i,j] = <t^2 phi_i phi_j, psi_2k> (exact Gauss-Hermite)."""
    if d in _cache: return _cache[d]
    N = 2*d + 4
    v, w = roots_hermite(N)               # weight e^{-v^2}; y = v / sqrt(2 pi)
    y = v/np.sqrt(2*np.pi)
    W = w*np.exp(v*v)/np.sqrt(2*np.pi)     # integral F(y) dy = sum W F(y_m)
    Ph = hfun(np.sqrt(np.pi)*y, d)         # phi_i(y_m)
    Ps = psi_even(y, d)                    # psi_2k(y_m)
    T = np.empty((d+1, d+1, d+1)); T2 = np.empty((d+1, d+1, d+1))
    for k in range(d+1):            # loop: avoids an O(N d^3) einsum intermediate (memory)
        T[k] = (Ph*(W*Ps[:, k])[:, None]).T @ Ph
        T2[k] = (Ph*(W*y*y*Ps[:, k])[:, None]).T @ Ph
    _cache[d] = (T, T2)
    return T, T2

def blocks(n):
    return [list(range(0, n+1, 2)), list(range(1, n+1, 2))]

def gram_expr(T, n, margin):
    """sum over parity blocks of <G_b, T[:, b, b]> ; returns (expr vector over k, list of PSD vars)."""
    exprs = []; vs = []
    for b in blocks(n):
        if not b: continue
        G = cp.Variable((len(b), len(b)), PSD=True); vs.append(G)
        Tb = T[:, b][:, :, b].reshape(T.shape[0], -1)
        exprs.append(Tb @ cp.vec(G, order='C'))
        if margin > 0: vs[-1] = G
    return sum(exprs), vs

def solve(d, s, lam, theta=THETA, solver='CLARABEL', verbose=False, margin=0.0, full=False):
    T, T2 = tensors(d)
    c = (s*lam)**2
    sg = np.array([(-1.0)**k for k in range(d+1)])
    a = cp.Variable(d+1)
    e1, G1 = gram_expr(T, d, margin)
    e2, G2 = gram_expr(T[:, :d, :d], d-1, margin)
    e3, G3 = gram_expr(T2[:, :d, :d] - c*T[:, :d, :d], d-1, margin)
    e2b, G2b = gram_expr(T, d, margin)
    cons = [a == e1, -cp.multiply(sg, a) == e2b + e3]
    if margin > 0:
        cons += [G - margin*np.eye(G.shape[0]) >> 0 for G in G1 + G2b + G3]
    p00 = psi_even(0.0, d)[0]
    cons.append(s*(sg*p00) @ a == 1)
    gx, gw = leggauss(400); yy = (gx+1)/2*theta/s; ww = gw*theta/s/2
    cobj = p00 + 2*s*s*(ww*yy) @ psi_even(yy, d)
    prob = cp.Problem(cp.Maximize(cobj @ a), cons)
    try:
        if solver == 'CLARABEL':
            prob.solve(solver=solver, verbose=verbose, tol_gap_abs=1e-11, tol_gap_rel=1e-11, tol_feas=1e-11, tol_ktratio=1e-9, max_iter=400)
        else:
            prob.solve(solver=solver, verbose=verbose)
    except Exception as ex:
        return float('nan'), None, 'FAIL'
    if prob.status not in ('optimal', 'optimal_inaccurate'):
        return float('nan'), None, prob.status
    return prob.value, a.value, prob.status

def check(d, s, lam, a, theta=THETA):
    """float sanity: min f on [0,15], max fhat on [lam,15] (30001-point grids, chunked for memory), fhat(0), J."""
    sg = np.array([(-1.0)**k for k in range(d+1)])
    x = np.linspace(0, 15, 30001); xi = np.linspace(lam, 15, 30001)
    f = np.concatenate([psi_even(c/s, d) @ a for c in np.array_split(x, 30)])
    fh = s*np.concatenate([psi_even(s*c, d) @ (sg*a) for c in np.array_split(xi, 30)])
    gx, gw = leggauss(400); yy = (gx+1)/2*theta; ww = gw*theta/2
    J = (psi_even(0.0, d)[0] @ a) + 2*(ww*yy) @ (psi_even(yy/s, d) @ a)
    fh0 = s*(psi_even(0.0, d)[0] @ (sg*a))
    return dict(minf=f.min(), maxfh=fh.max(), fh0=fh0, J=J)

def best_s(d, lam, lo=0.6, hi=1.5, it=12, **kw):
    g = (5**0.5 - 1)/2
    def F(s):
        r = solve(d, s, lam, **kw); return r[0] if r[2] == 'optimal' else -1e9
    a, b = lo, hi
    c, e = b - g*(b-a), a + g*(b-a); fc, fe = F(c), F(e)
    for _ in range(it):
        if fc > fe: b, e, fe = e, c, fc; c = b - g*(b-a); fc = F(c)
        else: a, c, fc = c, e, fe; e = a + g*(b-a); fe = F(e)
    return (c, fc) if fc > fe else (e, fe)

if __name__ == '__main__':
    d = int(sys.argv[1]); lam = float(sys.argv[2])
    t = time.time()
    if len(sys.argv) > 3:
        for s in map(float, sys.argv[3].split(',')):
            v, a, st = solve(d, s, lam)
            print(f"d={d} s={s} lam={lam}: J={v!r} {st} {check(d, s, lam, a) if a is not None else ''} ({time.time()-t:.1f}s)", flush=True)
    else:
        s, v = best_s(d, lam)
        print(f"d={d} lam={lam}: best s={s:.4f} J={v!r} ({time.time()-t:.1f}s)")
