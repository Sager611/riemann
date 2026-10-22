#!/usr/bin/env python3
"""Numerics for Riemann/Deep/Smale (arXiv 2610.21472).  Needs numpy, scipy.  Each mode < 60 s, < 200 MB.
  dieudonne        min |root| of the Dieudonne polynomials of H (pair 9/16, 3/8) over c in [-1,1]; zeros of H'
  cshape           degree-one pair W = e^{3.1z}: robust zero count of ((W-W(b))B_a)' in D, and 1/B(0)
  nmin             potential-theory count n_min = max(-sigma/omega) on Gamma_t = W({|B_a| = t}) (paper pair, C-shape)
  slp D TRIALS     best S(P) over degree D from random starts (P' = prod(1-w/c_j), SLP on min log|P(c_j)/c_j|)
  sym N K TRIALS   k-symmetric family P = w R(w^K), P' = T(w^K) (K = 0: limit G = u e^L, objective min Re L)
"""
import sys, numpy as np
from scipy.optimize import linprog
from scipy.special import roots_jacobi
Pn = np.polynomial.polynomial

def paper_pair():
    x, y = 9/16, 3/8; u = complex(-x, -y); a, b = 1 + u, 1 - 1j*u
    Q = Pn.polypow([-1, 1], 4).astype(complex); Q[0] -= u**4
    quo, _ = Pn.polydiv(Q, [-a, 1]); H = Pn.polymul([1, -np.conj(a)], quo)
    return Q, H, a, b

def dieudonne():
    Q, H, a, b = paper_pair(); worst = (9, None)
    for cc in np.linspace(-1, 1, 20001):
        U = [1, 2*cc, 4*cc*cc - 1, 8*cc**3 - 4*cc]
        co = [H[k + 1]*U[k] for k in range(4)]
        while abs(co[-1]) < 1e-13: co = co[:-1]
        m = np.abs(np.roots(co[::-1])).min()
        worst = min(worst, (m, cc))
    print('min |root| of P_c over c in [-1,1]: %.6f at c = %.4f' % worst)
    print('|zeros of H\'|:', np.round(np.sort(np.abs(Pn.polyroots(Pn.polyder(H)))), 5))
    eta = (b - a)/(1 - np.conj(a)*b); print('1/|eta| = %.10f' % (1/abs(eta)))

def winding(F, M=4096, depth=22):
    t = 2*np.pi*np.arange(M + 1)/M; v = F(np.exp(1j*t)); tot = 0.0; mn = np.abs(v).min()
    st = [(t[i], t[i + 1], v[i], v[i + 1], 0) for i in range(M)]
    while st:
        s0, s1, v0, v1, d = st.pop(); dd = np.angle(v1/v0)
        if abs(dd) < 0.15: tot += dd; continue
        if d >= depth: return None, mn
        sm = (s0 + s1)/2; vm = F(np.exp(1j*sm)); mn = min(mn, abs(vm))
        st += [(s0, sm, v0, vm, d + 1), (sm, s1, vm, v1, d + 1)]
    return round(tot/(2*np.pi), 6), mn

CS = dict(alpha=3.1, b=-0.1814 + 0.6221j, a=-0.0555 - 0.8659j)
def cshape():
    al, b, a = CS['alpha'], CS['b'], CS['a']
    W = lambda s: np.exp(al*s); dW = lambda s: al*np.exp(al*s)
    N = lambda s: dW(s)*(s - a)*(1 - np.conj(a)*s) + (W(s) - W(b))*(1 - abs(a)**2)
    w, mn = winding(N)
    print('W = e^{%.1f z}: zeros of ((W-W(b))B_a)\' in D = %s (min |N| on T = %.2e); 1/B(0) = %.5f'
          % (al, w, mn, abs(1 - np.conj(a)*b)/abs(b - a)))

def kress(gam, dg):
    M = len(gam); k = np.arange(M); tk = 2*np.pi*k/M; m = np.arange(1, M//2)
    Rk = -(4*np.pi/M)*(np.cos(np.outer(tk, m)) @ (1/m)) - (4*np.pi/M**2)*np.cos(M*tk/2)
    Iq = np.arange(M); D = (Iq[:, None] - Iq[None, :]) % M; Tm = tk[D]
    with np.errstate(divide='ignore', invalid='ignore'):
        L = np.log(np.abs(gam[:, None] - gam[None, :])) - 0.5*np.log(4*np.sin(Tm/2)**2)
    L[Iq, Iq] = np.log(np.abs(dg))
    A = 0.5*Rk[D] + (2*np.pi/M)*L
    K = np.zeros((M + 1, M + 1)); K[:M, :M] = A; K[:M, M] = -1; K[M, :M] = 2*np.pi/M
    return K

def nmin_curve(Wf, dWf, a, b, t, M=800):
    """Gamma = W(circle |B_a| = t) in w = (W - W(b))/N coordinates; model f = w B_a(z)/B_a(b)."""
    N = (1 - abs(b)**2)*dWf(b); eta = (b - a)/(1 - np.conj(a)*b)
    ct = a*(1 - t*t)/(1 - t*t*abs(a)**2); rt = t*(1 - abs(a)**2)/(1 - t*t*abs(a)**2)
    s = 2*np.pi*np.arange(M)/M; zz = ct + rt*np.exp(1j*s)
    gam = (Wf(zz) - Wf(b))/N; dW = dWf(zz); dg = dW/N*1j*rt*np.exp(1j*s)
    Ba = (zz - a)/(1 - np.conj(a)*zz); dBa = (1 - abs(a)**2)/(1 - np.conj(a)*zz)**2
    fp = (dW/N*Ba + gam*dBa)/(eta*dW/N)
    K = kress(gam, dg)
    rhs = np.zeros(M + 1); rhs[M] = 1; e = np.linalg.solve(K, rhs)[:M]
    rhs = np.zeros(M + 1); rhs[:M] = np.log(np.abs(fp)); sg = np.linalg.solve(K, rhs)[:M]
    return t/abs(eta), np.max(-sg/e), e.min()*2*np.pi

def nmin():
    Q, H, a, b = paper_pair(); dH = Pn.polyder(H)
    Wf = lambda zz: Pn.polyval(zz, H); dWf = lambda zz: Pn.polyval(zz, dH)
    eta = abs((b - a)/(1 - np.conj(a)*b))
    for At in (1.01, 1.02):
        print('paper pair  A_Gamma=%.4f  n_min=%.3g  min eq.density=%.2g' % nmin_curve(Wf, dWf, a, b, At*eta))
    al, b, a = CS['alpha'], CS['b'], CS['a']
    Wf = lambda zz: np.exp(al*zz); dWf = lambda zz: al*np.exp(al*zz)
    eta = abs((b - a)/(1 - np.conj(a)*b))
    for At in (1.003, 1.01):
        print('C-shape     A_Gamma=%.4f  n_min=%.3g  min eq.density=%.2g' % nmin_curve(Wf, dWf, a, b, At*eta))

def slp_maxmin(fun, x, delta=0.1, iters=600, dmin=1e-12):
    phi, J = fun(x); best = phi.min(); n = len(x)
    for _ in range(iters):
        c = np.zeros(n + 1); c[-1] = -1
        res = linprog(c, A_ub=np.hstack([-J, np.ones((len(phi), 1))]), b_ub=phi,
                      bounds=[(-delta, delta)]*n + [(None, None)], method='highs')
        if res.status != 0:
            delta *= 0.5
            if delta < dmin: break
            continue
        d = res.x[:n]; pred = res.x[-1] - best
        phin, Jn = fun(x + d); nb = phin.min()
        if not np.isfinite(nb) or nb - best <= 0.1*max(pred, 0) or pred <= 1e-12:
            delta *= 0.5
            if delta < dmin: break
            continue
        x, phi, J, best = x + d, phin, Jn, nb
        if nb - best > 0.75*pred: delta = min(2*delta, 10)
    return x, best

def make_obj(N, k):
    """Critical points c_1..c_N of P' (k=1), of T for P = w R(w^k) (k>1), or limit (k=0)."""
    nq = N//2 + 4
    if k == 0:
        xg, wg = np.polynomial.legendre.leggauss(nq); tau = (xg + 1)/2; wt = wg/2
    else:
        beta = 1.0/k - 1.0; xg, wg = roots_jacobi(nq, 0.0, beta); tau = (1 + xg)/2; wt = wg*0.5**(beta + 1)/k
    idx = np.arange(N)
    def fun(xv):
        c = xv[0::2] + 1j*xv[1::2]
        Fm = 1 - tau[:, None, None]*(c[:, None]/c[None, :])[None]
        Tv = np.prod(Fm, axis=2)
        with np.errstate(divide='ignore', invalid='ignore'): G = Tv[:, :, None]/Fm
        if k == 0:
            val = (wt[:, None]*(Tv - 1)/tau[:, None]).sum(0)
            D = (wt[:, None, None]*G*(c[:, None]/c[None, :]**2)[None]).sum(0)
            D[idx, idx] = -1/c + (1/c)*(wt[:, None]*G[:, idx, idx]).sum(0)
            phi, r = val.real, D
        else:
            val = (wt[:, None]*Tv).sum(0)
            D = (wt[:, None, None]*G*tau[:, None, None]*(c[:, None]/c[None, :]**2)[None]).sum(0)
            Tp = -(G/c[None, None, :]).sum(2)
            D[idx, idx] = (wt[:, None]*(tau[:, None]*Tp + G[:, idx, idx]*tau[:, None]/c[None, :])).sum(0)
            phi, r = np.log(np.abs(val)), D/val[:, None]
        J = np.zeros((N, 2*N)); J[:, 0::2] = r.real; J[:, 1::2] = -r.imag
        return phi, J
    return fun

def search(N, k, trials, seed=1):
    f = make_obj(N, k); rng = np.random.default_rng(seed); vals = []
    for _ in range(trials):
        c = rng.normal(size=N) + 1j*rng.normal(size=N); xv = np.zeros(2*N); xv[0::2] = c.real; xv[1::2] = c.imag
        vals.append(slp_maxmin(f, xv)[1])
    return np.array(sorted(vals))

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'dieudonne': dieudonne()
    elif mode == 'cshape': cshape()
    elif mode == 'nmin': nmin()
    elif mode == 'slp':
        d, tr = int(sys.argv[2]), int(sys.argv[3]); v = np.exp(search(d - 1, 1, tr))
        print('degree %d: best S = %.5f over %d starts; (d-1)/d = %.5f' % (d, v.max(), tr, (d - 1)/d))
    elif mode == 'sym':
        N, k, tr = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); v = search(N, k, tr)
        print('N=%d k=%d: best objective %.5f (%s); Smale extremal gives %.5f'
              % (N, k, v.max(), 'min Re L' if k == 0 else 'log S', -1/N if k == 0 else np.log(k*N/(k*N + 1))))
