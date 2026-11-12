# t2: at-scale falsification check of dpcbeyond's Lemma Z (PROOF.md) on actual zeros, extending R454 (T = 2000, 6000):
#   |R_T(t)| <= C (1 + On(t))^k  (RH range: W = On),  On(t) = sum_rho min(1, (L|gamma - t|)^-A),  k = floor(theta/eta)
# DPC rough coefficients c_n = mu(n) n^-1/2 f(log n/L) prod_{q|n} rho_eta(log q/L), rho_eta(u) = psi(u) omega(u/eta).
import sys, re, numpy as np
sys.path.insert(0, '.')
from t1 import rho, primes_upto
def S(x):
    x = np.clip(x, 1e-12, 1-1e-12); a = np.exp(-1/x); b = np.exp(-1/(1-x)); return a/(a+b)
def fsm(x, theta):  # smooth, f(0) = 1, f = 0 on [theta, inf)
    x = np.asarray(x, float); out = np.zeros_like(x); m = x < theta; out[m] = 1 - S(np.clip(x[m]/theta, 0, 1)); return out
def coeffs(T, theta, eta, b):
    L = np.log(T); N = int(np.exp(theta*L)); pmin = np.exp(eta*L)
    ps = primes_upto(N); ps = ps[ps > pmin].astype(np.int64)
    k = int(np.floor(theta/eta))
    ns, cs = [np.array([1], np.int64)], [np.array([1.0])]
    # products of j distinct primes, j <= k (k <= 2 or 3 here), n <= N
    cur = [(np.array([1], np.int64), np.array([1.0]), np.array([0], np.int64))]  # (n, weight, index of last prime)
    for j in range(1, k+1):
        nxt = []
        for n0, w0, i0 in cur:
            for n_, w_, i_ in zip(n0, w0, i0):
                lim = N // n_
                qs = ps[(ps <= lim)]
                if j > 1: qs = qs[qs > ps[i_]] if len(qs) else qs
                if len(qs) == 0: continue
                idx = np.searchsorted(ps, qs)
                nxt.append((n_*qs, w_*(-1)*rho(np.log(qs)/L, eta, theta, b), idx))
        if not nxt: break
        n1 = np.concatenate([a for a, _, _ in nxt]); w1 = np.concatenate([w for _, w, _ in nxt]); i1 = np.concatenate([i for _, _, i in nxt])
        ns.append(n1); cs.append(w1); cur = [(n1, w1, i1)]
    n = np.concatenate(ns).astype(float); w = np.concatenate(cs)
    c = w/np.sqrt(n)*fsm(np.log(n)/L, theta)
    keep = c != 0
    return n[keep], c[keep], k, L
def RT(ts, n, c):
    out = np.empty(len(ts), complex); ln = np.log(n); ch = max(1, int(1.5e6/len(n)))
    for i in range(0, len(ts), ch):
        out[i:i+ch] = (np.exp(-1j*np.outer(ts[i:i+ch], ln))*c).sum(axis=1)
    return out
def On(ts, zeros, L, A):
    out = np.zeros(len(ts))
    for i, t in enumerate(ts):
        g = zeros[np.abs(zeros - t) < 50/L]
        x = L*np.abs(g - t); out[i] = np.sum(np.minimum(1.0, np.where(x > 0, x, 1e-300)**(-A)))
    return out
def gram_windows(path):
    W = []; tc = None
    for line in open(path):
        if line.startswith('ZG'): tc = float(re.search(r'tc=([0-9.eE+-]+)', line).group(1))
        elif line.startswith('ZZ'):
            W.append((tc, tc + np.array([float(v) for v in line.split()[2:]])))
    return W
def report(tag, ratio, extra=''):
    print('%s: points %d  max %.3f  99.9%% %.3f  99%% %.3f  %s' % (tag, len(ratio), ratio.max(), np.quantile(ratio, 0.999), np.quantile(ratio, 0.99), extra), flush=True)
if __name__ == '__main__':
    theta = 0.55; eta = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2; b = theta + 0.1
    print('theta %.2f eta %.2f' % (theta, eta))
    import os
    R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..') + '/'
    z = np.load(R + 'Riemann/Deep/ExploreHelio/results/zeta_zeros.npy')
    for Tlab, T, ts, zs in [('T=1e4 (t in [1e4,1.99e4], step 0.02)', 1.0e4, np.arange(1.0e4, 1.99e4, 0.02), z)]:
        n, c, k, L = coeffs(T, theta, eta, b)
        for A in (2, 4):
            o = On(ts, zs, L, A); r = np.abs(RT(ts, n, c))
            report('%s terms=%d k=%d A=%d' % (Tlab, len(n), k, A), r/(1 + o)**k, 'max|R_T| %.3f' % r.max())
    for path, lab in [('Riemann/Deep/ExploreDPP/logs/zeros_1e8.txt', 't~1e8'), ('Riemann/Deep/ExploreDPP/logs/zeros_1e11.txt', 't~9.8e10')]:
        W = gram_windows(R + path)
        T = W[0][0]/1.5
        n, c, k, L = coeffs(T, theta, eta, b)
        rr = {2: [], 4: []}; mx = 0
        for tc, zz in W:
            ts = tc + np.arange(-3, 3, 0.05)
            r = np.abs(RT(ts, n, c)); mx = max(mx, r.max())
            for A in (2, 4):
                rr[A].append(r/(1 + On(ts, zz, L, A))**k)
        for A in (2, 4):
            report('%s T=%.3g windows=%d terms=%d k=%d A=%d' % (lab, T, len(W), len(n), k, A), np.concatenate(rr[A]), 'max|R_T| %.3f' % mx)
