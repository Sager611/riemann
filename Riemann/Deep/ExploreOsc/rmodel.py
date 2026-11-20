# rmodel.py (explore60): same spectrum, different phases. Riemann-Siegel main sums on [T, T+L]:
#   true: 2 sum_{n<=N} n^{-1/2} cos(theta(t) - t log n)              (zeta; DH analogue with its coefficients a_n)
#   ind : the same with independent uniform phases phi_n added       (spectral data only: Rice/Kac density)
#   mult: phases phi_n = sum_p v_p(n) phi_p, phi_p uniform            (random Euler product, Bohr-Jessen model)
# Output: sign changes / ((theta(T+L)-theta(T))/pi), and the Rice prediction sqrt(lam2/lam0)/theta'.
import sys, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreOsc')
import oscz as O
T = float(sys.argv[1]); L = float(sys.argv[2]); seeds = int(sys.argv[3]); h = 0.01
t = np.arange(T, T + L, h)


def spf(n):
    f, p = {}, 2
    while p*p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f


def rs(nm, coef, phi):
    th = O.theta(nm, t)
    n = np.arange(1, len(coef) + 1)
    out = np.zeros(len(t))
    for i in range(0, len(t), 1000):
        tt = t[i:i + 1000]
        ph = th[i:i + 1000, None] - tt[:, None]*np.log(n)[None, :] + phi[None, :]
        out[i:i + 1000] = 2*(np.abs(coef)[None, :]*n[None, :]**-0.5*np.cos(ph + np.angle(coef)[None, :])).sum(1)
    return out


def count(z):
    return int((np.sign(z[:-1])*np.sign(z[1:]) < 0).sum())


for nm in ('zeta', 'dh'):
    q = O.FUN[nm][0]
    N = int(np.sqrt(q*T/(2*np.pi)))
    a = np.array([1.0 if nm == 'zeta' else {0: 0.0, 1: 1.0, 2: O.KAP, 3: -O.KAP, 4: -1.0}[k % 5] for k in range(1, N + 1)])
    gc = (O.theta(nm, np.array([T + L]))[0] - O.theta(nm, np.array([T]))[0])/np.pi
    Om = (O.theta(nm, np.array([T + 1e-3]))[0] - O.theta(nm, np.array([T - 1e-3]))[0])/2e-3
    w = Om - np.log(np.arange(1, N + 1)); wt = a**2/np.arange(1, N + 1)
    rice = np.sqrt((wt*w**2).sum()/wt.sum())/Om
    res = {'true': [count(rs(nm, a, np.zeros(N)))/gc]}
    fac = [spf(k) for k in range(1, N + 1)]
    for s in range(seeds):
        rng = np.random.default_rng(1000 + s)
        res.setdefault('ind', []).append(count(rs(nm, a, rng.uniform(0, 2*np.pi, N)))/gc)
        php = {p: rng.uniform(0, 2*np.pi) for p in range(2, N + 1) if len(spf(p)) == 1 and spf(p).get(p) == 1}
        phm = np.array([sum(e*php[p] for p, e in f.items()) for f in fac])
        res.setdefault('mult', []).append(count(rs(nm, a, phm))/gc)
    print('%-4s T=%g L=%g N=%d theta-count=%.1f Rice=%.4f | ' % (nm, T, L, N, gc, rice) +
          '  '.join('%s %.4f +- %.4f (n=%d)' % (k, np.mean(v), np.std(v), len(v)) for k, v in res.items()), flush=True)
