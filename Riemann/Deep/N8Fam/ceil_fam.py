# Trial-state ceilings (as in 2610.08965 Prop. 7.3 / N8Opt ceil_range.py) for family windows at bandwidth lam.
# Units: z (mean zero spacing 1). A periodic state = gap sequence g (period p, length P), pair energy e per period over
# ordered pairs at index distance <= s.  Mixtures: sum w P = 1, sum w (e - n) <= -(2 - C); ceiling = min sum w n.
# Any certificate that only sees pairs at index distance <= s proves at most this (upper bound; floating point).
import numpy as np, sys, itertools
from scipy.optimize import linprog
sys.path.insert(0, '.')
from pbeta import pbeta
from cfam import C_model, cos_window, A2

def kernel_from_cells(f, e):
    # even step density, value f_i on +-[e_i, e_{i+1}]
    def K(t):
        t = np.asarray(t, float)
        out = np.zeros_like(t)
        small = np.abs(t) < 1e-12
        tt = np.where(small, 1.0, t)
        for fi, a, b in zip(f, e[:-1], e[1:]):
            out += fi * (np.sin(2 * np.pi * tt * b) - np.sin(2 * np.pi * tt * a)) / (np.pi * tt)
        out[small] = np.sum(f * 2 * np.diff(e))
        return out
    return K

def paper_kernel(lam):
    a = A2
    def K(t):
        t = lam * np.asarray(t, float)
        s_ = np.sinc(t)
        for j, aj in enumerate(a, 1):
            s_ = s_ + aj / 2 * (np.sinc(t - j) + np.sinc(t + j))
        return s_
    return K

def make_K2(K, tmax=200.0, H=1e-4):
    tg = np.arange(0, tmax, H); kg = K(tg) ** 2
    return lambda t: np.interp(np.abs(t), tg, kg, right=0.0)

def state_stats(G, K2, s):
    # G: (m, p) gap sequences; returns e (pair energy per period, ordered pairs, index distance <= s), n = p, P
    m, p = G.shape
    reps = s // p + 2
    Y = np.concatenate([np.zeros((m, 1)), np.cumsum(np.tile(G, reps), axis=1)], axis=1)
    e = np.zeros(m)
    for k in range(1, s + 1):
        e += K2(Y[:, k:k + p] - Y[:, :p]).sum(axis=1)
    return 2 * e, np.full(m, p, float), G.sum(axis=1)

def ceiling(E, Nn, P, C):
    r = linprog(Nn, A_ub=[E - Nn], b_ub=[-(2 - C)], A_eq=[P], b_eq=[1.0], bounds=[(0, None)] * len(E), method='highs')
    return (r.fun, r.x) if r.status == 0 else (np.inf, None)

def library(K2, s, lo, hi, step, rng, nrand=4000):
    out = []
    g1 = np.arange(lo, hi, step)
    out.append(state_stats(g1[:, None], K2, s) + (g1[:, None],))
    g2 = np.arange(lo, hi, 2 * step)
    a, b = np.meshgrid(g2, g2); m = a <= b
    G2 = np.stack([a[m], b[m]], axis=1)
    out.append(state_stats(G2, K2, s) + (G2,))
    for p in (3, 4, 5):
        G = rng.uniform(lo, hi, size=(nrand, p))
        out.append(state_stats(G, K2, s) + (G,))
    return out

def refine(G, K2, s, C, others, iters=60):
    # coordinate pattern search on one state's gaps, keeping the others fixed, to lower the LP value
    best_G = G.copy()
    def val(Gc):
        e, n, P = state_stats(Gc[None, :], K2, s)
        E = np.concatenate([others[0], e]); Nn = np.concatenate([others[1], n]); PP = np.concatenate([others[2], P])
        return ceiling(E, Nn, PP, C)[0]
    best = val(best_G); step = 0.02
    while step > 1e-5:
        moved = False
        for i in range(len(best_G)):
            for sg in (1, -1):
                tr = best_G.copy(); tr[i] += sg * step
                if tr[i] < 0.05: continue
                v = val(tr)
                if v < best - 1e-12:
                    best, best_G, moved = v, tr, True
        if not moved: step /= 2
    return best, best_G

def run(name, K, C, s_list=(7, 40), lo=0.30, hi=3.2, step=0.01, seed=1):
    rng = np.random.default_rng(seed)
    K2 = make_K2(K)
    print('== %s: plain 2-C = %.6f' % (name, 2 - C), flush=True)
    for s in s_list:
        lib = library(K2, s, lo, hi, step, rng)
        E = np.concatenate([l[0] for l in lib]); Nn = np.concatenate([l[1] for l in lib]); P = np.concatenate([l[2] for l in lib])
        Gs = [g for l in lib for g in l[3]]
        val, w = ceiling(E, Nn, P, C)
        idx = np.nonzero(w > 1e-9)[0]
        # refine the active states one at a time, against the other active states only (LP optimum uses <= 2 states)
        cur = val
        for it in range(2):
            for j in idx:
                mask = np.zeros(len(E), bool); mask[idx] = True; mask[j] = False
                v, Gn = refine(np.array(Gs[j]), K2, s, C, (E[mask], Nn[mask], P[mask]))
                if v < cur - 1e-12:
                    e1, n1, P1 = state_stats(Gn[None, :], K2, s)
                    E[j], Nn[j], P[j] = e1[0], n1[0], P1[0]; Gs[j] = Gn; cur = v
            val, w = ceiling(E, Nn, P, C); idx = np.nonzero(w > 1e-9)[0]
        print('  s=%d: ceiling <= %.6f (gain over plain <= %.6f); active states %s' %
              (s, val, val - (2 - C), [np.round(Gs[j], 3).tolist() for j in idx]), flush=True)

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('zeta', 'all'):   # sanity: paper window at lam = 1 should give <= 0.6735660 (s=7) and <= 0.6738216 (all)
        run('paper window, lam=1 (zeta check)', paper_kernel(1.0), C_model(cos_window(A2), 1.0, 'zeta'))
    if which in ('fam', 'all'):
        pv, f, e = pbeta(2.0, M=200)
        Cs = 2 - pv
        run('Johns-optimal window (F1, beta=2)', kernel_from_cells(f, e), Cs)
        lam = 1.999
        run('paper window dilated, lam=%.3f' % lam, paper_kernel(lam), C_model(cos_window(A2), lam, 'gue'))
    if which in ('mid', 'all'):
        for beta in (1.5, 4 / 3):
            pv, f, e = pbeta(beta, M=200)
            run('Johns-optimal window (F1, beta=%.4f)' % beta, kernel_from_cells(f, e), 2 - pv)
