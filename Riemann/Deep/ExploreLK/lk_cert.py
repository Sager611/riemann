#!/usr/bin/env python3
"""explore58 certificate (Arb, python-flint): exact cone membership in the period-N model.
1. numpy: CUE samples -> hard-core gaps >= delta+margin -> zeros rounded to rationals p/10^6 (gaps rechecked exactly in Fractions)
2. LP (HiGHS) picks a basic feasible mixture (<= rows columns)
3. Arb (prec 256): spectra |e_j(zeta)|^2 2^{-4n} (acb product), slopes, d; square solve; all weights > 0 certified.
mode 'disc': rows = flat discrete spectrum t_0 = 1/(N+1), t_j = 2/(N+1) (j=1..N/2), d = m1_target.
mode 'cont': rows = m_0..m_J with Z's continuum values pi^{2j}/(2j+1), d = pi^2/3; columns also scaled by s (rational grid).
"""
import sys, math, json
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
import flint
from flint import arb, acb, arb_mat
sys.path.insert(0, 'Riemann/Deep/ExploreLK')
from lk_cue import cue_angles, hardcore, stats
PI = math.pi
flint.ctx.prec = 256

def arb_column(xq, N):
    """xq: list of Fractions (zeros in [0,N)). returns (w[0..n] folded, m0, d) as arb (unnormalised)."""
    n2 = len(xq); n = n2 // 2
    pi = arb.pi()
    poly = [acb(1)]
    for x in xq:
        z = (-2 * pi * arb(x.numerator) / arb(x.denominator) / N)  # angle of zeta_k = exp(-2 pi i x/N)
        zk = acb(z.cos(), z.sin())
        new = [acb(0)] * (len(poly) + 1)
        for i, c in enumerate(poly):
            new[i] += c; new[i + 1] -= c * zk
        poly = new
    scale = arb(2) ** (-4 * n)
    p = [scale * (c.real ** 2 + c.imag ** 2) for c in poly]  # |c_{j-n}|^2, j = 0..2n
    w = [p[n]] + [p[n + j] + p[n - j] for j in range(1, n + 1)]
    m0 = sum(w, arb(0))
    d = arb(0)
    for k in range(n2):
        prod = arb(1)
        for l in range(n2):
            if l != k:
                a = pi * (arb((xq[k] - xq[l]).numerator) / arb((xq[k] - xq[l]).denominator)) / N
                prod *= a.sin()
        d += (pi / N * prod) ** 2
    d = d / N
    return w, m0, d


def Dfloat(q, N, s, alphas):
    xs = np.array([float(x) for x in q]); out = []
    t = np.arange(8 * N) * (N / (8 * N))
    lf = np.sum(np.log(np.abs(np.sin(PI * (t[:, None] - xs[None, :]) / N)) + 1e-300), axis=1); sh = np.max(lf)
    f = np.prod(np.sign(np.sin(PI * (t[:, None] - xs[None, :]) / N)), axis=1) * np.exp(lf - sh)
    m0 = np.mean(f ** 2)
    for a in alphas:
        pts = xs + s * float(a)
        lv = np.sum(np.log(np.abs(np.sin(PI * (pts[:, None] - xs[None, :]) / N)) + 1e-300), axis=1)
        out.append(s / N * np.sum(np.exp(2 * (lv - sh))) / m0)
    return out

def Darb(xq, N, s, a):
    """unnormalised D for g(x) = f(s x): (s/N) sum_k f(x_k + s a)^2, f = prod sin(pi (x - x_l)/N) (same scale as arb_column)."""
    pi = arb.pi(); tot = arb(0)
    sa = s * a
    for k in range(len(xq)):
        prod = arb(1)
        for l in range(len(xq)):
            u = xq[k] + sa - xq[l]
            prod *= (pi * arb(u.numerator) / arb(u.denominator) / N).sin()
        tot += prod ** 2
    return arb(s.numerator) / s.denominator / N * tot

def rational_zeros(xs, N, delta, den=10 ** 6):
    q = sorted(Fr(int(round(x * den)), den) % N for x in xs)
    g = [q[i + 1] - q[i] for i in range(len(q) - 1)] + [q[0] + N - q[-1]]
    return q, min(g)

if __name__ == '__main__':
    mode = sys.argv[1]; N = int(sys.argv[2]); S = int(sys.argv[3]); delta = Fr(sys.argv[4]); J = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    alphas = [Fr(a) for a in sys.argv[6].split(',')] if len(sys.argv) > 6 else []
    rng = np.random.default_rng(58)
    n = N // 2
    cand = []
    for _ in range(S):
        xs = cue_angles(N, rng) * N / (2 * PI)
        ys, g = hardcore(xs, float(N), float(delta) + 2e-5)
        q, gmin = rational_zeros(ys, N, delta)
        if gmin >= delta: cand.append(q)
    print(f'mode={mode} N={N} delta={delta} candidates={len(cand)}', flush=True)
    om = [2 * PI * j / N for j in range(n + 1)]
    svals = [Fr(1)] if mode in ('disc', 'cont1') else [Fr(k, 20) for k in range(10, 21)]
    cmode = mode in ('cont', 'cont1')
    cols = []; meta = []
    for ci, q in enumerate(cand):
        w, d, m1 = stats(np.array([float(x) for x in q]), float(N))
        for s in svals:
            sf = float(s)
            if not cmode:
                v = np.concatenate([w, [d]])
            else:
                v = np.array([np.sum(w * (sf * np.array(om)) ** (2 * j)) for j in range(J + 1)] + [d * sf ** 3] + Dfloat(q, N, sf, alphas))
            cols.append(v); meta.append((ci, s))
    if not cmode:
        tgt = np.array([1.0] + [2.0] * n) / (N + 1); m1t = float(np.sum(tgt * np.array(om) ** 2)); tgt = np.concatenate([tgt, [m1t]])
    else:
        tgt = np.array([PI ** (2 * j) / (2 * j + 1) for j in range(J + 1)] + [PI ** 2 / 3] + [1 - (math.sin(PI * float(a)) / (PI * float(a))) ** 2 for a in alphas])
    A = np.array(cols).T; R = 1.0 / tgt
    nr, nc = A.shape
    res = linprog(np.concatenate([np.zeros(nc), np.ones(2 * nr)]), A_eq=np.hstack([A * R[:, None], np.eye(nr), -np.eye(nr)]), b_eq=tgt * R,
                  bounds=[(0, None)] * (nc + 2 * nr), method='highs')
    x = res.x[:nc]
    print(f'LP infeasibility {res.fun:.3e}', flush=True)
    if res.fun > 1e-9: sys.exit(1)
    basis = [i for i in np.argsort(-x)[:nr]]  # top nr columns (basic solution has <= nr nonzeros)
    basis = [i for i in basis if x[i] > 0] + [i for i in basis if x[i] <= 0]
    # Arb square system
    M = arb_mat(nr, nr); b = arb_mat(nr, 1)
    pi = arb.pi()
    for c_, i in enumerate(basis):
        ci, s = meta[i]
        w, m0, d = arb_column(cand[ci], N)
        sA = arb(s.numerator) / s.denominator
        if not cmode:
            vec = [wj / m0 for wj in w] + [d / m0]
        else:
            vec = []
            for j in range(J + 1):
                acc = arb(0)
                for k in range(n + 1):
                    acc += w[k] * (sA * 2 * pi * k / N) ** (2 * j)
                vec.append(acc / m0)
            vec.append(d / m0 * sA ** 3)
            vec += [Darb(cand[ci], N, s, a) / m0 for a in alphas]
        for r_ in range(nr): M[r_, c_] = vec[r_]
    if not cmode:
        tv = [arb(1) / (N + 1)] + [arb(2) / (N + 1)] * n
        m1a = sum((tv[k] * (2 * pi * k / N) ** 2 for k in range(n + 1)), arb(0)); tv.append(m1a)
    else:
        tv = [pi ** (2 * j) / (2 * j + 1) for j in range(J + 1)] + [pi ** 2 / 3] + [1 - ((pi * arb(a.numerator) / a.denominator).sin() / (pi * arb(a.numerator) / a.denominator)) ** 2 for a in alphas]
    for r_ in range(nr): b[r_, 0] = tv[r_]
    try:
        sol = M.solve(b)
    except Exception as e:
        print('Arb solve failed:', e); sys.exit(2)
    ws = [sol[i, 0] for i in range(nr)]
    allpos = all(wv > 0 for wv in ws)
    lows = [float(wv.lower()) if hasattr(wv, 'lower') else float(wv.mid() - wv.rad()) for wv in ws]
    print(f'Arb: {nr}x{nr} system solved; all weights > 0 certified: {allpos}; min weight lower bound {min(lows):.3e}; max radius {max(float(wv.rad()) for wv in ws):.1e}', flush=True)
    out = dict(mode=mode, N=N, delta=str(delta), J=J, certified=bool(allpos), min_weight_lower=min(lows),
               columns=[dict(zeros=[str(z) for z in cand[meta[i][0]]], s=str(meta[i][1])) for i in basis],
               weights_mid=[float(wv.mid()) for wv in ws])
    json.dump(out, open(f'Riemann/Deep/ExploreLK/logs/cert_{mode}_N{N}_d{str(delta).replace("/", "_")}_J{J}_D{len(alphas)}.json', 'w'), indent=0)
