#!/usr/bin/env python3
"""c_reach.py -- explore15, adversarial test for transfer A/C (N6): the reach of Jensen criteria.

Toy 'Davenport-Heilbronn-type' function: Xi with ONE zero pair moved off the line,
    Theta~(w) = Theta(w) (1 - w/w~)(1 - w/conj w~) / (1 + w/a),  a = gamma_k^2,  w~ = (delta + i gamma_k)^2,
where 8 xi(1/2+z) = Theta(z^2) = sum gamma(n) z^{2n}/n! (Jensen9 normalisation).  Theta~ keeps the functional
equation and reality on the line, has no Euler product, and exactly one off-line pair (Re rho = 1/2 +- delta).
All its zeros lie in the sector |arg w - pi| <= th = 2 arctan(delta/gamma_k), so Chasse's sector lemma gives
J^{d,n} hyperbolic for all n and d <= D_sec = floor(1/sin^2 th) ~ gamma_k^2/(4 delta^2).
We find the true first failing degree d*(n) (J^{d,n} non-hyperbolic) by tracking the off-line pair of
J^{d,n}(X/d) from the Hurwitz limit down in d, and certify hyperbolicity just below d* by sign counts.

gamma(n) = 64 4^n n!/(2n)! b_n,  b_n = int_0^inf u^{2n} Phi(u) du  (positive integrand: relative precision
directly), by the trapezoidal rule (Phi even and entire; error e^{-2 pi a/h}); checked against the Arb
power series of xi for n <= 160 and by a second step size.
Run: uv run --no-project --with python-flint python3 -u c_reach.py
"""
import sys, time, math
from math import comb
from flint import arb, acb, arb_poly, acb_poly, ctx

def phi(u, M):
    pi = arb.pi(); e4 = (4 * u).exp(); e5 = (5 * u).exp(); e9 = (9 * u).exp()
    s = arb(0)
    for m in range(1, M + 1):
        m2 = m * m
        s += (2 * pi * pi * m2 * m2 * e9 - 3 * pi * m2 * e5) * (-pi * m2 * e4).exp()
    return s

def moments_trap(K, h, U, M):
    """b_n for n = 0..K by the trapezoidal rule on (0, U]"""
    acc = [arb(0)] * (K + 1)
    nmax = int(U / h)
    h = arb(h)                                  # the float h is exact in binary
    for i in range(1, nmax + 1):
        u = arb(i) * h
        f = phi(u, M)
        u2 = u * u
        p = f
        for n in range(K + 1):
            acc[n] += p
            p *= u2
    b = [a * h for a in acc]
    b[0] += h * phi(arb(0), M) / 2
    return b

def gammas_from_b(b):
    """int_0^inf Phi(u) cosh(zu) du = xi(1/2 + z/2)/8 (Polymath15 H_0), so gamma(n) = 64 4^n n!/(2n)! b_n"""
    out = []; r = arb(1)                      # r = 4^n n!/(2n)!
    for n in range(len(b)):
        if n: r = r * 4 * n / ((2 * n) * (2 * n - 1))
        out.append(64 * r * b[n])
    return out

def jensen(gs, d, n):
    return [comb(d, j) * gs[n + j] for j in range(d + 1)]

def horner(c, x):
    v = 0 * x; dv = 0 * x
    for a in reversed(c):
        dv = dv * x + v; v = v * x + a
    return v, dv

def newton(c, x, steps=300):
    """Newton on midpoints (no ball blow-up); c: list of acb midpoints"""
    x = acb(x.mid())
    tol = arb(2) ** (-(ctx.prec - 80))
    for _ in range(steps):
        v, dv = horner(c, x)
        dx = acb((v / dv).mid())
        x = acb((x - dx).mid())
        if abs(dx).mid() < abs(x).mid() * tol:
            break
    return x

def certified_real_count(c):
    """certified lower bound on the number of real zeros (sign changes at exact points between accurate
    approximations of the roots; ball arithmetic)"""
    P = arb_poly([arb(a) for a in c])
    tol = arb(2) ** (-(ctx.prec // 2))
    rr = acb_poly([acb(a.mid()) for a in c]).roots(tol=tol)
    xs = sorted([arb(r.real.mid()) for r in rr], key=lambda a: float(a))
    pts = [xs[0] * 2 - 1] + [arb(((xs[i] + xs[i + 1]) / 2).mid()) for i in range(len(xs) - 1)] + [arb(1)]
    sg = []
    for p in pts:
        v = P(p); sg.append(1 if v > 0 else (-1 if v < 0 else 0))
    return sum(1 for i in range(len(sg) - 1) if sg[i] * sg[i + 1] < 0)

def toy(gs, gk, delta, K):
    """gamma~(k), k <= K, of Theta~ for the zero w = -gk^2 of Theta replaced by w~, conj w~.
    G = Theta/(1 + w/a) by the stable backward recursion g_{k-1} = a (t_k - g_k) (t_k = gamma(k)/k!);
    then Theta~ = G (1 + q1 w + q2 w^2) with q1, q2 > 0 (no cancellation)."""
    Kx = len(gs) - 1
    gk = arb(gk)                                   # must be the zero to full working precision (else Theta - c)
    w = acb(arb(delta), gk) ** 2
    a = gk ** 2
    inv = 1 / w
    q1 = -2 * inv.real; q2 = abs(inv) ** 2              # Q(w) = (1 - w/w~)(1 - w/conj w~) = 1 + q1 w + q2 w^2
    facs = [arb(1)]
    for k in range(1, Kx + 1): facs.append(facs[-1] * k)
    t = [gs[k] / facs[k] for k in range(Kx + 1)]
    g = [arb(0)] * (Kx + 1)
    for k in range(Kx, 0, -1):
        g[k - 1] = a * (t[k] - g[k])
    out = []
    for k in range(K + 1):
        s = g[k] + (q1 * g[k - 1] if k >= 1 else 0) + (q2 * g[k - 2] if k >= 2 else 0)
        out.append(s * facs[k])
    # sanity: Theta~(0) = Theta(0); G(-a) finite
    return out, w

def deriv_zero(gt, w, n):
    """off-line zero of Theta~^{(n)}(w) = sum_k gt[n+k] w^k/k! near w~ (Newton on the truncated series)"""
    K = len(gt) - 1 - n
    c = []; f = arb(1)
    for k in range(K + 1):
        if k: f *= k
        c.append(acb((gt[n + k] / f).mid()))
    return newton(c, w)

def toy_multi(gs, pairs, K):
    """several zeros moved: pairs = [(gamma_k (arb), delta), ...]; returns (gamma~, list of w~)"""
    g = gs; ws = []
    for gk, delta in pairs:
        g, w = toy(g, gk, delta, len(g) - 1 - 2)
        ws.append(w)
    return g[:K + 1], ws

def first_failure(gt, w, n, dlo, dhi):
    """smallest d in (dlo, dhi] with J^{d,n} non-hyperbolic, by tracking the off-line pair from d = dhi down.
    Returns (d*, tracked zero at d*, tracked zero at d*-1 (real))."""
    # start: Hurwitz, J^{d,n}(X/d) -> Theta~^{(n)}(X); its off-line pair sits near w~ (n = 0) or w_n
    d = dhi
    wn = w if n == 0 else deriv_zero(gt, w, n)
    x = acb(wn.real / d, wn.imag / d)
    c = [acb(a.mid()) for a in jensen(gt, d, n)]
    x = newton(c, x)
    if abs(x.imag.mid()) < 1e-30 * abs(x.real.mid()):
        return None
    last_complex = d
    # step down in d, rescaling the zero (X ~ w/d)
    while d > dlo:
        step = max(1, (d - dlo) // 8) if abs(x.imag.mid()) > 1e-3 * abs(x.real.mid()) else 1
        d2 = d - step
        c = [acb(a.mid()) for a in jensen(gt, d2, n)]
        x2 = newton(c, acb(x.real * d / d2, x.imag * d / d2))
        if abs(x2.imag.mid()) < arb(10) ** (-40) * abs(x2.real.mid()):
            if step == 1:
                return d, x, x2
            # overshoot: bisect between d2 (real) and d (complex), restarting from the complex zero at d
            lo_d, xr = d2, x2
            while d - lo_d > 1:
                mid = (d + lo_d) // 2
                c = [acb(a.mid()) for a in jensen(gt, mid, n)]
                xm = newton(c, acb(x.real * d / mid, x.imag * d / mid))
                if abs(xm.imag.mid()) < arb(10) ** (-40) * abs(xm.real.mid()):
                    lo_d, xr = mid, xm
                else:
                    d, x = mid, xm
            return d, x, xr
        d, x = d2, x2
    return dlo, x, None

def main():
    t0 = time.time()
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 900
    digits = int(0.66 * K) + 80
    ctx.prec = int(digits * 3.33)
    M = int(math.sqrt(digits * 2.31 / math.pi)) + 3
    U = 2.0 if K <= 1200 else 2.15
    a = math.pi / 16
    h = 2 * math.pi * a / ((digits + 200) * math.log(10))
    print(f'K = {K}, {digits} digits ({ctx.prec} bits), M = {M}, U = {U}, h = {h:.3e} ({int(U / h)} nodes)')
    KX = K + 60                                   # extra terms for the backward recursion in toy()
    b1 = moments_trap(KX, h, U, M)
    g1 = gammas_from_b(b1)
    print(f'  moments, step h: {time.time() - t0:.0f} s')
    if 'nocheck' in sys.argv:
        b2 = b1
    else:
        b2 = moments_trap(KX, h * 0.85, U, M)
    g2 = gammas_from_b(b2)
    err = max(float((abs((g1[n] - g2[n]) / g1[n]).mid().log() / arb(10).log())) if g1[n] != g2[n] else -9999 for n in range(KX + 1))
    print(f'  second step 0.85h: max log10 rel diff {err:.1f}  ({time.time() - t0:.0f} s)')
    # check against the Arb power series for n <= 160 (Jensen9 method)
    sys.path.insert(0, '.')
    from b_fragility import gammas as gser
    prec_keep = ctx.prec
    gs = gser(160, 4400)
    ctx.prec = prec_keep
    err2 = max(float((abs((g1[n] - gs[n]) / gs[n]).mid() + arb(10) ** -1300).log() / arb(10).log()) for n in range(161))
    print(f'  vs Arb power series, n <= 160: max log10 rel diff {err2:.1f}; gamma(0) = {g1[0].str(20)}')
    zz = acb.zeta_zeros(1, 3)                     # 1/2 + i gamma_k at full working precision (Arb)
    zeros = {k + 1: zz[k].imag for k in range(3)}
    print('  gamma_1..3 =', ', '.join(z.str(25) for z in zeros.values()))
    print('toy: zero #k moved to 1/2 + delta + i gamma_k.  D_sec = floor(1/sin^2(2 arctan(delta/gamma_k))) (Chasse)')
    print('  k  delta  |  D_sec  |  d*(n=0)  d*/D_sec | certified: J^{d*-1,0} real zeros / degree | Im X at d*')
    res = []
    table = [(1, 0.45), (1, 0.35), (1, 0.25), (2, 0.45), (2, 0.35), (3, 0.45)] if 'shift' not in sys.argv else []
    extra = [a[4:] for a in sys.argv if a.startswith('cfg=')]
    if extra:
        table = []
        zz2 = acb.zeta_zeros(1, 25)
        for cfg in extra:
            pairs = [(int(p.split(':')[0]), float(p.split(':')[1])) for p in cfg.split('+')]
            gt, ws = toy_multi(g1, [(zz2[k - 1].imag, dl) for k, dl in pairs], K)
            Ds = [int(1 / math.sin(2 * math.atan(dl / float(zz2[k - 1].imag.mid()))) ** 2) for k, dl in pairs]
            j = Ds.index(min(Ds)); Dsec = Ds[j]; w = ws[j]
            dhi = min(K, int(1.25 * Dsec) + 2)
            r = first_failure(gt, w, 0, Dsec // 2, dhi)
            if r is None:
                print(f'  cfg {cfg}: D_sec = {Ds}; no complex pair at d = {dhi}   ({time.time() - t0:.0f} s)'); continue
            print(f'  cfg {cfg}: D_sec = {Ds}; first failing d = {r[0]} (tracked pair of the steepest zero), '
                  f'excess {r[0] - Dsec} = {100 * (r[0] - Dsec) / Dsec:.2f}%; Im X = {float(r[1].imag.mid()):.3e}   ({time.time() - t0:.0f} s)')
    for k, delta in table:
        gk = zeros[k]
        th = 2 * math.atan(delta / float(gk.mid()))
        Dsec = int(1 / math.sin(th) ** 2)
        if int(1.25 * Dsec) + 2 > K:
            print(f'  {k}  {delta:.2f}  | {Dsec:6d}  |  skipped (needs K >= {int(1.25 * Dsec) + 2})'); continue
        gt, w = toy(g1, gk, delta, K)
        dhi = min(K, int(1.25 * Dsec) + 2)
        r = first_failure(gt, w, 0, Dsec, dhi)
        if r is None:
            print(f'  {k}  {delta:.2f}  | {Dsec:6d}  |  no complex pair at d = {dhi}'); continue
        dstar, xs, xr = r
        cert = ''
        if dstar - 1 <= 420:
            cnt = certified_real_count(jensen(gt, dstar - 1, 0))
            cert = f'{cnt}/{dstar - 1}'
        print(f'  {k}  {delta:.2f}  | {Dsec:6d}  |  {dstar:6d}   {dstar / Dsec:6.3f}  | {cert:>12s} | {float(xs.imag.mid()):.3e}'
              f'   ({time.time() - t0:.0f} s)')
        res.append((k, delta, Dsec, dstar))
    # shift dependence for the first toy
    if extra:
        print(f'done ({time.time() - t0:.0f} s)'); return
    print('shift dependence (k=1, delta=0.45): first failing d for n = 0, 1, 2, 3, 5, 10; off-line zero of Theta~^(n)')
    gt, w = toy(g1, zeros[1], 0.45, K)
    Dsec = int(1 / math.sin(2 * math.atan(0.45 / float(zeros[1].mid()))) ** 2)
    for n in [0, 1, 2, 3, 5, 10]:
        if n + 2 * Dsec > K: break
        wn = w if n == 0 else deriv_zero(gt, w, n)
        thn = math.atan(float(abs(wn.imag.mid() / wn.real.mid())))
        if thn < 1e-30:
            cnt = certified_real_count(jensen(gt, 300, n))
            print(f'  n = {n:3d}: Newton from w~ lands on the REAL zero {float(wn.real.mid()):.4f} of Theta~^(n)'
                  f' (Polya-Wiman: derivatives lose the off-line pair); certified J^(300,{n}): {cnt}/300 real'
                  f'   ({time.time() - t0:.0f} s)')
            continue
        Dn = int(1 / math.sin(thn) ** 2)
        r = first_failure(gt, w, n, Dsec // 2, 2 * Dsec)
        print(f'  n = {n:3d}: zero of Theta~^(n) at {float(wn.real.mid()):.4f} + {float(wn.imag.mid()):.4f}i,'
              f' 1/sin^2(angle) = {Dn}; d*(n) = {r[0] if r else None}   ({time.time() - t0:.0f} s)')
    print(f'done ({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
