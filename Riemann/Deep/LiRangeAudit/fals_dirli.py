#!/usr/bin/env python3
"""fals_dirli.py -- LiRangeAudit falsification for DirLi. Imports Riemann/Deep/DirLi/certify.py read-only.
D1 corner pairs: Re[(1-z_rho^n)+(1-z_rho*^n)] >= G_n(|y|) - q_n(|y|), rho = 1/2+d+iy, rho* = 1-conj(rho), d = a edge, y = +-(t0+1e-6), +-1.001 t0.
D2 Re lambda_n(chi), n <= 200, by the Cauchy integral of d/dz log xi(1/(1-z), chi) (mpmath, |z| = 0.8, 256 points, 40 digits),
   chi mod 3 (real odd), chi mod 5 with chi(2) = i (complex odd), chi mod 5 Legendre (real even); compared with DirLi's certified lower bound
   max(gainC(n), gainA(n,n)) - loss(n) at t0 replaced by 50 and 100, a = 1/2 and 3/8; chain N(t0) for these small t0.
"""
import sys, time, math
sys.path.insert(0, 'Riemann/Deep/DirLi')
import certify as C
from fractions import Fraction
from flint import arb, acb
import mpmath as mp
t_start = time.time()
def D1():
    cnt = 0; worst = None
    for a in (C.R(1, 2), C.R(3, 8), C.R(1, 8)):
        for t0 in (C.R(50), C.R(200), C.R(100000000) / 3):
            for y in (t0 + C.R(1, 10 ** 6), t0 * C.R(1001, 1000)):
                for sg in (1, -1):
                    for d in (a, a / 2):
                        if d >= C.R(1, 2): d = C.R(1, 2) - C.R(1, 10 ** 30)
                        rho = acb(C.R(1, 2) + d, sg * y); z = 1 - 1 / rho; zs = 1 - 1 / (1 - rho.conjugate())
                        lz, lzs = z.log(), zs.log(); nmax = 80 * float(y) ** 2 / float(a)
                        for k in range(170):
                            n = int(nmax ** (k / 169.0)) + (k % 3)
                            lhs = (2 - (n * lz).exp() - (n * lzs).exp()).real
                            qq = C.qn(n, a, y); dl = lhs - (C.G(n, y) - qq); cnt += 1
                            if not (dl > 0):
                                print("D1 VIOLATION/UNDECIDED", a, t0, y, sg, d, n, dl); return
                            rel = float(dl / qq)
                            if worst is None or rel < worst[0]: worst = (rel, float(a), float(y) * sg, float(d), n)
    print("D1 corner pairs (both signs of gamma): %d certified cases; min slack/q_n = %.3g at a=%.3g gamma=%.8g delta=%.3g n=%d" % ((cnt,) + worst), flush=True)
def re_lambda(q, chi, par, NMAX=200, M=256, rad=0.8):
    mp.mp.dps = 40
    def f(z):
        s = 1 / (1 - z)
        L = sum(chi[b] * mp.zeta(s, mp.mpf(b) / q) for b in chi) * mp.power(q, -s)
        dL = -mp.log(q) * L + sum(chi[b] * mp.zeta(s, mp.mpf(b) / q, 1) for b in chi) * mp.power(q, -s)
        return s * s * (mp.log(mp.mpf(q) / mp.pi) / 2 + mp.digamma((s + par) / 2) / 2 + dL / L)
    zk = [rad * mp.expjpi(2 * mp.mpf(k) / M) for k in range(M)]; fk = [f(z) for z in zk]
    zinv = [1 / z for z in zk]; pw = [mp.mpf(1)] * M; out = [None]
    for n in range(1, NMAX + 1):
        out.append(float(mp.re(mp.fsum(fk[k] * pw[k] for k in range(M))) / M)); pw = [pw[k] * zinv[k] for k in range(M)]
    return out
def lower(n, p):
    best = None
    for gg in ((C.gainC(n, p) if n <= 400 else None), (C.gainA(n, n, p) if n >= 3 else None)):
        if gg is not None:
            v = gg - C.loss(n, p)
            if best is None or v.upper() > best.upper(): best = v
    return best
def D2():
    I = mp.mpc(0, 1)
    chars = [(3, {1: 1, 2: -1}, 1, 'chi_-3 (real odd)'), (5, {1: 1, 2: I, 4: -1, 3: -I}, 1, 'chi mod 5, chi(2)=i (complex odd)'),
             (5, {1: 1, 2: -1, 3: -1, 4: 1}, 0, 'chi mod 5 Legendre (real even)')]
    for q, chi, par, name in chars:
        lam = re_lambda(q, chi, par)
        mn = min(lam[1:]); print("D2 %s: Re lambda_1 = %.12f, Re lambda_100 = %.6f, Re lambda_200 = %.6f, min = %.6g; (n/2)(log n + log(q/2pi) - 1 + gamma) at 200: %.4f"
                                 % (name, lam[1], lam[100], lam[200], mn, 100 * (math.log(200) + math.log(q / (2 * math.pi)) - 1 + 0.5772156649)), flush=True)
        for t0 in (50, 100):
            for a in (Fraction(1, 2), Fraction(3, 8)):
                p = C.P(q, q, Fraction(t0), a, 'test'); bad = []; mr = 0.0
                for n in range(1, 201):
                    lb = lower(n, p)
                    if lb is None: continue
                    if not (float(lb.upper()) <= lam[n]): bad.append(n)
                    mr = max(mr, float(lb.upper()) / lam[n])
                try:
                    N, n0, steps, w = C.certify(p)
                except SystemExit as e:
                    print('   chain stopped:', e); N, n0 = 0, 0
                print("   t0=%d a=%s: lower bound <= Re lambda_n for all n<=200: %s (max ratio %.3f); chain N(t0) = %d = %.2f t0^2 (regime C n<=%d)"
                      % (t0, a, 'yes' if not bad else 'NO ' + str(bad[:5]), mr, N, N / t0 ** 2, n0), flush=True)
D1(); D2(); print("total %.0f s" % (time.time() - t_start))
