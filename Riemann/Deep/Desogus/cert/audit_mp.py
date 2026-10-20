"""Independent audit (mpmath only) of Desogus Thm 8.6 / (201).
Weil form on odd real v supported in (-a,a):
  W(v) = sum over ALL zeros rho=1/2+i*gam of |vhat(gam)|^2,  vhat(r)=int v(y) e^{iry} dy   (paper eq. (3))
explicit side (Guinand-Weil, psi(s) = -euler + int_0^inf (e^-t - e^-st)/(1-e^-t) dt, t->2u):
  W = -2 U^2 - 2 sum_{n<e^{2a}} Lam(n) n^{-1/2} c(log n)
      + c(0)(-euler - log pi - log(1-e^{-4a})) + int_0^{2a} 2(c(0)e^{-2t} - c(t)e^{-t/2})/(1-e^{-2t}) dt,
  c(t)=int v(y)v(y-t)dy, U=int v e^{y/2} dy.
Step 0 checks this against zeta zeros on smooth odd tests. Basis: odd monomials (y/a)^j."""
import sys, time
from fractions import Fraction as Fr
from math import comb
import mpmath as mp

mp.mp.dps = 50

def vm(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return mp.log(p) if m == 1 else None

def corr_poly(i, j):
    """int_{s-1}^1 x^i (x-s)^j dx as exact poly in s (list of Fractions), x = y/a, s = t/a, 0<=s<=2."""
    out = {}
    for l in range(j + 1):
        cl = comb(j, l) * (-1) ** (j - l)       # coefficient of x^l s^(j-l)
        p = i + l + 1
        # (1 - (s-1)^p)/p
        out[j - l] = out.get(j - l, Fr(0)) + Fr(cl, p)
        for r in range(p + 1):
            c = comb(p, r) * (-1) ** (p - r)    # (s-1)^p = sum C(p,r) s^r (-1)^(p-r)
            out[j - l + r] = out.get(j - l + r, Fr(0)) - Fr(cl * c, p)
    deg = max(out)
    return [out.get(m, Fr(0)) for m in range(deg + 1)]

class Window:
    def __init__(self, k1):          # window log(k1): a = log(k1)/2
        self.k1 = k1; a = self.a = mp.log(k1) / 2
        self.const = -mp.euler - mp.log(mp.pi) - mp.log(1 - mp.e ** (-4 * a))
        self.J = {}
        f0 = lambda t: 2 * (mp.expm1(-2 * t) - mp.expm1(-t / 2)) / (-mp.expm1(-2 * t))
        self.J[0] = mp.quad(f0, [0, a, 2 * a])
        self.primes = [(n, vm(n)) for n in range(2, k1) if vm(n) is not None]   # n < e^{2a} = k1

    def Jm(self, m):
        if m not in self.J:
            a = self.a
            self.J[m] = -mp.quad(lambda t: 2 * (t / a) ** m * mp.e ** (-t / 2) / (-mp.expm1(-2 * t)), [0, a, 2 * a])
        return self.J[m]

    def form_from_corr(self, C):
        """C = poly coeffs (in s=t/a) of symmetrised correlation (times a already). Returns arch - primes."""
        a = self.a
        arch = C[0] * (self.const + self.J[0]) + sum(C[m] * self.Jm(m) for m in range(1, len(C)))
        cval = lambda t: sum(C[m] * (t / a) ** m for m in range(len(C)))
        pr = sum(L / mp.sqrt(n) * 2 * cval(mp.log(n)) for n, L in self.primes)
        return arch - pr

    def matrix(self, degs):
        a = self.a; N = len(degs)
        U = [a * mp.quad(lambda x: x ** i * mp.e ** (a * x / 2), [-1, 0, 1]) for i in degs]
        M = mp.matrix(N, N)
        for p, i in enumerate(degs):
            for q, j in enumerate(degs):
                if q < p: M[p, q] = M[q, p]; continue
                gij = corr_poly(i, j); gji = corr_poly(j, i)
                L = max(len(gij), len(gji))
                fr = [(gij[m] if m < len(gij) else 0) + (gji[m] if m < len(gji) else 0) for m in range(L)]
                C = [a * mp.mpf(Fr(f).numerator) / Fr(f).denominator / 2 for f in fr]
                M[p, q] = self.form_from_corr(C) - 2 * U[p] * U[q]
        return M

def W_direct(win, coeffs_by_deg):
    """W(v) for v(y)=sum c_j (y/a)^j via the explicit side with quadrature of c(t) (no matrix)."""
    a = win.a; mp.mp.dps = 25
    v = lambda y: sum(c * (y / a) ** j for j, c in coeffs_by_deg.items())
    ct = lambda t: mp.quad(lambda y: v(y) * v(y - t), [-a + t, a])
    c0 = ct(0)
    arch = c0 * win.const + mp.quad(lambda t: 2 * (c0 * mp.e ** (-2 * t) - ct(t) * mp.e ** (-t / 2)) / (-mp.expm1(-2 * t)), [0, a, 2 * a])
    U = mp.quad(lambda y: v(y) * mp.e ** (y / 2), [-a, 0, a])
    pr = sum(L / mp.sqrt(n) * 2 * ct(mp.log(n)) for n, L in win.primes)
    mp.mp.dps = 50
    return arch - 2 * U ** 2 - pr

def W_zeros(win, coeffs_by_deg, nz):
    a = win.a
    v = lambda y: sum(c * (y / a) ** j for j, c in coeffs_by_deg.items())
    tot = 0; last = 0
    for n in range(1, nz + 1):
        g = mp.im(mp.zetazero(n))
        vh = 2 * mp.quad(lambda y: v(y) * mp.sin(g * y), [0, a / 2, a])    # vhat = i*this
        last = 2 * vh ** 2; tot += last
    return tot, last

def step0():
    mp.mp.dps = 30
    win = Window(8)
    for pw in (6, 10):
        # v = x (1-x^2)^pw, x=y/a: expand
        cf = {1 + 2 * r: mp.mpf(comb(pw, r) * (-1) ** r) for r in range(pw + 1)}
        degs = sorted(cf); M = win.matrix(degs)
        c = mp.matrix([cf[d] for d in degs])
        Wm = (c.T * M * c)[0]
        Wz, last = W_zeros(win, cf, 40)
        print(f"[conv] window log8, v=x(1-x^2)^{pw}: W(matrix)={mp.nstr(Wm, 12)}  2*sum_{{40 zeros}}|vhat|^2={mp.nstr(Wz, 12)}  (last term {mp.nstr(last, 3)})", flush=True)
    mp.mp.dps = 50

def ground_ratio(win, k, degs, M, which):
    a = win.a; ak = mp.log(k) / 2
    if which == 'exp':      # sgn(y) e^{-|y|/2} on collar  (odd extension of t^{-1/2})
        prof = lambda y: mp.e ** (-y / 2); n2 = 2 * (mp.e ** (-ak) - mp.e ** (-a))
    else:                   # sinh(y/2) on collar  (odd part of t^{-1/2}|collar ~ polar t_k)
        prof = lambda y: mp.sinh(y / 2); n2 = mp.sinh(a) - mp.sinh(ak) - (a - ak)
    e = mp.matrix([2 * mp.quad(lambda y: (y / a) ** i * prof(y), [ak, a]) / mp.sqrt(n2) for i in degs])
    x = mp.lu_solve(M, e)
    return 1 / (e.T * x)[0], x, e

def Q_true_e(k):
    """Q(e_k), e_k = normalised sgn(y)e^{-|y|/2} on a_k<|y|<a_{k+1}; c(t) in closed-form quadrature."""
    mp.mp.dps = 30
    win = Window(k + 1); a = win.a; ak = mp.log(k) / 2
    nrm2 = 2 * (mp.e ** (-ak) - mp.e ** (-a))
    pieces = [((ak, a), 1), ((-a, -ak), -1)]
    def ct(t):
        s = 0
        for (l1, h1), s1 in pieces:
            for (l2, h2), s2 in pieces:
                lo, hi = max(l1, l2 + t), min(h1, h2 + t)
                if hi > lo:
                    s += s1 * s2 * mp.quad(lambda y: mp.e ** (-(abs(y) + abs(y - t)) / 2), [lo, hi] if not (lo < t < hi) else [lo, t, hi])
        return s / nrm2
    c0 = ct(0)
    br = sorted(set([mp.mpf(0), a - ak, 2 * ak, a + ak, 2 * a]))
    arch = c0 * win.const + mp.quad(lambda t: 2 * (c0 * mp.e ** (-2 * t) - ct(t) * mp.e ** (-t / 2)) / (-mp.expm1(-2 * t)), br)
    U = ((a - ak) - (mp.e ** (-ak) - mp.e ** (-a))) / mp.sqrt(nrm2)
    pr = sum(L / mp.sqrt(n) * 2 * ct(mp.log(n)) for n, L in win.primes)
    mp.mp.dps = 50
    return arch - 2 * U ** 2 - pr, c0

if __name__ == '__main__':
    t0 = time.time()
    step0(); print(f"  t={time.time()-t0:.0f}s", flush=True)
    for k in (7, 12):
        Q, c0 = Q_true_e(k)
        print(f"[diag] k={k}: ||e_k||^2={mp.nstr(c0, 10)}  Q(e_k)={mp.nstr(Q, 10)}  (paper: B_k={mp.nstr(0.5 + mp.log(1 / mp.log(1 + mp.mpf(1) / k)), 6)}, g_k>=1.878)  t={time.time()-t0:.0f}s", flush=True)
    k = 7; win = Window(k + 1)
    for D in (9, 13, 19):
        degs = list(range(1, D + 1, 2)); M = win.matrix(degs)
        ev = mp.eigsy(M)[0]
        for which in ('exp', 'sinh'):
            r, x, e = ground_ratio(win, k, degs, M, which)
            msg = f"[gal] k=7 window log8 deg{D} e={which}: Q(v)/<e,v>^2 = {mp.nstr(r, 10)}  lam_min(M)={mp.nstr(min(ev), 4)}"
            if D == 19 and which == 'exp':
                Wd = W_direct(win, {d: x[i] for i, d in enumerate(degs)})
                ip = (e.T * x)[0]
                msg += f"  | direct-quad recheck: {mp.nstr(Wd / ip ** 2, 10)}"
            print(msg + f"  t={time.time()-t0:.0f}s", flush=True)
