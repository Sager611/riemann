#!/usr/bin/env python3
"""
sanity.py -- numerical sanity checks for certify.py (not part of the proof).
1. Normalisation: lambda_n of zeta from the Taylor series of log xi at s = 1 (lambda_1 = 0.0230957..., lambda_2 = 0.0923457...).
2. For every primitive chi mod q in {3,4,5,7,8,11,12}: Re lambda_n(chi), 1 <= n <= NMAX, from the Taylor series of
   log xi(s, chi) at s = 1, xi(s,chi) = (q/pi)^((s+a)/2) Gamma((s+a)/2) L(s,chi):
   log xi(1 + x) = sum a_k x^k  =>  coefficient of z^n in log xi(1/(1-z)) is c_n/n with c_n = n sum_k C(n-1,k-1) a_k,
   and Re c_n = Re lambda_n(chi) (an imaginary linear term i beta s of the Hadamard product only changes Im c_n).
   Check: the certified lower bound of certify.py (regimes C/A minus loss) is below Re lambda_n(chi) for every n.
3. For chi mod 3: zeros up to T = 300 (sign changes of Z(t, chi)); sum_{|gamma|<=T} (1 - cos(n theta)) plus the
   expected tail n^2 (log(qT/2pi)+1)/(2 pi T) reproduces lambda_n for n <= 5.
Run: uv run --no-project --with python-flint --with numpy python3 sanity.py
"""
import sys, time, math
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import certify as C
from flint import arb, acb, ctx, acb_series, arb_series, dirichlet_char

def binom_transform(co, NMAX):
    lam = [None]
    for n in range(1, NMAX + 1):
        s = acb(0); c = arb(1)                      # c = C(n-1, k-1)
        for k in range(1, n + 1):
            if k > 1:
                c = c * (n - k + 1) / (k - 1)
            s += c * co[k]
        lam.append(n * s)
    return lam

def zeta_check(NMAX=20):
    old = ctx.prec; ctx.prec = 200; ctx.cap = NMAX + 1
    x = acb_series([0, 1], prec=NMAX + 1)
    dz = acb_series.zeta(acb_series([1, 1], prec=NMAX + 1), deflate=True)   # zeta(1+x) - 1/x
    f = (acb_series([1]) + x * dz).log() + (acb_series([1]) + x).log() - acb_series([acb(2).log()]) \
        - (acb_series([1]) + x) * acb(arb.pi().log() / 2) + ((acb_series([1]) + x) / 2).lgamma()
    co = f.coeffs() + [acb(0)] * (NMAX + 2)
    lam = binom_transform(co, NMAX)
    ctx.prec = old
    return lam

def lambda_chi(chi, q, NMAX):
    ctx.cap = NMAX + 1                                # series truncation length (default 10 would truncate!)
    a = 0 if chi.parity() == 0 else 1
    s = acb_series([1, 1], prec=NMAX + 1)
    L = acb_series.dirichlet_l(s, chi)
    one = acb_series([1], prec=NMAX + 1)
    f = L.log() + ((s + a) / 2).lgamma() + (s + a) / 2 * acb((arb(q) / arb.pi()).log())
    co = f.coeffs() + [acb(0)] * (NMAX + 2)
    return binom_transform(co, NMAX)

def lower_bound(n, p):
    gs = []
    gA = C.gainA(n, n, p)
    if gA is not None:
        gs.append(gA)
    gs.append(C.gainC(n, p))
    g = max(gs, key=lambda v: float(v))
    return g - C.loss(n, p)

if __name__ == "__main__":
    t0 = time.time()
    lz = zeta_check()
    print("zeta: lambda_1 = %s, lambda_2 = %s, lambda_3 = %s" % (lz[1].real.str(12), lz[2].real.str(12), lz[3].real.str(12)))
    assert abs(float(lz[1].real) - 0.0230957089661210338) < 1e-12 and abs(float(lz[2].real) - 0.0923457352280611) < 1e-12
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    ctx.prec = int(1.2 * NMAX) + 120
    worst_ratio = 0.0
    for q in [3, 4, 5, 7, 8, 11, 12]:
        p = C.P(q, q, C.t0_of(q), C.A_HALF, 'q=%d' % q)
        old = ctx.prec; ctx.prec = 128
        LB = [None] + [lower_bound(n, p) for n in range(1, NMAX + 1)]
        ctx.prec = old
        for k in range(1, q + 1):
            if math.gcd(k, q) != 1:
                continue
            chi = dirichlet_char(q, k)
            if not chi.is_primitive():
                continue
            lam = lambda_chi(chi, q, NMAX)
            re = [None] + [lam[n].real for n in range(1, NMAX + 1)]
            mn = min(range(1, NMAX + 1), key=lambda n: float(re[n]))
            for n in range(1, NMAX + 1):
                if not (re[n] - LB[n] > 0):
                    raise SystemExit("FAIL: lower bound above Re lambda_n: q=%d chi=%d n=%d" % (q, k, n))
                worst_ratio = max(worst_ratio, float(LB[n] / re[n]))
            n = NMAX
            asym = n / 2 * (math.log(n) + 0.5772156649015329 - 1 + math.log(q / (2 * math.pi)))
            print("q=%-3d chi=%-3d parity=%d real=%d: Re lambda_1 = %s, min_{n<=%d} Re lambda_n at n=%d, "
                  "Re lambda_%d = %.4f (GRH main term (n/2)(log n + gamma - 1 + log(q/2pi)) = %.4f), Im lambda_1 = %s; "
                  "LB_1 = %s"
                  % (q, k, chi.parity(), chi.is_real(), re[1].str(10), NMAX, mn, NMAX, float(re[NMAX]), asym,
                     lam[1].imag.str(5), LB[1].str(5)), flush=True)
    print("certified lower bound <= Re lambda_n(chi) for all n <= %d, all primitive chi mod q in {3,4,5,7,8,11,12}; "
          "max ratio LB/lambda = %.3f" % (NMAX, worst_ratio))
    # zero-sum check for the real character mod 3
    ctx.prec = 64
    chi = dirichlet_char(3, 2)
    T = 300.0; h = 0.02; zeros = []
    t = 0.01; zp = chi.hardy_z(arb(t)).real if hasattr(chi.hardy_z(arb(t)), 'real') else chi.hardy_z(arb(t))
    def Z(t):
        v = chi.hardy_z(arb(t))
        return float(v.real if hasattr(v, 'real') else v)
    zp = Z(t)
    while t < T:
        t2 = t + h; z2 = Z(t2)
        if zp * z2 < 0:
            lo, hi = t, t2
            for _ in range(40):
                mid = (lo + hi) / 2
                if Z(mid) * Z(lo) <= 0: hi = mid
                else: lo = mid
            zeros.append((lo + hi) / 2)
        t, zp = t2, z2
    Nt = (T / math.pi) * math.log(3 * T / (2 * math.pi * math.e)) + 0.25
    print("chi mod 3: %d positive zeros found up to T = %g (both signs: %d; BMOR main term %.1f); first zeros %s"
          % (len(zeros), T, 2 * len(zeros), Nt, ["%.6f" % z for z in zeros[:3]]))
    ctx.prec = int(1.2 * 10) + 120
    lam = lambda_chi(chi, 3, 10)
    for n in range(1, 6):
        S = sum(2 * (1 - math.cos(n * 2 * math.atan(1 / (2 * g)))) for g in zeros)
        tail = n * n * (math.log(3 * T / (2 * math.pi)) + 1) / (2 * math.pi * T)
        print("  n=%d: Taylor lambda_n = %.8f, zero sum (|gamma|<=T) = %.8f, + tail estimate %.8f = %.8f"
              % (n, float(lam[n].real), S, tail, S + tail))
    print("runtime %.1f s" % (time.time() - t0))
