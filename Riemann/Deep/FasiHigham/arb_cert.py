#!/usr/bin/env python3
"""Arb (python-flint) certificates behind NOTE.md; optional, verify.py needs only mpmath.

  python3 arb_cert.py diag M          c_n(m,m) > 0 for all n, odd 3 <= m <= M
  python3 arb_cert.py near M1 M2      k = m-4,-2,-1,+1,+2,+4, odd M1 <= m <= M2
  python3 arb_cert.py col  M1 M2      every k < kmin(m) (verify.py table), odd M1 <= m <= M2
  python3 arb_cert.py kmin M1 M2 [K]  least k <= K with a negative c_n, n <= 4(k+m)+200

Method: 1/q(sz) as an Arb power series (rigorous balls) up to n1; certified root
balls of Q = b_0 q; D1 check; c_n >= (R_1 - sum_j |R_j| (R/|z_j|)^n) R^-n > 0 for
n >= n1, with R_j = -b_0/(z_j Q'(z_j)).  Needs: pip install python-flint.
"""
import sys, time
import flint
from math import comb
flint.ctx.cap = 10 ** 7

def bcoef(k, m):
    rising = [1]
    for i in range(1, m + 1): rising.append(rising[-1] * (k + i))
    return [(-1) ** j * comb(m, j) * rising[m - j] for j in range(m + 1)]

def series(k, m, L, s, prec):
    """Arb balls for c_n(k,m) s^n, n < L"""
    flint.ctx.prec = prec
    b = bcoef(k, m)
    S = flint.arb_series([flint.arb(flint.fmpq(b[j] * s ** j, b[0])) for j in range(m + 1)], prec=L)
    return (1 / S).coeffs()

def roots(k, m, prec):
    b = bcoef(k, m)
    flint.ctx.prec = prec
    P = flint.acb_poly([flint.acb(x) for x in b])
    rts = P.roots(tol=flint.arb(2) ** (-prec // 2), maxprec=400000)
    return sorted(rts, key=lambda z: abs(z).mid()), b, P

def tail(k, m):
    """n1 with c_n > 0 for all n >= n1, or None"""
    prec = max(256, 16 * m)
    while prec <= 64 * m + 4096:
        rts, b, P = roots(k, m, prec)
        dP = P.derivative(); r = rts[0]; ar = abs(r)
        if all(abs(z) > ar for z in rts[1:]):
            R1 = -b[0] / (r * dP(r))
            if R1.real > 0:
                terms = [(abs(-b[0] / (z * dP(z))) / abs(R1), ar / abs(z)) for z in rts[1:]]
                T = lambda n: sum((a * g ** n for a, g in terms), flint.arb(0))
                if T(0) < 1: return 0
                lo, hi = 0, 1
                while not (T(hi) < 1) and hi < 10 ** 8: hi *= 2
                if T(hi) < 1:
                    while hi - lo > 1:
                        mid = (lo + hi) // 2
                        if T(mid) < 1: hi = mid
                        else: lo = mid
                    return hi
        prec *= 2
    return None

def certify(k, m):
    """('POS', n1) if c_n(k,m) > 0 for every n; ('NEG', n) at the first negative n < n1; else ('UNKNOWN',)"""
    n1 = tail(k, m)
    if n1 is None: return ('UNKNOWN',)
    rts, _, _ = roots(k, m, 128)
    s = max(1, int(float(rts[0].real.mid()) + 0.5))
    L, prec = n1 + 1, 64 + 4 * m + (n1 + 1) // 8
    while prec < 10 ** 6:
        cs = series(k, m, L, s, prec)
        if all(c > 0 or c < 0 for c in cs):
            neg = [n for n, c in enumerate(cs) if c < 0]
            return ('NEG', neg[0]) if neg else ('POS', n1)
        prec *= 2
    return ('UNKNOWN',)

def first_neg(k, m, L):
    s = max(1, int(float(roots(k, m, 128)[0][0].real.mid()) + 0.5))
    prec = 64 + 2 * m + L // 4
    while True:
        cs = series(k, m, L, s, prec)
        if all(c > 0 or c < 0 for c in cs): break
        prec *= 2
    return next((n for n, c in enumerate(cs) if c < 0), None)

if __name__ == '__main__':
    t0 = time.time(); mode = sys.argv[1]; a = list(map(int, sys.argv[2:]))
    if mode == 'diag':
        for m in range(3, a[0] + 1, 2): print(m, certify(m, m), '%.0fs' % (time.time() - t0), flush=True)
    elif mode == 'near':
        for m in range(a[0], a[1] + 1, 2):
            print(m, [(d, certify(m + d, m)[0]) for d in (-4, -2, -1, 1, 2, 4)], '%.0fs' % (time.time() - t0), flush=True)
    elif mode == 'col':
        import verify
        for m in range(a[0], a[1] + 1, 2):
            km = verify.KMIN[m][0]
            bad = [k for k in range(km) if certify(k, m)[0] != 'POS']
            print(m, 'kmin', km, 'all smaller k positive' if not bad else bad, '%.0fs' % (time.time() - t0), flush=True)
    elif mode == 'kmin':
        K = a[2] if len(a) > 2 else 1000
        for m in range(a[0], a[1] + 1, 2):
            hit = next(((k, n) for k in range(K + 1) for n in [first_neg(k, m, 4 * (k + m) + 200)] if n is not None), None)
            print(m, hit, '%.0fs' % (time.time() - t0), flush=True)
