#!/usr/bin/env python3
"""odlyzko31 independent checks (float / mpmath; not part of the certificate).
(1) Weil-Poitou identity for K = Q (n = r1 = 1, d = 1) with F = S(x/b) sech(w x): 0 = Gamma terms - E + zero sum + prime sum,
    zero sum over the first NZ Riemann zeros (mpmath.zetazero), tail estimated.
(2) Frequency-domain check of the Gamma terms: (1/pi) int Fhat(t) (log 2pi - Re psi(1/2+it)) dt = 2 (gamma + log 8pi - Ic) and
    (1/(2pi)) int Fhat(t) (log pi - Re psi(1/4+it/2)) dt = gamma + pi/2 + log 8pi - Ic - Ir.
(3) Reproduction of Odlyzko 1976 Tables 3/4 (A, B, E per b) and Tables 1/2 (per n) from the same formula."""
import math, sys, time
import numpy as np
import mpmath as mp
from kern import consts, SHAPES, sech, EG, L8P

def fhat(shape, b, w, ts, m=6000):
    # Fhat(t) = int_R F(x) e^{itx} dx = 2 int_0^{2b} F(x) cos(tx) dx, composite Gauss-Legendre
    T, W = np.polynomial.legendre.leggauss(20)
    br = np.linspace(0, 2 * b, m // 20 + 1)
    xs = np.concatenate([(c - a) / 2 * T + (c + a) / 2 for a, c in zip(br[:-1], br[1:])])
    ws = np.concatenate([(c - a) / 2 * W for a, c in zip(br[:-1], br[1:])])
    F = SHAPES[shape](xs / b) * sech(w * xs)
    out = np.empty(len(ts))
    for i in range(0, len(ts), 200):
        tt = np.asarray(ts[i:i + 200])[:, None]
        out[i:i + 200] = 2 * (np.cos(tt * xs[None, :]) * (F * ws)[None, :]).sum(axis=1)
    return out

def check_Q(shape, b, w, NZ):
    lA, lB, E = consts(shape, b, w)
    gam = lA  # n = r1 = 1: r1 logA + 2 r2 logB = logA
    # primes
    X = 2 * b; P = 0.0
    for p in range(2, int(math.exp(X)) + 2):
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            m = 1
            while m * math.log(p) < X:
                x = m * math.log(p)
                P += 2 * math.log(p) * p ** (-m / 2) * float(SHAPES[shape](np.array([x / b]))[0] * sech(np.array([w * x]))[0])
                m += 1
    t0 = time.time()
    mp.mp.dps = 20
    zs = [float(mp.zetazero(k).imag) for k in range(1, NZ + 1)]
    fh = fhat(shape, b, w, zs)
    Z = 2 * float(np.sum(fh))          # each pair rho, conj(rho) contributes 2 Fhat(gamma)
    # tail: fit |Fhat| ~ C t^-k on the last zeros, density log(t/2pi)/(2pi)
    tail_terms = np.abs(fh[-50:]) * np.array(zs[-50:]) ** 4
    C = float(np.max(tail_terms)); T = zs[-1]
    tail = 2 * C * (math.log(T / (2 * math.pi)) + 1) / (2 * math.pi) / (3 * T ** 3)
    resid = gam - E + Z + P
    print(f'(1) K=Q {shape} b={b} w={w}: Gamma {gam:.10f}  -E {-E:.10f}  zeros({NZ}, gamma<={T:.1f}) {Z:.10f}  primes {P:.10f}')
    print(f'    identity residual {resid:.3e}  (zero-sum tail bound ~{tail:.1e}; zetazero time {time.time()-t0:.1f}s)')
    return resid, tail

def check_gamma(shape, b, w):
    lA, lB, E = consts(shape, b, w)
    Ic = EG + L8P - lB; Ir = lB + math.pi / 2 - lA
    T, W = np.polynomial.legendre.leggauss(30)
    br = np.concatenate([np.linspace(0, 50, 201), np.geomspace(50, 1500, 160)[1:]])
    ts = np.concatenate([(c - a) / 2 * T + (c + a) / 2 for a, c in zip(br[:-1], br[1:])])
    ws = np.concatenate([(c - a) / 2 * W for a, c in zip(br[:-1], br[1:])])
    fh = fhat(shape, b, w, ts, m=int(40 * 1500 * 2 * b / (2 * math.pi)) + 2000)  # >= 40 nodes per period at t = 1500
    mp.mp.dps = 15
    psi_c = np.array([float(mp.re(mp.digamma(mp.mpc(0.5, t)))) for t in ts])
    psi_r = np.array([float(mp.re(mp.digamma(mp.mpc(0.25, t / 2)))) for t in ts])
    lhs_c = 2 / math.pi * float(np.sum(ws * fh * (math.log(2 * math.pi) - psi_c)))  # even integrand: int_R = 2 int_0^inf
    lhs_r = 1 / math.pi * float(np.sum(ws * fh * (math.log(math.pi) - psi_r)))
    print(f'(2) {shape} b={b} w={w}: Gamma_C check {lhs_c:.9f} vs {2*(EG+L8P-Ic):.9f} (diff {lhs_c-2*(EG+L8P-Ic):.1e});'
          f' Gamma_R check {lhs_r:.9f} vs {EG+math.pi/2+L8P-Ic-Ir:.9f} (diff {lhs_r-(EG+math.pi/2+L8P-Ic-Ir):.1e})')

def check_odlyzko(csv='data/odlyzko1976_tables34.csv'):
    rows = {3: [], 4: []}
    for line in open(csv):
        if line[0].isdigit():
            t, b, A, B, E = line.strip().split(',')
            rows[int(t)].append((float(b), float(A), float(B), float(E)))
    for t, shape, w in ((3, 'G', 0.0), (4, 'H', 0.5)):
        mx, mn, cnt, low_e = 0, 1e9, 0, []
        for b, A, B, E in rows[t]:
            lA, lB, Eo = consts(shape, b, w)
            if math.exp(lA) < A - 1e-9 or math.exp(lB) < B - 1e-9:
                print('   MISMATCH A/B', t, b, A, B, math.exp(lA), math.exp(lB))
            if Eo > E * (1 + 1e-12):
                low_e.append((b, E, Eo))
            mx = max(mx, math.exp(lA) / A - 1, math.exp(lB) / B - 1); mn = min(mn, math.exp(lA) / A - 1, math.exp(lB) / B - 1); cnt += 1
        print(f'(3) Odlyzko 1976 Table {t}: {cnt} rows; ours/his - 1 for A, B in [{mn:.2e}, {mx:.2e}] (his A, B are lower estimates)')
        for b, E, Eo in low_e:
            print(f'    his E at b = {b}: {E:.5e} < exact {Eo:.5e} (his own closed form gives the exact value; rounding/typo in the table)')

def s_family(n, r1):
    """Stark's single-s inequality (s > 1): log|d| >= r1 (log pi - psi(s/2)) + 2 r2 (log 2pi - psi(s)) - 2/s - 2/(s-1)."""
    best = (-1e9, None)
    for s in np.exp(np.linspace(math.log(1.0005), math.log(4.0), 3000)):
        v = (r1 * (math.log(math.pi) - float(mp.digamma(s / 2))) + (n - r1) * (math.log(2 * math.pi) - float(mp.digamma(s))) - 2 / s - 2 / (s - 1)) / n
        if v > best[0]: best = (v, s)
    return math.exp(best[0]), best[1]

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
    NZ = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    if mode == 's_only':
        for n in (2, 8, 100, 1000, 10000, 10 ** 6):
            a, sa = s_family(n, 0); c, sc = s_family(n, n)
            print(f's-family (Stark, s > 1): n={n}: all signatures rd > {a:.3f} (s={sa:.4f}); totally real rd > {c:.3f} (s={sc:.4f})')
        print('s-family limits: 2 pi e^gamma =', 2 * math.pi * math.exp(EG), ' 4 pi e^gamma =', 4 * math.pi * math.exp(EG))
        sys.exit()
    if mode == 'gamma_only':
        check_gamma('H', 3.0, 15 / 32); check_gamma('G', 1.5, 0.3733819666132076); sys.exit()
    check_odlyzko()
    check_gamma('H', 3.0, 15 / 32)
    check_gamma('G', 1.5, 0.3733819666132076)
    check_Q('H', 2.0, 15 / 32, NZ)
    check_Q('G', 1.6, 0.0, NZ)
    for n in (8, 100, 10000, 10 ** 6):
        a, sa = s_family(n, 0); c, sc = s_family(n, n)
        print(f'(4) s-family (Stark, single s > 1): n={n}: all signatures rd > {a:.3f} (s={sa:.4f}); totally real rd > {c:.3f} (s={sc:.4f})')
