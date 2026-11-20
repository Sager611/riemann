#!/usr/bin/env python3
"""cert.py -- lidens certificate (node N3, records only): lambda_n > 0 for all n <= N with Kadiri-Lumley-Ng's zero density
inserted into Li200k's loss sum.  Ball arithmetic (python-flint/Arb); every comparison is a certified sign.

Bound (REPORT.md sec. 2):  lambda_n >= M_n - Q_n,  Q_n <= Qall_n(a*) + sum_j S_j(n) + R_n, where
  Qall_n(a*) = Li200k's Q_hi with a = a* = s* - 1/2 (every zero above T0 charged psi_n(gamma, a*));
  S_j(n)     = f_j(t_K) U_j(t_K) + sum_k U_j(t_{k+1}) (f_j(t_k) - f_j(t_{k+1})),  f_j = psi_n(., a_{j+1}) - psi_n(., a_j),
               U_j(t) = min(2 KLN(s_j, t), Nplus(t)) >= #{T0 < gamma <= t : |beta - 1/2| > a_j},  a_j = s_j - 1/2,
               s_0 = s* < s_1 < ... < s_J = B (dyadic), grid T0 = t_0 < ... < t_K = 2 T0 (integers);
  R_n        = 2 psi_n(t_K, a_B) E(t_K) + int_{t_K}^inf psi_n(t, a_B) w(t) dt  (zeros above t_K charged at a_B, total count).
  KLN(s, t)  = (t-H) log t/(2 pi d) log(1 + 2 pi d A(s) (log kt)^{2s} (log t)^{4(1-s)} t^{8(1-s)/3}/(t-H)) + B(s)(log t)^2
               (KLN Thm 1.1 / Lemma 4.14, H = H0 - 1, H0 = 3e12; '--plain': the simplified form A(s)(log kt)^{2s}(log t)^{5-4s}t^{8(1-s)/3} + B(s)(log t)^2;
               '--nomin': U_j = 2 KLN only).
  Nplus(t)   = Ntilde(t) - Ntilde(T0) + E(t) + E(T0) >= N(t) - N(T0-)  (Hasanalizade-Shen-Wong, E as in Li200k).
Usage: cert.py [--B 141|144] [--sstar 0.83 --J 190 | --mesh 0.80:0.002,0.84:0.0001] [--fine] [--plain] [--nomin] [--delta1]
               [--test c1,c2] [--falsify] [--small]   (all runs: code/runs.sh all)
"""
import sys, os, time, importlib.util, random
from fractions import Fraction
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel)); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m
L = load('li200k', 'Riemann/Deep/Li200k/verify.py')      # read-only: E, Ntilde, Q_hi, tail_bound, M_lo, certify
C = load('ing4cert', 'Riemann/Deep/Ingham4/cert.py')     # read-only: KLN rows (class KP, build_kln_H0)
from flint import arb, acb, ctx
ctx.prec = 160
t_start = time.time()
def opt(k, d):
    return sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
PI = arb.pi(); T0 = L.T0; T0i = 3000175332800
def dyadic_up(fr, bits=64):                    # exact dyadic >= the Fraction fr
    m = -((-fr.numerator * 2 ** bits) // fr.denominator); return Fraction(m, 2 ** bits)
def A_(fr): return arb(fr.numerator) / fr.denominator if fr.denominator & (fr.denominator - 1) == 0 else None
BSRC = opt('--B', '141')
if BSRC == '141':   # theta_141 = 0.873345856106718576801184796... (QRH decimal, PR 141/143); ...797 is an upper bound
    Bfr = dyadic_up(Fraction('0.873345856106718576801184797'))
else:               # QuasiRiemannTracker PR 144 (pending): 436666397/500000000
    Bfr = dyadic_up(Fraction(436666397, 500000000))
SS = Fraction(opt('--sstar', '0.83')); J = int(opt('--J', '190'))
Sfr = dyadic_up(SS, 24)
sig_fr = [Sfr + (Bfr - Sfr) * j / J for j in range(J)] + [Bfr]
if '--mesh' in sys.argv:    # non-uniform slabs: --mesh s*:h1,s2:h2  (step h1 from s* to s2, then h2 up to B)
    parts = [tuple(Fraction(v) for v in x.split(':')) for x in opt('--mesh', '').split(',')]
    Sfr = dyadic_up(parts[0][0], 24); sig_fr = []; s = Sfr
    for i, (s0, h) in enumerate(parts):
        stop = parts[i + 1][0] if i + 1 < len(parts) else Bfr
        while s < stop - h / 1000: sig_fr.append(s); s = s + h
        s = stop if i + 1 < len(parts) else s
    sig_fr.append(Bfr); J = len(sig_fr) - 1
sig_fr = [dyadic_up(s, 40) if j < J else s for j, s in enumerate(sig_fr)]     # exact dyadics, increasing
assert all(sig_fr[j] < sig_fr[j + 1] for j in range(J)) and sig_fr[0] == Sfr
SIG = [A_(s) for s in sig_fr]; AJ = [s - arb(1) / 2 for s in SIG]
aB = AJ[-1]; astar = AJ[0]
assert astar > 0          # a_0 >= 0 is needed for the reflection count and for f_j >= 0, decreasing (reviewer)
PLAIN = '--plain' in sys.argv; NOMIN = '--nomin' in sys.argv
H = arb(3 * 10 ** 12) - 1
DELTA1 = '--delta1' in sys.argv     # robustness: KLN Thm 1.1 as printed (delta >= 1), delta = 1 in every row
ROWS = [p for p in C.build_kln_H0('3e12', '1' if DELTA1 else None) if p.name in ('0.85', '0.86', '0.87', '0.88')]
ctx.prec = 160
for p in ROWS: assert p.ok and p.k == 1 and p.smin < SIG[0]
def kln(p, s, t):
    A = p.logA(s).exp(); Bc = p.Bc(s); lt = t.log(); lkt = (p.k * t).log()
    if PLAIN:
        return A * lkt ** (2 * s) * lt ** (5 - 4 * s) * t ** (arb(8) / 3 * (1 - s)) + Bc * lt ** 2
    v = t - H
    Y = 2 * PI * p.d * A * lkt ** (2 * s) * lt ** (4 * (1 - s)) * t ** (arb(8) / 3 * (1 - s))
    return v * lt / (2 * PI * p.d) * (Y / v).log1p() + Bc * lt ** 2
def Nplus(t): return L.Ntilde(t) - L.Ntilde(T0) + L.E(t) + L.E(T0)
# t-grid (exact integers): u = t - T0 geometric 1.05 from 1e8 to 0.02 T0, then t x1.01 to 1.3 T0, then x1.03 to 2 T0
FINE = '--fine' in sys.argv          # finer grid: u x1.02, t x1.005
ru, rt = (102, 1005) if FINE else (105, 101)
G = [T0i]; u = 10 ** 8
while u < T0i // 50: G.append(T0i + u); u = -(-u * ru // 100)
t = G[-1]
while t < 13 * T0i // 10: t = -(-t * rt // (1000 if FINE else 100)); G.append(t)
while t < 2 * T0i: t = min(-(-t * 103 // 100), 2 * T0i); G.append(t)
TG = [arb(x) for x in G]; K = len(G) - 1
# count bounds U[j][k] (n-independent), stored as exact upper endpoints
U = []
for j in range(J):
    s = SIG[j]
    p = min(ROWS, key=lambda p: float(kln(p, s, T0)))
    row = []
    for k in range(K + 1):
        if k == 0: row.append(arb(0)); continue
        u2 = arb((2 * kln(p, s, TG[k])).upper())             # exact, >= 2 KLN(s_j, t_k)
        if not NOMIN:
            np_ = arb(Nplus(TG[k]).upper())                   # exact, >= Nplus(t_k)
            if np_ < u2: u2 = np_                             # exact comparison of exact numbers
        row.append(u2)
    U.append(row)
L.VARIANTS['astar'] = L.a_const(astar); L.VARIANTS['aB'] = L.a_const(aB)
def psi(n, t, a):
    return 4 * (n * a / (2 * t * t)).sinh() ** 2 + 4 * n * a * a / t ** 3
def Q(m, parts=False):
    """upper bound for sum_{gamma > T0} psi_m(gamma, |beta - 1/2|) (m a Python int)."""
    n = arb(m)
    qa = L.Q_hi(m, 'astar')
    P = [[psi(n, TG[k], a) for k in range(K + 1)] for a in AJ]
    S = arb(0)
    for j in range(J):
        Pj, Pj1, Uj = P[j], P[j + 1], U[j]
        f = [Pj1[k] - Pj[k] for k in range(K + 1)]
        S += f[K] * Uj[K]
        for k in range(K): S += Uj[k + 1] * (f[k] - f[k + 1])
    R = 2 * psi(n, TG[K], aB) * L.E(TG[K]) + L.tail_bound(n, aB, TG[K])
    tot = qa + S + R
    return (tot, qa, S, R) if parts else tot

def falsify(trials=400, seed=7):
    """per-zero inequality at the real scale: rho = 1/2 + d + i y, y in [T0, 3 T0], 0 < d <= a_B, n <= 120 T0^2 (prec 480)."""
    old = ctx.prec; ctx.prec = 480; rnd = random.Random(seed); worst = None
    for i in range(trials):
        y = T0 * (1 + 2 * arb(rnd.random())) if i % 4 else T0
        d = aB * (arb(1) - arb(rnd.random()) ** 3) if i % 3 else aB
        n = int(float(T0) ** 2 * rnd.uniform(1, 120)) + 1
        rho = acb(arb(1) / 2 + d, y)
        z = 1 - 1 / rho; zs = 1 - 1 / (1 - rho.conjugate())
        lhs = (2 - z ** n - zs ** n).real
        th = 2 * (1 / (2 * y)).atan()
        rhs = 2 - 2 * (n * th).cos() - psi(arb(n), y, d)
        sl = lhs - rhs
        if not (sl > -arb(10) ** -60 * (1 + abs(rhs))):
            raise SystemExit('FAIL per-zero inequality: y=%s d=%s n=%d slack=%s' % (y.str(15), d.str(15), n, sl.str(10)))
        rel = float(sl / (abs(rhs) + 1)); worst = rel if worst is None or rel < worst else worst
    # monotonicity spot check: psi decreasing in t, f = psi(a') - psi(a) decreasing in t (analytic proof in REPORT)
    for i in range(200):
        n = arb(int(float(T0) ** 2 * rnd.uniform(1, 120))); t1 = T0 * (1 + arb(rnd.random())); t2 = t1 * (1 + arb(rnd.random()) / 100)
        a1 = aB * arb(rnd.random()); a2 = a1 + (aB - a1) * arb(rnd.random())
        assert psi(n, t1, a2) - psi(n, t2, a2) > 0
        assert (psi(n, t1, a2) - psi(n, t1, a1)) - (psi(n, t2, a2) - psi(n, t2, a1)) > 0
    ctx.prec = old
    print('falsify: per-zero inequality holds in %d real-scale cases (min relative slack %.3g); 200 monotonicity checks pass' % (trials, worst), flush=True)

def chain(n, Qf):
    steps = 0; evals = 0
    while True:
        adv = False
        for rr in L.RATIOS:
            m = int(n * rr)
            if m <= n: continue
            evals += 1
            if L.M_lo(n, m) - Qf(m) > 0:
                n = m; steps += 1; adv = True; break
        if not adv: return n, steps, evals

if __name__ == '__main__':
    T0sq = T0 * T0
    print('lidens cert: B=%s (dyadic %s), s*=%s, J=%d slabs, K=%d t-steps to %s, plain=%s, nomin=%s, delta1=%s; prec %d'
          % (BSRC, float(Bfr), float(Sfr), J, K, G[-1], PLAIN, NOMIN, DELTA1, ctx.prec), flush=True)
    print('rows: ' + ' '.join('%s:A(s0)=%s,B(s0)=%s' % (p.name, p.logA(p.s0).exp().str(6, radius=False), p.Bc(p.s0).str(6, radius=False)) for p in ROWS))
    print('count bound at T0+1e10: U_0=%s U_J-1=%s (Nplus %s); setup %.1f s' % (U[0][G.index(min(G, key=lambda x: abs(x - T0i - 10**10)))].str(5),
          U[J-1][G.index(min(G, key=lambda x: abs(x - T0i - 10**10)))].str(5), Nplus(arb(T0i + 10**10)).str(5), time.time() - t_start), flush=True)
    if '--falsify' in sys.argv: falsify()
    if '--small' in sys.argv: L.part_small()      # n <= 88 via the on-line zero near 14.1347 (Li200k part 1, all variants incl. a*, a_B)
    if '--test' in sys.argv:
        for c in opt('--test', '110').split(','):
            m = int(float(c) * T0i * T0i); t1 = time.time()
            tot, qa, S, R = Q(m, True); Mm = L.M_lo(m, m)
            print('test c=%s: M_lo=%s Q=%s (Qall(a*)=%s, dens=%s, tail=%s) gap>0: %s  [%.2f s]' % (c, Mm.str(6), tot.str(6), qa.str(4), S.str(4), R.str(3), bool(Mm - tot > 0), time.time() - t1), flush=True)
        sys.exit(0)
    t1 = time.time()
    NF, st0 = L.certify('aB')
    print('chain 1 (Li200k bound, a = a_B): lambda_n > 0 for 89 <= n <= %d = %.6f T0^2 (%d steps, %.1f s)' % (NF, float(arb(NF) / T0sq), st0, time.time() - t1), flush=True)
    t1 = time.time()
    N, st, ev = chain(NF, Q)
    print('chain 2 (density bound) from %d: lambda_n > 0 for all n <= N = %d  (N/T0^2 = %.6f; %d steps, %d evaluations, %.1f s)'
          % (NF, N, float(arb(N) / T0sq), st, ev, time.time() - t1), flush=True)
    tot, qa, S, R = Q(N, True)
    print('at N: M_lo(N,N) = %s, Q = %s (Qall(a*) %s, density %s, tail %s)' % (L.M_lo(N, N).str(8), tot.str(8), qa.str(5), S.str(5), R.str(3)))
    print('total runtime %.1f s' % (time.time() - t_start))
