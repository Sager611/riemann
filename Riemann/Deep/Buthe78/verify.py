#!/usr/bin/env python3
"""
verify.py -- rigorous certificate (ball arithmetic, python-flint/arb) for the enlarged range of
Buthe's partial-RH bounds (Buthe 2016, arXiv 1410.7015, Theorem 2; Johnston 2022, arXiv 2109.02249).

Inputs (all cited in REPORT.md):
  * RH for 0 < Im rho <= T = 3e12 (Platt-Trudgian, arXiv 2004.09765);
  * no zero with Re s > 7/8 (OpenAI, Oct 2026; Lean-checked in Riemann/Deep/QRH)   [modes M0, M2 only];
  * Buthe 2016: explicit formula (Prop 1 + 0.024 bound => Theta(2)), Prop 3 (tail, (e:rem1),(e:rem2)),
    Lemma 3 (sums of 1/gamma), Prop 4 (|psi - psi_{c,eps}| for 0 <= alpha < 1);
  * psi(x)-theta(x) < x^{1/2}+x^{1/3}+(log x/log 2) x^{1/4}  (theta(t) < t for t <= 1.39e17, Platt-Trudgian 2016).

For every y in [2.169e25, X'] (log-subintervals, adaptive bisection) it encloses, divided by sqrt(y) log y,
  bound(y) = alpha*eps*e^{alpha eps} y + sqrt(x) S(c) + Z1(x) + 2 + A(x,c,eps,alpha),   x in y e^{[-alpha eps, alpha eps]},
with c = theta0*log y + D, eps = c/T (so c/eps = T), and checks  bound(y) < sqrt(y) log y (log y - C)/(8 pi).
Modes:  B0 = Buthe's own split (e:rem2)+(e:zsum3-bound), alpha = 0, tail with Re rho < 1 only (re-optimized Buthe/Johnston);
        M0 = Buthe's own split (e:rem2)+(e:zsum3-bound) with a^2 = 2k/c, alpha = 0, tail with Re rho <= 7/8;
        M1 = step zero sum + alpha > 0, tail with Re rho < 1 only (no 7/8);
        M2 = step zero sum + alpha > 0, tail with Re rho <= 7/8   (headline).
Statements: psi (C=0), theta (C = C_theta), Pi (C=2), pi (C = 2 + C_theta).
Below 2.169e25 everything is Buthe's published theorem with T = 3e12.

Runtime: 5-10 seconds single-threaded (measured 4.1 s and 7.6 s on an Apple laptop, python-flint 0.9.0); the last line prints it.
Run:  python3 -m venv v && v/bin/pip install python-flint && v/bin/python verify.py   (optional filter: B0 M0 M1 M2)
"""
import sys, time
from flint import arb, ctx

ctx.prec = 160
t0 = time.time()
PI = arb.pi()
T = arb(3 * 10**12)
TI = 3 * 10**12
LOG2 = arb(2).log()

def iv(lo, hi):
    return arb(lo).union(arb(hi))

# ---------- step zero sum (Buthe 2016 Lemma 3 (i),(ii); ell decreasing on [0, c/eps]) ----------
T1 = 10**6           # >= 5000 for Lemma 3 (ii), >= 14 for Lemma 3 (i)
NBLK = 400
grid = [T1]
r = (TI / T1) ** (1.0 / NBLK)
for i in range(1, NBLK):
    g = int(T1 * r**i)
    if g > grid[-1] and g < TI: grid.append(g)
grid.append(TI + 1)  # last block [t_n, T+1) contains all zeros with t_n <= gamma <= T
LG = [(arb(t) / (2 * PI)).log() for t in grid]
BLK = []
for i in range(len(grid) - 1):
    ti = arb(grid[i])
    BLK.append((LG[i + 1]**2 - LG[i]**2) / (4 * PI) + 5 * ti.log() / ti)
S_FIRST = LG[0]**2 / (4 * PI)
SV = [(1 - (arb(grid[i]) / T)**2).sqrt() for i in range(len(grid) - 1)]   # v_i = sqrt(1 - (t_i/T)^2), t_i < T

def S_step(c):
    """upper bound for sum_{0<|gamma|<=T} ell(gamma)/|gamma| (both signs); valid for all c' >= c (ell decreasing in c)."""
    acc = S_FIRST
    e2c = (-2 * c).exp()
    for v, b in zip(SV, BLK):
        # ell(t) = sinh(c v)/(v sinh c) = e^{c(v-1)} (1 - e^{-2cv}) / (v (1 - e^{-2c}))
        ell = (c * (v - 1)).exp() * (1 - (-2 * c * v).exp()) / (v * (1 - e2c))
        acc += ell * b
    return 2 * acc

# ---------- nu_c(alpha), mu_c(alpha) by the exact series of the Bessel kernel ----------
NSER = 220
def JK(al):
    al = arb(al); w = 1 - al * al
    J = 1 - al; out = []
    for n in range(NSER + 1):
        if n > 0:
            J = (2 * n * J - al * w**n) / (2 * n + 1)
        K = w**(n + 1) / (2 * (n + 1))
        out.append((J, K - al * J))
    return out

def nu_mu(c, jk):
    q = c * c / 4
    term = arb(1); smu = arb(0); snu = arb(0)
    for n in range(NSER + 1):
        if n > 0: term = term * q / (n * n)
        J, KJ = jk[n]
        smu += term * J; snu += term * KJ
    # tail n > NSER: term_n <= term_NSER * rho^{n-NSER}, rho = q/(NSER+1)^2 <= 1/2; J, K-alJ <= 1
    assert (q / (NSER + 1)**2 < arb(1) / 2)
    tail = arb(0, (2 * term).upper())
    f = c / (2 * c.sinh())
    return f * (snu + tail), f * (smu + tail)   # |nu|, mu   (upper-enclosed tails)

# ---------- normalized margin on an L-interval ----------
def margin(L, mode, D, al, k, C):
    th0 = arb(3) / 8 if mode in ('M0', 'M2') else arb(1) / 2
    th = th0 + arb(1) / 2
    D = arb(D); al = arb(al)
    c = th0 * L + D
    eps = c / T
    sy = (L / 2).exp()                    # sqrt(y)
    epsy = eps * sy                       # eps*sqrt(y)
    dmax = (al * eps).upper()
    delta = arb(0, dmax) if mode not in ('M0', 'B0') else arb(0)
    # side conditions
    assert eps.upper() < arb('1e-4') and c.lower() > 3 and L.lower() > 5
    # tail |gamma| > c/eps = T  (Buthe Lemma 2 + pairing rho, 1 - conj(rho); Re rho in [1-th, th])
    z1 = arb('0.16') * (2 * (-D).exp() / (1 - (-2 * c).exp())) * ((th * delta).exp() + (-2 * th0 * L + (1 - th) * delta).exp()) \
         * (arb('0.71') * (c * eps).sqrt()).exp() * (3 * c).log() * T.log() / L
    # zeros 0 < |gamma| <= T
    if mode in ('M0', 'B0'):
        k = arb(k); a = (2 * k / c).sqrt()
        assert (a * T).lower() > 5000 and a.upper() < 1
        u = (c * c - 2 * k * c).sqrt()
        coshr = (-2 * k * c / (u + c)).exp() * (1 + (-2 * u).exp()) / (1 - (-2 * c).exp())   # cosh(u)/sinh(c)
        z2 = (1 + 11 * c * eps) / (PI * 2 * k) * T.log() * coshr / L
        z3 = (a * T / (2 * PI)).log()**2 / (2 * PI) / L
        zs = (delta / 2).exp() * (z2 + z3)
    else:
        zs = (delta / 2).exp() * S_step(c.lower()) / L
    # |psi(x) - psi_{c,eps}(x)| <= A(x,c,eps,alpha)  (Buthe Prop 4)
    jk = JK_CACHE[(str(al.mid()), mode)]
    nu, mu = nu_mu(c, jk)
    if mode in ('M0', 'B0'): mu = arb(1) / 2
    logB = eps.log() + L + delta - eps + (nu / (2 * mu)).log()
    assert logB.lower() > 0
    pre = (2 * eps).exp() * (L + delta + eps) / L
    a1 = pre * 2 * epsy * delta.exp() * nu / logB
    a2 = pre * (arb('2.01') * eps * (delta / 2).exp() + (LOG2 + 2 * (L + delta)).log() / 2 / sy)
    shift = al * epsy * (al * eps).exp() / L
    const = 2 / (sy * L)
    target = (L - C(L)) / (8 * PI)
    return target - (z1 + zs + a1 + a2 + shift + const), dict(z1=z1, zs=zs, a1=a1, shift=shift)

JK_CACHE = {}

def C_theta(L):
    return 8 * PI * (1 + (-L / 6).exp() + (L / LOG2) * (-L / 4).exp()) / L

CS = {'psi': lambda L: arb(0), 'theta': C_theta, 'Pi': lambda L: arb(2), 'pi': lambda L: 2 + C_theta(L)}

def certify(mode, stmt, D, al, k, Xp, minw=1e-7, w0=0.02):
    key = (str(arb(al).mid()), mode)
    if key not in JK_CACHE: JK_CACHE[key] = JK(al)
    Lst = float((arb(2169) * arb(10)**22).log().lower()) - 1e-9
    Len = float((arb(Xp)).log().upper()) + 1e-9
    stack = []
    n = int((Len - Lst) / w0) + 1
    for i in range(n):
        lo = Lst + (Len - Lst) * i / n; hi = Lst + (Len - Lst) * (i + 1) / n
        stack.append((lo, hi))
    stack.reverse()
    cnt = 0; worst = None
    while stack:
        lo, hi = stack.pop()
        m, _ = margin(iv(lo, hi), mode, D, al, k, CS[stmt])
        cnt += 1
        if m > 0:
            lb = m.lower()
            if worst is None or lb < worst[0]: worst = (lb, lo, hi)
            continue
        if hi - lo < minw:
            return False, cnt, (lo, hi, m)
        mid = (lo + hi) / 2
        stack.append((mid, hi)); stack.append((lo, mid))
    return True, cnt, worst

CASES = [
    # mode, statement, D, alpha, k, X'
    ('B0', 'psi',   '4.5182', '0', '0.8819', '1.99e26'),
    ('B0', 'theta', '4.5596', '0', '0.8819', '1.86e26'),
    ('B0', 'Pi',    '4.7391', '0', '0.8820', '1.40e26'),
    ('B0', 'pi',    '4.7938', '0', '0.8820', '1.28e26'),
    ('M0', 'psi',   '4.1951', '0', '0.8791', '2.60e26'),
    ('M0', 'theta', '4.2365', '0', '0.8791', '2.43e26'),
    ('M0', 'Pi',    '4.4172', '0', '0.8791', '1.82e26'),
    ('M0', 'pi',    '4.4721', '0', '0.8791', '1.67e26'),
    ('M1', 'psi',   '4.3125', '0.1019', None, '3.98e26'),
    ('M1', 'theta', '4.3438', '0.1022', None, '3.77e26'),
    ('M1', 'Pi',    '4.4813', '0.1025', None, '3.01e26'),
    ('M1', 'pi',    '4.5219', '0.1025', None, '2.82e26'),
    ('M2', 'psi',   '3.9875', '0.1156', None, '5.19e26'),
    ('M2', 'theta', '4.0187', '0.1159', None, '4.92e26'),
    ('M2', 'Pi',    '4.1594', '0.1163', None, '3.92e26'),
    ('M2', 'pi',    '4.2000', '0.1163', None, '3.68e26'),
]

def consts_5000():
    # |Pi(5000)-li(5000)-(psi(5000)-5000)/log 5000| and |pi(5000)-li(5000)-(theta(5000)-5000)/log 5000| < sqrt(5000)/(4 pi)
    N = 5000
    sieve = bytearray([1]) * (N + 1); sieve[0] = sieve[1] = 0
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]: sieve[p*p::p] = bytearray(len(sieve[p*p::p]))
    primes = [p for p in range(2, N + 1) if sieve[p]]
    th = sum((arb(p).log() for p in primes), arb(0)); ps = arb(0); Pi = arb(0)
    for p in primes:
        pk, m = p, 1
        while pk <= N:
            ps += arb(p).log(); Pi += arb(1) / m; pk *= p; m += 1
    li = arb(N).li(); lg = arb(N).log(); bound = arb(N).sqrt() / (4 * PI)
    cPi = abs(Pi - li - (ps - N) / lg); cpi = abs(len(primes) - li - (th - N) / lg)
    return cPi, cpi, bound

if __name__ == '__main__':
    cPi, cpi, bnd = consts_5000()
    print('5000-constants: |Pi..| =', cPi.str(6), ' |pi..| =', cpi.str(6), ' must be < sqrt(5000)/(4pi) =', bnd.str(6),
          'OK' if (cPi < bnd and cpi < bnd) else 'FAIL')
    allok = (cPi < bnd and cpi < bnd)
    only = sys.argv[1:]  # optional filter, e.g. M2
    for mode, stmt, D, al, k, Xp in CASES:
        if only and mode not in only: continue
        ok, cnt, info = certify(mode, stmt, D, al, k, Xp)
        allok &= ok
        if ok:
            print('%s %-5s x in [2.169e25, %s]: CERTIFIED  (D=%s alpha=%s%s; %d boxes; min margin lower bound %.3e at L=%.5f)'
                  % (mode, stmt, Xp, D, al, (' k=' + k) if k else '', cnt, float(info[0]), info[1]))
        else:
            print('%s %-5s up to %s: FAILED at L in [%.8f, %.8f], margin %s' % (mode, stmt, Xp, info[0], info[1], info[2].str(5)))
        sys.stdout.flush()
    print('ALL OK' if allok else 'SOME FAILED', ' runtime %.1f s' % (time.time() - t0))
