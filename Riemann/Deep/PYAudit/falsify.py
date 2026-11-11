# PYAudit falsification tests that complement round R438 (which tested Lemma 1.4 k=3,4,5 on log/polynomial phases,
# Lemma 2.3 on log blocks and Lemma 2.2). Here: (weyl) Lemma 1.2 + the r-sum inequality; (lem21) Lemma 2.1 incl. the alpha=1/2 corner;
# (phi5) Lemma 1.4 k=5 on the actual phase phi_r of Lemma 3.4 in Yang's main case N >= 2337; (block) the certificate's per-block
# inequality (display l. 826-831 + repairs 3a/3b, constants from SubWeyl/certify.py) against exact block sums at small t.
# uv run --no-project --with mpmath --with numpy --with python-flint python3 falsify.py <test>
import sys, math, os, time
import numpy as np
rng = np.random.default_rng(20261111)
TWO_PI = 2*math.pi
def e(x): return np.exp(2j*np.pi*x)
out = []
def log(s):
    print(s); sys.stdout.flush()

def test_weyl():
    worst = 0.0; n = 0
    for trial in range(400):
        N = int(rng.integers(2, 1500)); a = int(rng.integers(0, 10**6))
        kind = trial % 4
        x = np.arange(a + 1, a + N + 1, dtype=np.float64)
        if kind == 0: f = lambda y: (rng_t/TWO_PI)*np.log(y)
        if kind == 1: f = lambda y: c3*y**3 + pq*y**2
        if kind == 2: f = lambda y: pq*y**2
        if kind == 3: f = lambda y: np.sqrt(y)*c3*1e6
        rng_t = float(rng.uniform(1e3, 1e7)); c3 = float(10**rng.uniform(-12, -6)); pq = float(rng.integers(1, 13))/float(rng.integers(1, 13))
        z = e(f(x)); S = abs(z.sum())
        for q in sorted(set([1, 2, 3, int(rng.integers(1, 3*N)), int(rng.integers(1, N + 1)), N, N + 5])):
            acc = 0.0
            for r in range(1, q):
                if r < N: acc += (1 - r/q)*abs((z[r:]*np.conj(z[:-r])).sum())
            rhs = (N - 1 + q)*(N/q + 2/q*acc)
            worst = max(worst, S*S/rhs); n += 1
    log('weyl: %d (f,N,q) cases, max |S|^2/RHS = %.12f (must be <= 1)' % (n, worst))
    # r-sum inequality sum_{r<q}(1-r/q) r^s <= q^{1+s}/((1+s)(2+s)) for the exponents used
    ss = [11/30, 61/120, -1/2, 1/3, 0.0, -19/30, -11/30] + [1/(2*J - 2) for J in (4, 8, 16)] + [-1/(2*J - 2) for J in (4, 8, 16)] + [1.0, -0.99]
    wr = 0.0
    for s in ss:
        for q in list(range(1, 3000)) + [10**4, 10**5]:
            r = np.arange(1, q, dtype=np.float64)
            lhs = ((1 - r/q)*r**s).sum() if q > 1 else 0.0
            wr = max(wr, lhs/(q**(1 + s)/((1 + s)*(2 + s))))
    log('r-sum inequality: %d exponents, q<=1e5: max LHS/RHS = %.12f (must be <= 1)' % (len(ss), wr))

def gl_integral(phase, a, b, width=0.02, order=12):
    xg, wg = np.polynomial.legendre.leggauss(order)
    nP = max(1, int(math.ceil((b - a)/width)))
    edges = np.linspace(a, b, nP + 1)
    mid = (edges[1:] + edges[:-1])/2; hw = (edges[1:] - edges[:-1])/2
    X = (mid[:, None] + hw[:, None]*xg[None, :]).ravel(); W = (hw[:, None]*wg[None, :]).ravel()
    return X, W

def test_lem21():
    from mpmath import mp
    worst = (0, None); n = 0; t0 = time.time()
    for trial in range(240):
        kind = trial % 4
        if kind == 0:   # f = (t/2pi) log x, f' decreasing
            a = float(rng.uniform(5, 400)); L = float(rng.uniform(0.5, 80)); b = a + L
            span = float(rng.uniform(0.05, 25)); t = span*TWO_PI/(1/a - 1/b)
            f = lambda y: (t/TWO_PI)*np.log(y); fp = lambda y: t/(TWO_PI*y)
        elif kind == 1:  # quadratic -c x^2 + mu x
            a = float(rng.uniform(-50, 50)); L = float(rng.uniform(0.5, 80)); b = a + L
            c = float(rng.uniform(0.001, 12))/L; mu = float(rng.uniform(-5, 5))
            f = lambda y: -c*y*y + mu*y; fp = lambda y: -2*c*y + mu
        elif kind == 2:  # cubic with small leading term, f' decreasing on [a,b] (a>0)
            a = float(rng.uniform(1, 60)); L = float(rng.uniform(0.5, 60)); b = a + L
            c = float(10**rng.uniform(-6, -2)); mu = float(rng.uniform(-3, 3))
            f = lambda y: -c*y**3 + mu*y; fp = lambda y: -3*c*y*y + mu
        else:            # alpha pinned to the corners -1/2, +1/2 (exactly and nearby), short and long ranges
            a = float(rng.uniform(0, 30)); L = float(rng.choice([0.7, 1.3, 2.0, 5.0, 17.0, 60.0])); b = a + L
            c = float(rng.uniform(0.0005, 3))/L
            target = float(rng.choice([0.5, -0.5, 0.4999, -0.4999, 0.5001, 1.5, -1.5]))
            mu = target + 2*c*b   # f'(b) = -2cb + mu = target
            f = lambda y: -c*y*y + mu*y; fp = lambda y: -2*c*y + mu
        al, be = fp(b), fp(a)
        if not (be > al): continue
        if be - al > 40: continue
        ns = np.arange(math.floor(a) + 1, math.floor(b) + 1, dtype=np.float64)
        ns = ns[(ns > a) & (ns <= b)]
        S = e(f(ns)).sum()
        ms = [m for m in range(math.floor(al - 0.5) - 1, math.ceil(be + 0.5) + 2) if al - 0.5 < m < be + 0.5]
        X, W = gl_integral(None, a, b, width=min(0.02, 0.25/(1 + abs(be) + abs(al))), order=12)
        base = e(f(X))*W
        I = sum((base*e(-m*X)).sum() for m in ms)
        lhs = abs(S - I); rhs = 3/math.pi*math.log(be - al + 2) + 4
        n += 1
        if lhs/rhs > worst[0]: worst = (lhs/rhs, dict(kind=kind, a=a, b=b, alpha=al, beta=be, lhs=lhs, rhs=rhs, nm=len(ms)))
    log('lem21: %d cases, max LHS/RHS = %.4f at %s  (%.0f s)' % (n, worst[0], worst[1], time.time() - t0))
    # 30-digit recheck of the worst case (kind 1/3 quadratic or 0 log) by mpmath quad on the same panels
    return worst

def phi_r_and_d5(t, r, nu, mp):
    x = (mp.sqrt(r*r + 2*t*r/(mp.pi*nu)) - r)/2
    g = -(t/(2*mp.pi))*mp.log(1 + r/x)
    ph = g - nu*x
    X = mp.pi*r*nu/t
    d5 = 3/(2*mp.sqrt(mp.pi))*(8*X**3 + 36*X**2 + 60*X + 35)/(X + 2)**mp.mpf(3.5)*mp.sqrt(t*r)/nu**mp.mpf(4.5)
    return ph, d5

LT0, LT1 = 60.0, 95.0   # log(t r): Yang's main case N >= 2337, lambda_5 >= lambda_0(5) for every eta in the grid
def test_phi5():
    from mpmath import mp
    mp.dps = 40
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SubWeyl'))
    sys.dont_write_bytecode = True   # keep SubWeyl/ free of __pycache__
    import certify as C
    from flint import arb
    # closed form of phi^(5) checked against numerical differentiation once
    t, r, nu = mp.mpf(10)**22, mp.mpf(7), mp.mpf(123456)
    num = mp.diff(lambda v: phi_r_and_d5(t, r, v, mp)[0], nu, 5)
    log('phi5 closed form check: |phi^(5)| num %s vs closed %s' % (mp.nstr(abs(num), 15), mp.nstr(phi_r_and_d5(t, r, nu, mp)[1], 15)))
    etas = [0.01*1.25**k for k in range(30)]
    worst = (0, None); n = 0; t0 = time.time()
    for trial in range(36):
        N = int(rng.choice([2400, 3000, 5000, 8000, 12000]))
        H = float(rng.uniform(1.005, 1.1))
        A = int(N/(H*H - 1))
        lt = float(rng.uniform(LT0, LT1)); tr = mp.e**lt
        r = int(rng.integers(1, 40)); t = tr/r
        nus = [mp.mpf(A + k) for k in range(1, N + 1)]
        tot = mp.mpc(0)
        for v in nus:
            ph, _ = phi_r_and_d5(t, r, v, mp)
            tot += mp.expjpi(2*(ph - mp.floor(ph)))
        S = abs(tot)
        d5a = phi_r_and_d5(t, r, mp.mpf(A + 1), mp)[1]; d5b = phi_r_and_d5(t, r, mp.mpf(A + N), mp)[1]
        lam = min(d5a, d5b); h = max(d5a, d5b)/lam
        best = None
        for eta in etas:
            A5, B5 = C.AkBk(5, arb(eta), arb(float(h)))
            A5 = float(A5.mid()); B5 = float(B5.mid())
            bnd = A5*float(h)**(1/8)*N*float(lam)**(1/30) + B5*N**(7/8)*float(lam)**(-1/30)
            best = bnd if best is None else min(best, bnd)
        # Yang's main case flag: N >= 2337 and lambda >= lambda_0(5)
        main = [N >= 2337 and float(lam) >= (9*math.pi/1024*eta)**(-2 + 2/16)*N**(-4 + 4/16) for eta in (etas[0], etas[-1])]
        n += 1
        ratio = float(S)/best
        if ratio > worst[0]: worst = (ratio, dict(N=N, A=A, logtr=lt, r=r, lam=float(lam), h=float(h), S=float(S), bound=best, main=main))
    log('phi5: %d cases (k=5 on phi_r, min over %d etas in [0.01, 6.5]), max |S|/bound = %.5f at %s  (%.0f s)' % (n, len(etas), worst[0], worst[1], time.time() - t0))

def test_block():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SubWeyl'))
    sys.dont_write_bytecode = True   # keep SubWeyl/ free of __pycache__
    import certify as C
    from flint import arb
    worst = (0, None); n = 0; t0 = time.time()
    rcache = {}
    for trial in range(160):
        t = float(10**rng.uniform(6, 9.5))
        amax = math.sqrt(t/TWO_PI)
        a = int(10**rng.uniform(math.log10(50), math.log10(amax)))
        H = float(rng.choice([1.004, 1.012, 1.03, 1.0513, 1.2, 1.6, 1.99]))
        if H not in rcache:
            rcache[H] = C.floats_of(C.region_consts(arb(H)))
        rf = rcache[H]
        b = int(math.floor(H*a))
        if b <= a: continue
        ns = np.arange(a + 1, b + 1, dtype=np.float64)
        z = np.exp(-1j*((t % (TWO_PI*1e6))*np.log(ns) + (t - t % (TWO_PI*1e6))*np.log(ns) % TWO_PI))
        z = np.exp(-1j*np.mod(t*np.log(ns), TWO_PI))
        Smax2 = float(np.max(np.abs(np.cumsum(z)))**2)
        E6 = float(C.E6.mid())
        best = None
        for eta in list(rf['aba'].keys())[::3]:
            d = rf['aba'][eta]
            for q in sorted(set([1, 2, 3, 5, 8, 13, 21, 34, 55, max(1, int(1.0*a**(36/41)*t**(-11/41))), max(1, int(3*a**(36/41)*t**(-11/41)))])):
                F1 = (H - 1) + q/a
                F2 = (H - 1)*a*a/q + d['c2']*a**0.8*(q*t)**(11/30) + d['c3']*a**0.45*(q*t)**(61/120) + rf['c4']*a**2.5*(t*q)**-0.5 \
                    + (rf['c5a']*a**(-1/3) + rf['c5b'])*a**(2/3)*(t*q)**(1/3) + E6*a + d['c7a']*a**2.8*t**(-19/30)*q**(-19/30) + d['c7b']*a**2.2*t**(-11/30)*q**(-11/30)
                v = F1*F2
                best = v if best is None else min(best, v)
        n += 1
        ratio = Smax2/best
        if ratio > worst[0]: worst = (ratio, dict(t=t, a=a, H=H, Smax2=Smax2, rhs=best))
    log('block: %d blocks (display l.826-831 + repairs, min over q and eta), max max_L|S_f(a,L-a)|^2 / RHS = %.5f at %s (%.0f s)' % (n, worst[0], worst[1], time.time() - t0))

def test_poissoncorner():
    # Lemma 2.3 (imp_poisson) at the configuration where the S_3 step of its proof fails: nu_1 - alpha -> 1/2+, beta - nu_N -> 1/2+,
    # with lambda_2 of order 1 (so the constant 6 and 4.686/sqrt(lambda_2) are both tested). Quadratic phases (lambda_3 = 0 allowed:
    # |f'''| = 0 <= h_3 lambda_3 for every lambda_3 > 0, so T_3 -> 0) and cubic phases (exact x_nu).
    worst = (0, None); n = 0
    for trial in range(3000):
        lam = float(10**rng.uniform(-2.5, 0.3))          # |f''| = lam (quadratic)
        M = int(rng.integers(1, 12))                     # number of interior nu
        eps1 = float(10**rng.uniform(-9, -1)); eps2 = float(10**rng.uniform(-9, -1))
        nu1 = int(rng.integers(-20, 20))
        alpha = nu1 - 0.5 - eps1                         # nu1 - alpha = 1/2 + eps1
        beta = nu1 + M - 1 + 0.5 + eps2                  # beta - nu_N = 1/2 + eps2
        if trial % 3 == 1: alpha, beta = nu1 + eps1, nu1 + M + 1 - eps2        # alpha, beta just inside integers
        if trial % 3 == 2: alpha, beta = nu1 - eps1, nu1 + M + eps2           # alpha, beta just outside integers
        L = (beta - alpha)/lam                           # b - a, since f' = beta - lam (x - a)
        a = float(rng.uniform(-100, 100)); b = a + L
        if L > 4e5: continue
        # f(x) = beta (x - a) - lam (x - a)^2 / 2 + phi0 ; f' decreasing from beta to alpha
        phi0 = float(rng.uniform(0, 1))
        f = lambda x: beta*(x - a) - lam*(x - a)**2/2 + phi0
        ns = np.arange(math.floor(a) + 1, math.floor(b) + 1, dtype=np.float64); ns = ns[(ns > a) & (ns <= b)]
        if len(ns) == 0: continue
        # exact phases: reduce (x - a) products carefully in float64 (|values| < 1e7 here)
        S = np.exp(2j*np.pi*np.mod(f(ns), 1.0)).sum()
        nus = np.arange(math.floor(alpha) + 1, math.floor(beta) + 1, dtype=np.float64); nus = nus[(nus > alpha) & (nus <= beta)]
        xn = a + (beta - nus)/lam
        D = (np.exp(2j*np.pi*np.mod(f(xn) - nus*xn - 0.125, 1.0))/math.sqrt(lam)).sum()
        lhs = abs(S - D)
        rhs = 4.686/math.sqrt(lam) + 5/math.pi*math.log(beta - alpha + 2) + 6     # T_3 = 0 for quadratics
        n += 1
        if lhs/rhs > worst[0]: worst = (lhs/rhs, dict(lam=lam, M=M, alpha=alpha, beta=beta, L=L, lhs=lhs, rhs=rhs))
    log('poissoncorner: %d quadratic cases at the S_3-failure configuration, max LHS/RHS = %.4f at %s' % (n, worst[0], worst[1]))

if __name__ == '__main__':
    {'weyl': test_weyl, 'lem21': test_lem21, 'phi5': test_phi5, 'block': test_block, 'poissoncorner': test_poissoncorner}[sys.argv[1]]()
