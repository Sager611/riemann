"""Arb certificate for Bui-Milinovich's large-gap criterion (arXiv 1410.3635v3) with an optimised amplifier.

Reduction (BM l.164-202, l.321-347): Lambda > kappa if h1 = c0/(kappa^2 c1) > 1, where (Thm smoothed, l.300-313; r = 1)
  c_j = int_{V} e^{i Phi} a b (th d12 + s)(th d34 + s) (v - S)^{2j} P(w1) P(w2),
  V = {x, x1..x4, t1..t4 in [0,1]: x+x1+x2 <= 1, x+x3+x4 <= 1}, a = 1-th(x1+x3), b = 1-th(x2+x4),
  d12 = x1-x2, d34 = x3-x4, s = a t1 - b t2, S = th(x+x1+x2+x3+x4) + a t1 + b t2,
  w1 = 1-x-x1-x2, w2 = 1-x-x3-x4, Phi = th k pi(x2-x4-d34 t3+d12 t4) - k pi s (t3-t4).
Step 1 (exact): int dt3 dt4 e^{i Phi} = e^{i th k pi (x2-x4)} E(w3) E(w4), w3 = -k pi (th d34 + s), w4 = k pi (th d12 + s),
  E(w) = (e^{iw}-1)/(iw); times (th d12 + s)(th d34 + s) this is e^{i th k pi(x2-x4)}(e^{i w3}-1)(e^{i w4}-1)/(k pi)^2.
  Four terms: A = e^{i th k pi (x1-x3)}, B = e^{i th k pi (x2-x4)}, C = -e^{-i k pi (t1-t2)} e^{i l^C.x}, D = -e^{i k pi (t1-t2)} e^{i l^D.x},
  l^C = k pi th (t1, 1-t2, -(1-t1), -t2), l^D = k pi th (1-t1, t2, -t1, -(1-t2)) on (x1, x2, x3, x4).
Step 2: A, B: t-integrals exact (polynomial). C, D: Gauss-Legendre in (t1, t2) with a rigorous Bernstein-ellipse remainder.
Step 3 (x-integrals, for fixed real l): expand e^{i(l1 u1 + l2 u2)} to total order N (remainder |z|^{N+1}/(N+1)!, |z| <= k pi th),
  Dirichlet: int_{u1+u2<=R} u1^al u2^be (R-u1-u2)^c = al! be! c!/(al+be+c+2)! R^{al+be+c+2}; then int_0^1 x^e (1-x)^m dx exactly.
All arithmetic in Arb/acb balls (outward rounding) from exact rational inputs.
Run: uv run --no-project --with python-flint python3 cert_bm.py"""
import sys, time, resource
from math import comb, factorial
from flint import fmpq, fmpq_mpoly_ctx, arb, acb, acb_poly, ctx

ctx.prec = 192
PI = arb.pi()
A_ = lambda q: arb(q.p) / arb(q.q)
F = [fmpq(factorial(i)) for i in range(400)]


def shleg_mono(n):
    return [fmpq((-1) ** (n + k) * comb(n, k) * comb(n + k, k)) for k in range(n + 1)]


def dec(s):
    neg = s.startswith('-'); s = s.lstrip('+-'); a, b = (s.split('.') + [''])[:2]
    q = fmpq(int(a + b), 10 ** len(b)); return -q if neg else q


def leg_to_mono(pl):
    P = [fmpq(0)] * len(pl)
    for m, cm in enumerate(pl):
        for k, cf in enumerate(shleg_mono(m)):
            P[k] += cm * cf
    return P


def sup_abs_bound(P, M=400):
    """rigorous upper bound for max_{[0,1]} |P| (exact rationals): on each [k/M, (k+1)/M] with midpoint m and
    half-width h, |P(m+d)| <= sum_k |P^{(k)}(m)|/k! h^k (finite Taylor expansion, exact); also <= max |Bernstein coeff|."""
    d = len(P) - 1
    bern = max(abs(sum((fmpq(comb(k, j), comb(d, j)) * P[j] for j in range(k + 1)), fmpq(0))) for k in range(d + 1))
    ders = [P]
    for k in range(d):
        prev = ders[-1]
        ders.append([prev[i] * i for i in range(1, len(prev))])
    h = fmpq(1, 2 * M)
    best = fmpq(0)
    for k in range(M):
        m = fmpq(2 * k + 1, 2 * M)
        s = fmpq(0)
        for o, Dk in enumerate(ders):
            val = fmpq(0)
            for c in reversed(Dk):
                val = val * m + c
            s += abs(val) / F[o] * h ** o
        if s > best:
            best = s
    return min(best, bern)


class Cert:
    def __init__(self, th, kap, v, P, N=40, nt=14, label=''):
        self.th, self.kap, self.v, self.P, self.N, self.nt, self.label = th, kap, v, P, N, nt, label
        self.kp = A_(kap) * PI                     # kappa*pi
        self.lam = self.kp * A_(th)                # kappa*pi*theta
        self.deg = len(P) - 1
        self.Pnorm = sup_abs_bound(P)
        # polynomial parts Pi_j(x, x1..x4, t1, t2)
        cx = fmpq_mpoly_ctx.get(('x', 'x1', 'x2', 'x3', 'x4', 't1', 't2'), 'lex')
        x, x1, x2, x3, x4, t1, t2 = cx.gens()
        a = 1 - th * (x1 + x3); b = 1 - th * (x2 + x4)
        S = th * (x + x1 + x2 + x3 + x4) + a * t1 + b * t2
        self.Pis = [a * b, a * b * (v - S) ** 2]
        self.split = []
        for Pi in self.Pis:
            d = {}
            for mon, cf in zip(Pi.monoms(), Pi.coeffs()):
                d.setdefault(tuple(int(z) for z in mon[:5]), []).append(((int(mon[5]), int(mon[6])), fmpq(cf)))
            self.split.append(d)
        # PiL[L][c] = P_c c!/(L+c+2)!  (exact -> arb)
        self.amax = max(max(k[1] + k[2], k[3] + k[4]) for d in self.split for k in d)
        self.Lmax = self.amax + N
        self.PiL = [[A_(P[c] * F[c] / F[L + c + 2]) for c in range(self.deg + 1)] for L in range(self.Lmax + 1)]
        self.jdeg = self.Lmax + 2 + self.deg
        self.Beta = {e: [A_(F[e] * F[m] / F[e + m + 1]) for m in range(2 * self.jdeg + 1)] for e in range(3)}
        self.facs = [A_(F[i]) for i in range(N + self.amax + 2)]
        self.invf = [A_(1 / F[i]) for i in range(N + 1)]
        # truncation error of exp-Taylor per group: eps = lam^{N+1}/(N+1)!
        self.eps = self.lam ** (N + 1) / A_(F[N + 1])

    def Jpolys(self, l1, l2):
        """dict (al, be) -> list of acb coefficients in R of int_{u1+u2<=R} u1^al u2^be P(R-u1-u2) e^{i(l1u1+l2u2)} (truncated)."""
        N = self.N
        I = acb(0, 1)
        p1 = [acb(1)]; p2 = [acb(1)]
        z1 = I * l1; z2 = I * l2
        for m in range(1, N + 1):
            p1.append(p1[-1] * z1); p2.append(p2[-1] * z2)
        p1 = [p1[m] * self.invf[m] for m in range(N + 1)]
        p2 = [p2[m] * self.invf[m] for m in range(N + 1)]
        out = {}
        for al in range(self.amax + 1):
            f1 = acb_poly([p1[m] * self.facs[al + m] for m in range(N + 1)])
            for be in range(self.amax + 1 - al):
                f2 = acb_poly([p2[m] * self.facs[be + m] for m in range(N + 1)])
                sK = f1 * f2
                coef = [acb(0)] * (self.jdeg + 1)
                for K in range(N + 1):
                    sk = sK[K]
                    L = al + be + K
                    for c in range(self.deg + 1):
                        coef[L + 2 + c] += sk * self.PiL[L][c]
                out[(al, be)] = coef
        return out

    def Xint(self, Qs, l):
        """[int_{x-polytope} Q_j P(w1) P(w2) e^{i l.x} for each Q_j in Qs]; Q_j: dict (e,a1,a2,a3,a4) -> acb."""
        J1 = self.Jpolys(l[0], l[1]); J2 = self.Jpolys(l[2], l[3])
        res = []
        cache = {}
        for Q in Qs:
            tot = acb(0); absq = arb(0)
            groups = {}
            for (e, a1, a2, a3, a4), q in Q.items():
                groups.setdefault((a1, a2, a3, a4), []).append((e, q))
                absq += abs(q)
            for (a1, a2, a3, a4), eqs in groups.items():
                key = (a1, a2, a3, a4)
                if key not in cache:
                    cache[key] = acb_poly(J1[(a1, a2)]) * acb_poly(J2[(a3, a4)])
                pr = cache[key]
                for e, q in eqs:
                    Be = self.Beta[e]
                    s = acb(0)
                    for m in range(pr.degree() + 1):
                        s += pr[m] * Be[m]
                    tot += q * s
            # truncation remainder: |e^{iz1}e^{iz2} - T1T2| <= 2eps+eps^2, |Q| <= sum|q|, int |P P| <= ||P||^2/20
            err = (2 * self.eps + self.eps ** 2) * absq * A_(self.Pnorm) ** 2 / 20
            res.append(tot + acb(arb(0, err), arb(0, err)))
        return res

    def run(self):
        t0 = time.time()
        th, kp, lam = self.th, self.kp, self.lam
        # terms A, B: exact t-integration of Pi_j
        Qbar = []
        for d in self.split:
            Q = {}
            for xm, lst in d.items():
                Q[xm] = acb(sum((A_(cf / ((m + 1) * (k + 1))) for (m, k), cf in lst), arb(0)))
            Qbar.append(Q)
        XA = self.Xint(Qbar, [lam, arb(0), -lam, arb(0)])
        XB = self.Xint(Qbar, [arb(0), lam, arb(0), -lam])
        # terms C, D: Gauss-Legendre in (t1, t2)
        nodes = []
        for k in range(self.nt):
            r_, w_ = arb.legendre_p_root(self.nt, k, weight=True)
            nodes.append(((1 + r_) / 2, w_ / 2))
        SC = [acb(0), acb(0)]; SD = [acb(0), acb(0)]
        for (u1, w1) in nodes:
            for (u2, w2) in nodes:
                Qs = []
                for d in self.split:
                    Q = {}
                    for xm, lst in d.items():
                        Q[xm] = acb(sum((A_(cf) * u1 ** m * u2 ** k for (m, k), cf in lst), arb(0)))
                    Qs.append(Q)
                XC = self.Xint(Qs, [lam * u1, lam * (1 - u2), -lam * (1 - u1), -lam * u2])
                XD = self.Xint(Qs, [lam * (1 - u1), lam * u2, -lam * u1, -lam * (1 - u2)])
                ph = acb(0, 1) * kp * (u1 - u2)
                eC = (-ph).exp(); eD = ph.exp()
                for j in range(2):
                    SC[j] += w1 * w2 * eC * XC[j]
                    SD[j] += w1 * w2 * eD * XD[j]
        # rigorous remainder of the (t1,t2) tensor Gauss rule, per term and per variable (Bernstein ellipse rho):
        # |E| <= 4 M rho^{-(2nt-1)}/(rho-1), M <= e^{k pi b_rho} * sup|Pi_j| * ||P||^2/20, b_rho = (rho-1/rho)/4,
        # sup|Pi_0| <= 1, sup|Pi_1| <= (v + 2 th + tmax + 1)^2, tmax = 1/2 + (rho+1/rho)/4 (|t| on the ellipse).
        best = None
        for rho10 in range(15, 400):
            rho = arb(rho10) / 10
            brho = (rho - 1 / rho) / 4; tmax = arb(1) / 2 + (rho + 1 / rho) / 4
            Smax = A_(abs(self.v)) + 2 * A_(th) + tmax + 1
            base = 4 * (kp * brho).exp() * A_(self.Pnorm) ** 2 / 20 * rho ** (-(2 * self.nt - 1)) / (rho - 1)
            eb = [2 * 2 * base, 2 * 2 * base * Smax ** 2]          # 2 variables x 2 terms (C, D)
            if best is None or eb[1].upper() < best[1][1].upper():
                best = (rho, eb)
        rho, eb = best
        c = []
        for j in range(2):
            Ej = arb(0, eb[j].upper())
            cj = (XA[j] + XB[j] - SC[j] - SD[j] + acb(Ej, Ej)) / kp ** 2
            c.append(cj)
        self.rho, self.eb, self.c = rho, eb, c
        self.time = time.time() - t0
        return c

    def report(self, need=True):
        c0, c1 = self.c
        k2 = A_(self.kap) ** 2
        H = c0.real - k2 * c1.real
        h1 = c0.real / (k2 * c1.real)
        ok_imag = c0.imag.contains(0) and c1.imag.contains(0)
        ok_pos = c0.real > 0 and c1.real > 0
        ok_H = H > 0
        print('[%s] theta=%s kappa=%s v=%s deg P=%d N=%d nt=%d rho=%s' % (self.label, self.th, self.kap, self.v, self.deg, self.N, self.nt, self.rho.str(3)))
        print('   ||P||_inf <= %s ; Taylor eps = %s ; t-Gauss remainder (c0, c1 before /(k pi)^2) <= %s, %s' % (
            self.Pnorm, self.eps.str(3), self.eb[0].str(3), self.eb[1].str(3)))
        print('   c0 = %s' % c0.str(12))
        print('   c1 = %s' % c1.str(12))
        print('   h1 = c0/(kappa^2 c1) in %s ; H = c0 - kappa^2 c1 in %s' % (h1.str(12), H.str(8)))
        print('   checks: Im c0, Im c1 contain 0: %s ; c0, c1 > 0: %s ; H > 0: %s   (%.0fs, peak RSS %.0f MB)' % (
            ok_imag, ok_pos, ok_H, self.time, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20), flush=True)
        return ok_imag and ok_pos and (ok_H or not need), h1


# ---- certificate data (exact decimals; P in the shifted-Legendre basis psi_m on [0,1]) ----
BM_P = ['1', '-5.8', '6.4']                       # Bui-Milinovich's P (monomial), l.333
REC_THETA, REC_KAPPA, REC_V = '0.2499', '3.2648', '1.310630'
REC_PLEG = ['1', '2.295083', '5.711135', '2.475600', '13.800476', '-2.126714', '25.669320', '-11.008222',
            '38.228923', '-19.964301', '45.417070', '-22.792755', '40.630752', '-14.895337', '21.791921']   # degree 14


def full():
    ok_all = True
    print('Bui-Milinovich criterion (1410.3635v3): Lambda > kappa if h1 = c0/(kappa^2 c1) > 1 (l.327-330); r = 1.')
    print('BM P (monomial, l.333) =', ' '.join(BM_P))
    # (a) their published numerics: theta = 1/4 (boundary value they used), kappa = 3.18, v = 1.26, deg-2 P: h1 > 1.0002 (l.331-337)
    C = Cert(fmpq(1, 4), dec('3.18'), dec('1.26'), [dec(s) for s in BM_P], N=30, nt=20, label='BM l.333 numerics, theta=1/4')
    C.run(); ok, h1 = C.report()
    ok = ok and h1 > A_(dec('1.0002'))
    print('   reproduces BM l.337 (h1 > 1.0002): %s' % ok); ok_all &= bool(ok)
    # (b) the same P with theta = 2499/10000 < 1/4 (inside Thm smoothed): certified Lambda > 3.18
    C = Cert(dec('0.2499'), dec('3.18'), dec('1.26'), [dec(s) for s in BM_P], N=30, nt=20, label='BM P, theta=2499/10000')
    C.run(); ok, h1 = C.report(); ok_all &= bool(ok)
    # (c) record: optimised P
    pl = [dec(s) for s in REC_PLEG]
    P = leg_to_mono(pl)
    print('record P: Legendre coefficients p_m (exact decimals) =', ' '.join(REC_PLEG))
    print('record P: monomial coefficients (exact) =', [str(c) for c in P])
    C = Cert(dec(REC_THETA), dec(REC_KAPPA), dec(REC_V), P, N=30, nt=20, label='RECORD')
    C.run(); ok, h1 = C.report(); ok_all &= bool(ok)
    assert dec(REC_THETA) < fmpq(1, 4)
    print('theta = %s < 1/4 (BBLR range, l.246, l.259): True' % dec(REC_THETA))
    if ok_all:
        print('ALL CHECKS PASS: Lambda > %s unconditionally (critical zeros).' % REC_KAPPA)
    else:
        print('SOME CHECK FAILED')
    return ok_all


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'full'
    if mode == 'full':
        sys.exit(0 if full() else 1)
    if mode == 'quick':
        C = Cert(fmpq(1, 4), fmpq(318, 100), fmpq(126, 100), [fmpq(1), dec('-5.8'), dec('6.4')], N=int(sys.argv[2]), nt=int(sys.argv[3]), label='BM P quick')
        C.run(); C.report(need=False)
    if mode == 'try':
        # try.py-style: th kap v N nt p0 p1 ... (Legendre coefficients)
        th, kap, v = dec(sys.argv[2]), dec(sys.argv[3]), dec(sys.argv[4]); N, nt = int(sys.argv[5]), int(sys.argv[6])
        pl = [dec(s) for s in sys.argv[7:]]
        C = Cert(th, kap, v, leg_to_mono(pl), N=N, nt=nt, label='try')
        C.run(); C.report(need=False)
