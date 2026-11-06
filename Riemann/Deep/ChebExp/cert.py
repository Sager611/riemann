#!/usr/bin/env python3
"""Rigorous certificate (python-flint / Arb ball arithmetic) for the explicit polylogarithmic
least-prime bound in every Chebotarev class.  See DERIVATION.md for every inequality used.

  uv run --no-project --with python-flint python3 cert.py [chains.json] > logs/cert_output.txt

Checked here, for every chain in chains.json (one chain = one constant c on one range Lambda >= Lambda_0,
Lambda = log Delta_L):
  * per interval [La, Lb] with a fixed kernel (h, s*, s0) and x = (c (La + log 3))^(1/delta):
      LHS(x; Lb, n) > 0 for n = 2 and n = Lb / log(pi/2)   (DERIVATION (5.1), affine in n, decreasing in Lambda)
  * per chain end [LN, oo) with a fixed kernel: the monotone tail bound (DERIVATION (5.3)) at LN;
  * the intervals are contiguous;
  * every constant (Phi(1), Phi(0), ||k||, A1, A2, prime-power sums, digamma) is an Arb enclosure and
    is used in the safe direction (lower bounds for the main term, upper bounds for errors).
Also checked: the direct ranges. For delta = 1/32: Delta_L <= 12 (only quadratic fields, by Minkowski), computed
here. For delta = 1 - B_2: Delta_L <= 257 (degree <= 4 by Minkowski), read from logs/direct_range.txt, which is
produced by direct.gp (PARI nflist). Run from this directory.
"""
import json, math, sys, time
from flint import arb, acb, ctx

ctx.prec = 120
PI = arb.pi()
LOG2 = arb(2).log()
LOG3 = arb(3).log()
LP2 = (PI / 2).log()                     # log Delta_L >= n_L log(pi/2)  (Minkowski, n_L >= 2)
B2 = (777 - arb(1081).sqrt()) / 852      # Liu's half-plane (Hexagon 2611.00082)
DELTAS = {'1/32': arb(1) / 32, '1-(777-1081**0.5)/852': 1 - B2}
I = acb(0, 1)
XEXACT = 2.0e12                          # exact prime-power sum when x <= XEXACT
THETA_C = arb('1.01624')                 # theta(y) < 1.01624 y for y > 0 (Rosser-Schoenfeld 1962, Thm 9)

def up(a):  # upper endpoint as arb (exact)
    return a.upper()

def lo(a):
    return a.lower()

def sieve(n):
    b = bytearray([1]) * (n + 1); b[0] = b[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if b[i]: b[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if b[i]]

PRIMES = sieve(int(XEXACT ** 0.5) + 10)
LOGP = None

# ---------------------------------------------------------------- kernel
class Kernel:
    """F(u) = e^{ss u} B3(u; h) on [0, 3h], B3 = 1_[0,h]^{*3};  k(t) = F(log 1/t);
       Phi(s) = int_0^1 k(t) t^{s-1} dt = ((1 - e^{-h(s-ss)})/(s-ss))^3."""
    def __init__(self, h, ss, s0):
        self.hf, self.ssf, self.s0f = h, ss, s0
        self.h, self.ss, self.s0 = arb(h), arb(ss), arb(s0)
        self.U = 3 * self.h
        assert h > 0 and s0 > 1

    def Phi(self, s):
        z = s - self.ss
        w = self.h * z / 2
        return (self.h * (-w).exp() * (I * w).sinc()) ** 3   # (1-e^{-hz})/z = h e^{-w} sinh(w)/w

    @staticmethod
    def M3piece(j, v):
        # cardinal quadratic B-spline (Irwin-Hall density, 3 summands) on piece [j, j+1]
        if j == 0: return v * v / 2
        if j == 1: return (-2 * v * v + 6 * v - 3) / 2
        return (3 - v) * (3 - v) / 2

    def F_up(self, u):
        """rigorous upper bound (an exact arb) for F on the ball u; F = e^{ss u} h^2 M3(u/h) >= 0,
        evaluated with the polynomial of each knot piece the ball meets (no cancellation)."""
        v = u / self.h
        vl, vu = float(lo(v).mid()), float(up(v).mid())
        best = arb(0)
        for j in range(3):
            if vu < j - 1e-15 or vl > j + 1 + 1e-15:
                continue
            val = up((self.ss * u).exp() * self.h * self.h * self.M3piece(j, v))
            if val > best: best = val
        return best

    def consts(self):
        self.P1 = lo(self.Phi(acb(1)).real)          # main term: lower bound
        self.P0 = up(self.Phi(acb(0)).real)          # Phi(0) = int F: upper bound
        # ||k|| = max F on [0, 3h]: cover by N balls
        self.kn = self.maxF(arb(0))                  # ||k|| = max F
        self.S2 = self.maxF(arb(1) / 2)              # max_u F(u) e^{-u/2}
        self.S3 = self.maxF(arb(1) / 3)
        self.Ph2 = up(self.Phi(acb(arb(1) / 2)).real)  # Phi(1/2) = int F e^{-u/2}
        self.Ph3 = up(self.Phi(acb(arb(1) / 3)).real)
        return self

    def maxF(self, shift):
        """rigorous upper bound for max_u e^{-shift u} F(u) = max_u e^{(ss-shift) u} B3(u)"""
        N = 2000; best = arb(0); a = self.ss - shift
        for j in range(3):                            # piece [j h, (j+1) h], N closed sub-balls each
            for i in range(N):
                v = j + arb(2 * i + 1, 1) / (2 * N)
                val = up((a * self.h * v).exp() * self.h * self.h * self.M3piece(j, v))
                if val > best: best = val
        return best

    def R(self, beta, gam, sig):   # |Phi(beta + i gam)| / Re(1/(sig - beta - i gam))
        s = acb(beta, gam)
        d = sig - beta
        return abs(self.Phi(s)) * (d * d + gam * gam) / d

    def sup_line(self, beta, tol):
        """rigorous upper bound for sup_{gamma >= 0} R(beta, gamma; s0)"""
        sig = self.s0; d = sig - beta
        e = (self.h * (self.ss - beta)).exp()
        def tail(G):
            G = arb(G)
            return up((1 + e) ** 3 * (d * d + G * G) / (d * G ** 3))
        R0 = lo(self.R(beta, arb(0), sig))
        G = 1.0
        while tail(G) > R0 / 4:
            G *= 2
        # branch and bound on [0, G]
        n0 = 256
        work = [(G * i / n0, G * (i + 1) / n0) for i in range(n0)]
        best_lo = R0; best_up = tail(G); count = 0
        while work:
            nxt = []
            for (a, b) in work:
                count += 1
                ub = up(self.R(beta, arb((a + b) / 2, (b - a) / 2), sig))
                if ub <= best_lo * (1 + tol):
                    if ub > best_up: best_up = ub
                    continue
                pv = lo(self.R(beta, arb((a + b) / 2), sig))
                if pv > best_lo: best_lo = pv
                if b - a < 1e-12:
                    if ub > best_up: best_up = ub
                    continue
                m = (a + b) / 2
                nxt += [(a, m), (m, b)]
            # re-test pending intervals against the improved lower bound happens in the next round
            work = nxt
            if count > 4_000_000:
                raise RuntimeError('sup_line: too many intervals')
        return best_up, best_lo, count

# ---------------------------------------------------------------- arithmetic pieces
def digamma(x):
    return arb(x).digamma()

def g_gamma(s0):
    """lower bound for g(s0) = min(a, b/2): gamma_L'/gamma_L(s0) = -r1 a - r2 b <= -n_L g(s0)."""
    a = (PI.log() - digamma(s0 / 2)) / 2
    b = (2 * PI).log() - digamma(s0)
    return lo(arb.min(a, b / 2))

def pow_sum(K, lx):
    """upper bound for P_k(x) = sum_{p^j, j>=2} log p * k(p^j/x), x = e^lx (lx an arb ball)"""
    xf = math.exp(float(lx.mid()))
    if xf * (1 + 1e-9) <= XEXACT:
        tot = arb(0); j = 2
        Uf = float(K.U.mid()) * (1 + 1e-9) + 1e-9
        lxf = float(lx.mid())
        while 2.0 ** j <= xf * 1.000001:
            pmax = math.exp(lxf / j) * (1 + 1e-9) + 1
            pmin = math.exp((lxf - Uf) / j) * (1 - 1e-9) - 1
            for p in PRIMES:
                if p > pmax: break
                if p < pmin: continue
                lp = arb(p).log()
                u = lx - j * lp
                if lo(u) > K.U or up(u) < 0:
                    continue
                tot += lp * K.F_up(u)
            j += 1
        return up(tot), 'exact'
    # smooth (DERIVATION 4.3): sum_p log p phi(p) <= 1.01624 (sup_y y phi(y) + int phi), phi(y) = k(y^j/x) unimodal:
    # P_k <= 1.01624 sum_{j=2}^{J} x^{1/j} (S_j + Phi(1/j)/j),  S_j <= ||k||, Phi(1/j) <= Phi(0) for j >= 4
    J = math.floor(float(up(lx).mid()) / math.log(2)) + 1
    tot = (lx / 2).exp() * (K.S2 + K.Ph2 / 2) + (lx / 3).exp() * (K.S3 + K.Ph3 / 3) \
        + max(J - 3, 0) * (lx / 4).exp() * (K.kn + K.P0 / 4)
    return up(THETA_C * tot), 'smooth'

def triv_ratio(K, lx, J1=30):
    """upper bound for tau_x := sum_{j>=1} x^{-j} Phi(-j) / Phi(0)   (DERIVATION 3.1');
       Phi(-j) in closed form for j <= J1, then Phi(-j) <= e^{jU} Phi(0) (geometric tail)."""
    tot = arb(0)
    for j in range(1, J1 + 1):
        tot += (-j * lx).exp() * K.Phi(acb(-j)).real
    q = (K.U - lx).exp()
    tot += K.P0 * q ** (J1 + 1) / (1 - q)
    return up(tot / K.P0)

def ram_factor(K, lx):
    """R_f = min(||k|| log x, ||k||(U + log 2), Phi(0) + ||k|| log 2) / log 2   (DERIVATION 4.2)"""
    a = up(K.kn * lx); b = up(K.kn * (K.U + LOG2)); cc = up(K.P0 + K.kn * LOG2)
    m = a if a < b else b
    m = m if m < cc else cc
    return up(m / LOG2)

def lhs_interval(delta, c, La, Lb, K, A1, A2):
    """lower bound of LHS(x; Lambda, n) for Lambda in [La, Lb], all n in [2, Lb/log(pi/2)],
       at x = (c (La + log 3))^(1/delta)."""
    ell_a = arb(La) + LOG3
    lx = (arb(c) * ell_a).log() / delta
    assert lo(lx) > up(K.U)  # x > e^U: trivial-zero series converges
    x = lx.exp()
    Mx = up(arb.max(A1, up(((2 * delta - 1) * lx).exp() * A2)))
    xd = up((-delta * lx).exp())
    eps = triv_ratio(K, lx)
    Rf = ram_factor(K, lx)
    pw, mode = pow_sum(K, lx)
    g = g_gamma(K.s0)
    cB = up(1 / K.s0 + 1 / (K.s0 - 1))
    Lb_ = arb(Lb)
    res = []
    for n in (arb(2), Lb_ / LP2):
        Bb = Lb_ / 2 + cB - n * g
        val = K.P1 - xd * Mx * Bb - ((n - 1) * K.P0 + n * K.P0 * eps + 2 * Lb_ * Rf + n * pw) / x
        res.append(lo(val))
    return min(res[0], res[1]), mode, float(lx.mid()), float(Mx.mid()), float(K.P1.mid())

def lhs_tail(delta, c, LN, K, A1, A2):
    """lower bound, valid for all Lambda >= LN at x = (c(Lambda+log 3))^(1/delta) (DERIVATION (5.3))."""
    ellN = arb(LN) + LOG3
    cc = arb(c)
    assert lo(cc * ellN) > arb(1).exp()
    lx = (cc * ellN).log() / delta
    assert lo(lx) > up(K.U)  # x > e^U: trivial-zero series converges
    x = lx.exp()
    Mx = up(arb.max(A1, up(((2 * delta - 1) * lx).exp() * A2)))
    eps = triv_ratio(K, lx)                       # decreasing in x
    g = g_gamma(K.s0)
    assert g > 0                                  # used as -n g <= -2 g
    cB = up(1 / K.s0 + 1 / (K.s0 - 1))
    add = cB - LOG3 / 2 - 2 * g
    add = up(add) if up(add) > 0 else arb(0)
    nmax = arb(LN) / LP2
    # smooth powers (DERIVATION 4.3): P_k/x <= 1.01624 [x^{-1/2}(S2+Phi(1/2)/2) + x^{-2/3}(S3+Phi(1/3)/3) + J x^{-3/4}(||k||+Phi(0)/4)]
    J = math.floor(float(up(lx).mid()) / math.log(2)) + 1     # J <= log x/log 2 + 1 <= 2 log x for log x >= 2
    powx = THETA_C * ((-lx / 2).exp() * (K.S2 + K.Ph2 / 2) + (-2 * lx / 3).exp() * (K.S3 + K.Ph3 / 3)
                      + 2 * lx * (-3 * lx / 4).exp() * (K.kn + K.P0 / 4))
    # monotone decrease in Lambda (x = (c ell)^{1/delta}) of Lambda/x, Lambda log x/x, Lambda/sqrt x, Lambda x^{-2/3},
    # Lambda log x x^{-3/4}: sufficient conditions (DERIVATION 5.3)
    assert lo((1 - delta) * ellN) > up(LOG3)
    assert up(2 * delta * ellN / (ellN - LOG3)) < 1
    assert up(delta * (ellN / (ellN - LOG3) + 1 / (cc * ellN).log())) < arb(3) / 4
    val = K.P1 - Mx * (1 / (2 * cc) + add / (cc * ellN)) \
        - (nmax * K.P0 * (1 + eps) + 2 * K.kn * arb(LN) * lx / LOG2) / x - nmax * up(powx)
    return lo(val), float(lx.mid()), float(Mx.mid()), float(K.P1.mid())

# ---------------------------------------------------------------- direct range Delta_L <= 12
def kron(d, p):
    # Kronecker symbol (d/p) for p prime
    if p == 2:
        if d % 2 == 0: return 0
        return 1 if d % 8 in (1, 7) else -1
    r = pow(d % p, (p - 1) // 2, p)
    return 0 if d % p == 0 else (1 if r == 1 else -1)

def small_fields():
    # all number fields with 3 <= Delta_L <= 12 are quadratic (Minkowski: cubic >= 12.49, quartic >= 43.3, ...)
    ds = [-3, -4, 5, -7, 8, -8, -11, 12]
    out = []
    for d in ds:
        split = next(p for p in PRIMES if kron(d, p) == 1)
        inert = next(p for p in PRIMES if kron(d, p) == -1)
        ram = next(p for p in PRIMES if kron(d, p) == 0)
        # K = Q: classes {1} (split), {sigma} (inert); K = L: least degree-one prime of L = min(split, ram)
        out.append((d, abs(d), split, inert, min(split, ram)))
    return out

EXPECTED = {'C2': 157, 'C3': 3, 'S3': 32, 'C4': 1, 'V4': 3, 'D4': 2, 'A4': 0, 'S4': 2}

def direct_file(path='logs/direct_range.txt'):
    """rows (group, degree, Delta_L, description, least p) from direct.gp (PARI): every field with Delta_L <= 257,
       every subfield K with L/K Galois, every class.  Consistency checks: field counts per Galois group and the
       number of classes per field."""
    rows = [l.split() for l in open(path) if l.strip()]
    fields = {}
    for r in rows:
        key = (r[0], int(r[2]), r[1])
        fields.setdefault(key, []).append(r)
    # several fields can share (group, Delta): count K=L lines instead
    cnt = {}
    for r in rows:
        if 'K=L' in r: cnt[r[0]] = cnt.get(r[0], 0) + 1
    for g, n in EXPECTED.items():
        assert cnt.get(g, 0) == n, ('field count', g, cnt.get(g, 0), n)
    for (g, D, deg), rs in fields.items():
        nL = sum(1 for r in rs if 'K=L' in r)
        nQ = sum(1 for r in rs if r[3].startswith('K=Q'))
        nq = sum(1 for r in rs if r[3].startswith('K=quad'))
        if g in ('C2', 'C3', 'C4', 'V4'): assert nQ == nL * int(deg), ('classes over Q', g, D)
        if g in ('C4', 'D4'): assert nq == 2 * nL, ('quadratic subfield classes', g, D)
        if g == 'V4': assert nq == 6 * nL, ('V4 subfields', D)
        assert all(int(r[-1]) > 1 for r in rs)
    return [(r[0], int(r[2]), ' '.join(r[3:-1]), int(r[-1])) for r in rows]

def main():
    t0 = time.time()
    path = sys.argv[1] if len(sys.argv) > 1 else 'chains.json'
    chains = json.load(open(path))
    print('B2 =', B2.str(25), ' delta_Liu =', (1 - B2).str(25), ' 1/delta_Liu =', (1 / (1 - B2)).str(20))
    print('log(pi/2) =', LP2.str(12))
    # ---- direct range
    sf = small_fields()
    print('\nDirect range (all fields with Delta_L <= 12 are quadratic):')
    print('   d  Delta  least split  least inert  least deg-1 prime of L')
    for row in sf: print('  %3d  %5d  %11d  %11d  %d' % row)
    for dexpr, cdirect in [('1/32', None), ('1-(777-1081**0.5)/852', None)]:
        delta = DELTAS[dexpr]
        # c_min = max over fields of N^delta / log(3 Delta)
        cmin = max(((arb(max(r[2], r[3], r[4])).log() * delta).exp() / (arb(3 * r[1])).log()) for r in sf)
        print(f'  delta={dexpr}: smallest c allowed by the direct range = {up(cmin).str(8)} (attained by Q(sqrt(-3)), N = 7)')
    allok = True
    summary = []
    for ch in chains:
        name, dexpr, c, L0 = ch['name'], ch['delta_expr'], ch['c'], ch['Lam0']
        delta = DELTAS[dexpr]
        print(f'\n=== chain {name}: delta = {dexpr} = {delta.str(15)}, c = {c}, Lambda >= {L0} ===')
        if ch.get('uniform'):
            D0 = ch['direct']
            # covers every L != Q: analytic chain from log(D0 + 1), direct range Delta_L <= D0
            assert arb(L0) <= arb(D0 + 1).log(), 'chain must start at or below log(D0+1)'
            if D0 == 12:
                okd = all(arb(max(r[2], r[3], r[4])) <= (arb(c) * (arb(3 * r[1])).log()) ** (1 / delta) for r in sf)
                worst = max(((arb(max(r[2], r[3], r[4])).log() * delta).exp() / (arb(3 * r[1])).log()) for r in sf)
                nchk = 3 * len(sf)
            else:
                assert D0 == 257
                rows = direct_file()
                okd = all(arb(r[3]) <= (arb(c) * (arb(3 * r[1])).log()) ** (1 / delta) for r in rows)
                worst = max(((arb(r[3]).log() * delta).exp() / (arb(3 * r[1])).log()) for r in rows)
                nchk = len(rows)
            print(f'  uniform chain: analytic part starts at {L0} <= log({D0 + 1}); direct range Delta_L <= {D0} '
                  f'({nchk} (L, K, C) cases): all N <= (c log 3Delta)^(1/delta): {okd}; largest N^delta/log(3Delta) = {up(worst).str(8)}')
            allok = allok and okd
        ivs = ch['intervals']
        prev_b = L0; ok = True; worst = None
        for k, iv in enumerate(ivs):
            assert iv['Lam_a'] <= prev_b + 1e-15, ('gap', name, k)
            K = Kernel(iv['h'], iv['ss'], iv['s0']).consts()
            A1, A1lo, n1 = K.sup_line(1 - delta, 1e-3)
            A2, A2lo, n2 = K.sup_line(delta, 0.5)
            if iv['kind'] == 'int':
                v, mode, lx, Mx, P1 = lhs_interval(delta, c, iv['Lam_a'], iv['Lam_b'], K, A1, A2)
                rel = v / K.P1
                good = v > 0
                print(f'  [{iv["Lam_a"]:.6f}, {iv["Lam_b"]:.6f}] h={iv["h"]:.5f} s*={iv["ss"]:.5f} s0={iv["s0"]:.5f} log x={lx:.3f} '
                      f'P1={P1:.5g} A1={float(A1.mid()):.5g} (lo {float(A1lo.mid()):.5g}, {n1} boxes) M={Mx:.5g} pow:{mode} '
                      f'LHS/P1 >= {float(rel.mid()):.5f} {"OK" if good else "FAIL"}', flush=True)
                prev_b = iv['Lam_b']
            else:
                v, lx, Mx, P1 = lhs_tail(delta, c, iv['Lam_a'], K, A1, A2)
                rel = v / K.P1
                good = v > 0
                print(f'  tail [{iv["Lam_a"]:.6f}, oo) h={iv["h"]:.5f} s*={iv["ss"]:.5f} s0={iv["s0"]:.5f} log x>={lx:.3f} '
                      f'P1={P1:.5g} A1={float(A1.mid()):.5g} M={Mx:.5g} LHS/P1 >= {float(rel.mid()):.5f} {"OK" if good else "FAIL"}', flush=True)
                prev_b = float('inf')
            ok = ok and good
            worst = rel if worst is None or rel < worst else worst
        ok = ok and prev_b == float('inf')
        allok = allok and ok
        summary.append((name, dexpr, c, L0, ok, float(worst.mid())))
        print(f'  chain {name}: {"CERTIFIED" if ok else "NOT CERTIFIED"} (min LHS/P1 = {float(worst.mid()):.5f})')
    print('\nSUMMARY')
    for s in summary:
        print('  %-8s delta=%-22s c=%-7s Lambda>=%-8s %s  min margin %.5f' % (s[0], s[1], s[2], s[3], 'CERTIFIED' if s[4] else 'FAILED', s[5]))
    print('ALL CERTIFIED' if allok else 'SOME FAILED', f'({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
