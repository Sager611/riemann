#!/usr/bin/env python3
"""
tail.py -- the sieve for omega(p-1) > W1 (DERIVATION.md sec 5): Buchstab identity at w0, Selberg upper bounds for
every prime q in [w0, z], analytic bounds for all z >= Z1.  Imported by cert.py (section 'tail').

Weights on divisors of P(z) = prod_{p <= z} p (z = p_w, the shifted sieve of sec 2):
    lambda = mu * 1_{d | P(w0)}  -  sum_{w0 <= q <= z} q * (Selberg Lambda^2 weights for the primes < q, level xi_q),
    xi_q = (z^s / q)^{1/2}.
Main term  T = V(w0) - sum_q 1/(q G_q(xi_q)) = V+(z) - sum_q (V(q)/q)(Phi_q - 1),  Phi_q = 1/(V(q) G_q(xi_q)).
Cost       C - 1 <= 2^{pi(w0)} + c_sig^2 z^s (K_A + K_B).
Mertens inputs (Rosser-Schoenfeld 1962, Thms 4, 5, 7), used only for x >= 286 (resp. x >= 41 for theta):
    R1/R2  |sum_{p<=x} 1/p - loglog x - M| < 1/(2 log^2 x)                 (x >= 286)
    R3/R4  e^-g/log x (1 - 1/(2log^2 x)) < prod_{p<=x}(1-1/p) < e^-g/log x (1 + 1/(2log^2 x))   (x >= 286)
    R5     theta(x) > x (1 - 1/log x)                                           (x >= 41)
"""
import math
from fractions import Fraction as Fr
from flint import arb, ctx

ctx.prec = 120
EG = arb.const_euler().exp()          # e^gamma

def A(x):
    if isinstance(x, Fr):
        return arb(x.numerator) / x.denominator
    return arb(x)

def lo(x):
    return arb(x.lower())

def up(x):
    return arb(x.upper())

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]

def ein_plus(c):
    """sum_{k>=1} c^k/(k k!) = int_0^1 (e^{cv}-1) dv/v  (c > 0), upper bound (series + geometric tail)"""
    c = A(c); tot = arb(0); term = arb(1)
    for k in range(1, 80):
        term = term * c / k
        tot += term / k
    # tail k >= 80: term_k <= c^k/k! decays geometrically with ratio <= c/81 < 1/2 for c <= 40
    tail = term * c / 80 * 2
    return tot + tail

class Tail:
    def __init__(self, s, w0=50, Q2=10000, alphas=None, cs=None, acap=None):
        self.s = A(s); self.sF = Fr(s)
        self.w0, self.Q2 = w0, Q2
        P = primes_upto(max(Q2, 300000))
        self.Pall = P
        self.small = [p for p in P if p <= 283]
        self.PA = [q for q in P if w0 <= q <= Q2]
        self.pi_w0 = sum(1 for q in P if q < w0)
        # V(q) = prod_{p<q}(1-1/p) for q <= Q2 (exact balls)
        self.V = {}
        v = arb(1)
        for p in P:
            if p > Q2:
                break
            self.V[p] = v
            v = v * (1 - arb(1) / p)
        self.alphas = alphas or [Fr(3, 100), Fr(5, 100), Fr(8, 100), Fr(12, 100), Fr(18, 100), Fr(25, 100), Fr(35, 100), Fr(1, 2)]
        # Pi_alpha(q) = prod_{p<q}(1 + (p^alpha - 1)/p)  (exact balls, q <= Q2)
        self.Pia = {}
        for al in self.alphas:
            alA = A(al); pr = arb(1); d = {}
            for p in P:
                if p > Q2:
                    break
                d[p] = pr
                pr = pr * (1 + (arb(p) ** alA - 1) / p)
            self.Pia[al] = d
        self.LQ = arb(Q2).log()
        self.epsM = 1 / (2 * self.LQ ** 2)                     # R3 at q > Q2
        self.etaB = 1 / self.LQ ** 2 + arb(1) / Q2
        self.cs = cs or [Fr(1, 2), Fr(1), Fr(3, 2), Fr(2), Fr(5, 2), Fr(3), Fr(4), Fr(5), Fr(6)]
        L288 = arb(288).log()
        small = [p for p in P if p <= 283]
        self.RB = {}
        for c in self.cs:
            cA = A(c)
            ex = sum(((arb(p) ** (cA / self.LQ)) - 1) / p for p in small)
            self.RB[c] = ex + ein_plus(c) + (cA.exp() - 1) * (1 / (2 * self.LQ ** 2) + 1 / (2 * L288 ** 2))
        # acap: a little below the argmax of h_Q(a) = 2a - a log a - 1 - etaB (a-1) - (log a + etaB)/log Q2
        best = None
        for k in range(1000, 3000):
            a = k / 1000.0
            h = 2 * a - a * math.log(a) - 1 - float(self.etaB.upper()) * (a - 1) - (math.log(a) + float(self.etaB.upper())) / float(self.LQ.lower())
            if best is None or h > best[1]:
                best = (a, h)
        self.acap = acap if acap is not None else Fr(int(best[0] * 1000) - 1, 1000)
        # c_sigma = prod_p (1 + 2/(p(p-1))) <= prod_{p < N} (...) * exp(2/(N-1))
        cs_ = arb(1)
        for p in P:
            cs_ = cs_ * (1 + arb(2) / (p * (p - 1)))
        N = P[-1] + 1
        self.c_sig = cs_ * (arb(2) / (N - 1)).exp()

    # ---------------- Part A: per-prime lower bounds for G_q(xi), q in [w0, Q2]  (sec 5.2)
    def G_lo(self, q, lx):
        """certified lower bound (exact point) for G_q(xi), log xi = lx (arb, lx > 0)"""
        lq = arb(q).log()
        # (a) G(xi) >= G(min(xi, q)) >= log min(xi, q)   (Lemma T3a)
        la, lb = lx.lower(), lq.lower()
        cands_arb = [la if la < lb else lb]
        # (b) H-bound, xi > q   (Lemma T3b): exact sum over primes l in [q, min(xi, 288)), Mertens bound on [288 or q, xi)
        if lx > lq:
            xi = lx.exp()
            ex = arb(0)
            for l in self.small:
                if l < q:
                    continue
                if A(l) > xi:            # l certainly >= xi: stop; undecided l are included (upper bound)
                    break
                ex += (1 + lx - arb(l).log()) / l
            if not (lx < arb(293).log()):    # xi may exceed 293 (no primes in [288, 293))
                x0 = max(q, 288)
                L0 = arb(x0).log()
                a0 = lx / L0
                eta0 = 1 / L0 ** 2 + arb(1) / x0
                mer = a0.log() + eta0 * (1 + lx - L0) + L0 * (a0 * a0.log() - a0 + 1)
            else:
                mer = arb(0)
            gb = lx - ex - mer
            cands_arb.append(lo(gb))
        # (c) Rankin
        V = self.V[q]
        for al in self.alphas:
            r = (-A(al) * lx).exp() * self.Pia[al][q]
            gc = (1 - r) / V
            cands_arb.append(lo(gc))
        best = cands_arb[0]
        for c in cands_arb[1:]:
            if c > best:
                best = c
        return best

    def partA(self, u, uq_max):
        """sum_{w0 <= q <= min(Q2, e^uq_max)} (1/q)(1/G_lo(q, xi_q(u)) - V(q)),  xi_q(u) = exp((s u - log q)/2)"""
        tot = arb(0)
        for q in self.PA:
            if q > self.Q2:
                break
            lq = arb(q).log()
            if lq > uq_max:
                break
            lx = (self.s * u - lq) / 2
            g = self.G_lo(q, lx)
            tot += (1 / g - self.V[q]) / q
        return up(tot)

    def K_A(self, u1):
        tot = arb(0)
        for q in self.PA:
            lq = arb(q).log()
            if lq > u1:
                break
            lx = (self.s * u1 - lq) / 2
            g = self.G_lo(q, lx)
            tot += 1 / (q * g * g)
        # primes q in (e^{u1}, Q2] (if any) enter only for larger z, where G is larger: bound them at their own u = log q
        for q in self.PA:
            lq = arb(q).log()
            if not (lq > u1):
                continue
            lx = (self.s * lq - lq) / 2
            g = self.G_lo(q, lx)
            tot += 1 / (q * g * g)
        return up(tot)

    # ---------------- Part B: uniform bound psi(a) >= Phi_q - 1 for q > Q2  (sec 5.3)
    def psi(self, a):
        """a: exact Fraction (branches decided exactly); returns an exact upper endpoint"""
        aF = Fr(a)
        a = A(aF)
        vals = []
        if aF <= 1:
            vals.append(EG / ((1 - self.epsM) * a) - 1)
        else:
            vals.append(EG / (1 - self.epsM) - 1)
            am = A(min(aF, Fr(self.acap)))
            h = 2 * am - am * am.log() - 1 - self.etaB * (am - 1) - (am.log() + self.etaB) / self.LQ
            if h > 0:
                vals.append(EG / ((1 - self.epsM) * h) - 1)
        for c in self.cs:
            r = (-A(c) * a + self.RB[c]).exp()
            if r < 1:
                vals.append(r / (1 - r))
        best = up(vals[0])
        for v in vals[1:]:
            if up(v) < best:
                best = up(v)
        return best

    def Psi(self, t):
        return self.psi((self.sF / Fr(t) - 1) / 2)

    def integrals(self, n=4000, tmin=Fr(1, 50)):
        """I = int_0^1 Psi/t^2, J = int_0^1 Psi/t^4, K3 = sup Psi/t^3, Psi1 = Psi(1); Psi increasing in t (sec 5.3)"""
        s = self.s
        # Rankin-only majorant on (0, tmin]: Psi(t) <= r(t)/(1 - r(tmin)), r(t) = exp(c/2 + RB - c s/(2t))
        bestc = None
        for c in self.cs:
            rmin = (A(c) / 2 + self.RB[c] - A(c) * s / (2 * A(tmin))).exp()
            if rmin < arb(1) / 2:
                k = A(c) * s / 2
                if bestc is None or k > bestc[1]:
                    bestc = (c, k, rmin)
        if bestc is None:
            raise SystemExit("tail: tmin too large for Rankin")
        c, k, rmin = bestc
        pref = (A(c) / 2 + self.RB[c]).exp() / (1 - rmin)
        W = 1 / A(tmin)
        I = pref * (-k * W).exp() / k
        J = pref * (-k * W).exp() * (W * W / k + 2 * W / k ** 2 + 2 / k ** 3)
        if not (W >= 3 / k):
            raise SystemExit("tail: W < 3/k")
        K3 = pref * W ** 3 * (-k * W).exp()
        # grid on [tmin, 1]; cumulative-from-the-top arrays give the integrals over [t_i, 1] (restricted Part B)
        ts = [tmin + (1 - tmin) * Fr(i, n) for i in range(n + 1)]
        dI = []; dJ = []; dK = []
        for i in range(n):
            t0, t1 = ts[i], ts[i + 1]
            ps = self.Psi(t1)
            dI.append(up(ps * (1 / A(t0) - 1 / A(t1))))
            dJ.append(up(ps * (1 / A(t0) ** 3 - 1 / A(t1) ** 3) / 3))
            dK.append(up(ps / A(t0) ** 3))
        Itop = [arb(0)] * (n + 1); Jtop = [arb(0)] * (n + 1); Ktop = [arb(0)] * (n + 1)
        for i in range(n - 1, -1, -1):
            Itop[i] = up(Itop[i + 1] + dI[i]); Jtop[i] = up(Jtop[i + 1] + dJ[i])
            Ktop[i] = dK[i] if dK[i] > Ktop[i + 1] else Ktop[i + 1]
        self.grid = (ts, Itop, Jtop, Ktop)
        I = up(I + Itop[0]); J = up(J + Jtop[0]); K3 = up(K3); K3 = K3 if K3 > Ktop[0] else Ktop[0]   # exact points: total order
        return up(I), up(J), up(K3), self.Psi(1)

    def B_rel_restricted(self, u_hi, IJK):
        """Part B over q in (Q2, z], z <= e^{u_hi}, u >= u_lo... : t = log q/u >= t2 = log Q2/u >= log Q2/u_hi.
        Uses the integrals over [t_i, 1] with t_i <= log Q2/u_hi (grid point below). Returns (I_r, J_r, K_r)."""
        ts, Itop, Jtop, Ktop = self.grid
        t2 = self.LQ / A(u_hi)
        i = 0
        while i + 1 < len(ts) and A(ts[i + 1]) < t2:
            i += 1
        if not (A(ts[i]) <= t2):
            return None
        return Itop[i], Jtop[i], Ktop[i]

    def B_rel(self, u, IJK, restricted=None):
        I, J, K3, P1 = IJK
        if restricted is not None:
            I, J, K3 = restricted
        u = A(u)
        epsV = (1 + 1 / (2 * self.LQ ** 2)) / (1 - arb(1) / self.Q2)
        return up(epsV * (I + (2 * P1 + 4 * J + K3) / (2 * u * u)) / (1 - 1 / (2 * u * u)))

    def K_B(self):
        gB = min(Fr(1), (self.sF - 1) / 2)
        return up((1 / A(gB) ** 2) * (1 / (2 * self.LQ ** 2)) * (1 + 2 / self.LQ ** 2))
