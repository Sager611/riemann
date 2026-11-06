#!/usr/bin/env python3
"""Numeric constants for the smoothed prime-ideal count (Riemann/Deep/KoiranAM/DERIVATION.md, Lemma A-num).
Reuses the Arb kernel class of Riemann/Deep/ChebExp/cert.py with L = E = K, H = 1 (no Chebotarev condition).
  uv run --no-project --with python-flint python3 count_cert.py [h,s*,s0 ...]
Prints rigorous enclosures (upper endpoints for errors, lower for main terms):
  P1 = Phi(1) (lower), P0 = Phi(0) (upper), kn = max F (upper), A1 = sup_R on Re s = 31/32, A2 = on Re s = 1/32,
  g = g(s0) (lower), cs = 1/s0 + 1/(s0-1) (upper), S2, S3, Ph2, Ph3 (upper)."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ChebExp'))
from flint import arb, acb, ctx
import cert
ctx.prec = 120
delta = arb(1) / 32
def F(x):
    return float(x.mid())
def U(x):
    return float(cert.up(x).mid())
def Lo(x):
    return float(cert.lo(x).mid())
kernels = [(1.0, 0.0, 1.3), (1.0, 0.0, 1.5), (1.5, 0.0, 1.4), (1.0, 0.25, 1.4), (0.75, 0.0, 1.4)]
if len(sys.argv) > 1:
    kernels = [tuple(map(float, a.split(','))) for a in sys.argv[1:]]
t0 = time.time()
for (h, ss, s0) in kernels:
    K = cert.Kernel(h, ss, s0).consts()
    A1, A1lo, c1 = K.sup_line(1 - delta, 1e-3)
    A2, A2lo, c2 = K.sup_line(delta, 0.5)
    g = cert.g_gamma(K.s0)
    cs = cert.up(1 / K.s0 + 1 / (K.s0 - 1))
    print("kernel h=%g s*=%g s0=%g (U = 3h = %g)" % (h, ss, s0, 3 * h))
    print("  P1>=%.7g  P0<=%.7g  kn<=%.7g  A1<=%.7g (point value %.7g, %d boxes)  A2<=%.7g (%d boxes)" % (
        Lo(K.P1), U(K.P0), U(K.kn), U(A1), Lo(A1lo), c1, U(A2), c2))
    print("  g(s0)>=%.6g  1/s0+1/(s0-1)<=%.7g  S2<=%.6g S3<=%.6g Ph2<=%.6g Ph3<=%.6g" % (
        Lo(g), U(cs), U(K.S2), U(K.S3), U(K.Ph2), U(K.Ph3)))
    print("  P1/kn>=%.6g   A1/kn<=%.6g   A1/P1<=%.6g" % (Lo(K.P1) / U(K.kn), U(A1) / Lo(K.kn) if False else U(A1) / U(K.kn) * (U(K.kn) / Lo(K.kn) if Lo(K.kn) > 0 else 1), U(A1) / Lo(K.P1)))
    sys.stdout.flush()
print("(%.1f s)" % (time.time() - t0))


# ---------------------------------------------------------------- Lemma A-num for the chosen kernel (h, s*, s0) = (1.5, 0, 1.4)
def lemma_num():
    K = cert.Kernel(1.5, 0.0, 1.4).consts()
    A1, _, _ = K.sup_line(1 - delta, 1e-3)
    A2, _, _ = K.sup_line(delta, 0.5)
    cs = cert.up(1 / K.s0 + 1 / (K.s0 - 1))
    lx = arb(40)                                   # x >= e^40
    # M_x = max(A1, x^{2 delta - 1} A2) = A1 once x^{15/16} >= A2/A1
    assert (A2 / A1) < (lx * 15 / 16).exp()
    # P_k(x) <= 1.01624 [x^{1/2}(S2+Ph2/2) + x^{1/3}(S3+Ph3/3) + (J-3) x^{1/4}(kn+P0/4)], J <= log x/log 2 + 1;
    # (J-3) x^{-1/4} and x^{-1/6} decrease for x >= e^40, so their values at e^40 bound them
    J = lx / cert.LOG2 + 1
    pk = cert.THETA_C * ((K.S2 + K.Ph2 / 2) + (-lx / 6).exp() * (K.S3 + K.Ph3 / 3) + (J - 3) * (-lx / 4).exp() * (K.kn + K.P0 / 4))
    triv = 2 * K.P0 * (-lx / 2).exp()               # T(x) <= e P0 (1 + eps_x) <= 2 e P0 (eps_x <= 1 for x >= 2e^U); relative to e x^{1/2}
    cP = cert.up(pk + triv)                          # e*(P0(1+eps_x) + P_k(x)) <= cP * e * x^{1/2}
    a0 = cert.lo(K.P1 / K.kn)                        # main coefficient (lower)
    a1 = cert.up(A1 / 2 / K.kn)                      # coefficient of (Lambda + 2 cs) x^{-1/32}
    a2 = cert.up(2 * cs)
    a3 = cert.up(cP / K.kn)
    print("Lemma A-num (x >= e^40): #{deg-1 primes, x e^-4.5 <= Np <= x} >= (x/log x)[a0 - a1 (log Delta_K + a2) x^(-1/32) - a3 e x^(-1/2)]")
    print("  a0 >= %.7g   a1 <= %.7g   a2 <= %.7g   a3 <= %.7g   (cP <= %.7g, A1 <= %.7g, kn <= %.7g, P1 >= %.7g)" % (
        Lo(a0), U(a1), U(a2), U(a3), U(cP), U(A1), U(K.kn), Lo(K.P1)))
    # rounded constants used in the text, checked in the safe direction
    R0, R1, R2, R3 = arb('0.2772'), arb('0.4953'), arb('6.43'), arb('0.7233')
    print("  rounded 0.2772 <= a0:", bool(R0 <= a0), " 0.4953 >= a1:", bool(R1 >= a1), " 6.43 >= a2:", bool(R2 >= a2), " 0.7233 >= a3:", bool(R3 >= a3))
    # x0 conditions: bracket >= R0/2 if R1 (logD + R2) x^{-1/32} <= R0/4 and R3 e x^{-1/2} <= R0/4
    c32 = 4 * R1 / R0; c2 = 4 * R3 / R0
    print("  x^(1/32) >= %.6g (log D + 6.43)  and  x^(1/2) >= %.6g e  give bracket >= %.6g" % (U(c32), U(c2), Lo(R0 / 2)))
    # then pi_S(x) >= (R0/2) x/(e log x) - omega(a a_B) >= B0 if x/log x >= (2/R0) e (B0 + omega); with x/log x >= sqrt(x) (x >= e^40)
    print("  sqrt(x) >= %.6g e B'  (B' = B0 + log2|a a_B|), i.e. x >= %.6g (e B')^2" % (U(2 / R0), U((2 / R0) ** 2)))
lemma_num()
