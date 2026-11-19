# sgfaudit falsification / consistency tests for small-gap-floor (float64 + Arb). One process, single-threaded.
import re, sys, math, resource
import numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, 'Riemann/Deep/SGFAudit')
import indep_verify as iv
from flint import arb, ctx
ctx.prec = 256
def cert(case):
    name, h, theta, lam, k0, eta, K, M2, M3 = iv.setup(case)
    return float(h), float(theta), float(lam), k0, np.array([float(e) for e in eta]), h, theta, eta
def P_np(al, h, eta):
    k = np.arange(1, len(eta) + 1)
    return (1 - h) + 2 * h * (np.cos(2 * np.pi * np.outer(al, k * h)) @ eta)
GX, GW = np.polynomial.legendre.leggauss(120)   # small rule (leggauss(4000) needs a 128 MB eigenproblem)
def gl(f, a, b, pieces=60):
    e = np.linspace(a, b, pieces + 1); tot = 0.0
    for lo, hi in zip(e[:-1], e[1:]):
        tot += (hi - lo) / 2 * np.sum(GW * f((hi - lo) / 2 * GX + (lo + hi) / 2))
    return tot
sinc = lambda x: np.sinc(x)          # numpy sinc = sin(pi x)/(pi x)
def lhs_sigma(fhat, h, eta, N=2_000_000, chunk=200_000):
    tot = 0.0
    for s0 in range(0, N, chunk):
        k = np.arange(s0 + 1, s0 + chunk + 1); w = np.full(chunk, h)
        if s0 == 0: w[:len(eta)] = h * (1 + eta)
        tot += 2 * np.sum(w * fhat(k * h))
    return tot
# F1: Poisson identity with a C^2 kernel whose support contains the dual-lattice points m/h (tests the m != 0 term)
B3 = lambda x: np.where(np.abs(x) <= 1, (4 - 6 * x**2 + 3 * np.abs(x)**3) / 6, np.where(np.abs(x) <= 2, (2 - np.abs(x))**3 / 6, 0.0))
for case in ['B', 'B1']:
    h, th, lam, k0, eta, *_ = cert(case)
    s = 5.0
    f = lambda a: B3(a / s); fhat = lambda u: s * sinc(s * u)**4
    L = lhs_sigma(fhat, h, eta)
    mlist = [m for m in range(1, 20) if m / h < 2 * s]
    R = f(0.0) + 2 * sum(float(f(m / h)) for m in mlist) + 2 * sum(gl(lambda a: f(a) * (P_np(a, h, eta) - 1), lo, hi) for lo, hi in [(0, 5), (5, 10)])
    print(f'F1 {case}: sum f^ dsigma = {L:.12f}; f(0)+sum_(m!=0) f(m/h)+int f(P-1) = {R:.12f}; diff {L-R:.1e}; m/h inside supp: {[round(m/h,3) for m in mlist]}')
# F2: Lamzouri's own kernel (l.452-469) f(u) = lam_f Rhat(lam_f u), f^(xi) = R(xi/lam_f), at lam_f = 0.7892 (> 0.775), against CertB
lf = 0.7892
Rhat = lambda u: np.where(np.abs(u) <= 1, 1 - np.abs(u) + np.sin(2 * np.pi * np.abs(u)) / (2 * np.pi), 0.0)
Rsel = lambda x: sinc(x)**2 / (1 - x**2)
f = lambda a: lf * Rhat(lf * a); fhat = lambda u: Rsel(u / lf)
for thL in [4 / 7 - 1e-5, 4 / 7]:
    Jq = f(0.0) + 2 * gl(lambda a: f(a) * a, 0, thL)
    Jc = lf + lf * thL**2 - 2 / 3 * lf**2 * thL**3 - thL * math.cos(2 * math.pi * lf * thL) / (2 * math.pi**2) + math.sin(2 * math.pi * lf * thL) / (4 * math.pi**3 * lf)
    print(f'F2 J_theta(f) at theta={thL:.8f}: quadrature {Jq:.12f}, closed form l.468 {Jc:.12f}')
h, th, lam, k0, eta, *_ = cert('B')
L = lhs_sigma(fhat, h, eta)
J = f(0.0) + 2 * gl(lambda a: f(a) * a, 0, th)
slack = 2 * gl(lambda a: f(a) * (P_np(a, h, eta) - a), 0, th) + 2 * gl(lambda a: f(a) * P_np(a, h, eta), th, 1 / lf)
print(f'F2 CertB identity: sum f^ dsigma = {L:.12f}; (J-1) + int f(P-|a|1_th) = {J-1:.12f} + {slack:.12f} = {J-1+slack:.12f}; diff {L-(J-1+slack):.1e}')
w7 = h * (1 + eta[6])
print(f'F2 atom at u = 7h = {7*h:.6f} (weight {w7:.6f}, f^ = {float(fhat(7*h)):.6f}): contributes {2*w7*float(fhat(7*h)):.6f}; rest of sigma {L-2*w7*float(fhat(7*h)):.6f}')
# F3: Theorem A identity S_a(alpha) = 1 - sum_m tri(alpha - m/a) vs (1-a) - 2a sum_k sinc^2(ka) cos(2 pi alpha k a)
for a in [7 / 11, 0.5]:
    N = 200_000; k = np.arange(1, N + 1); c = sinc(k * a)**2
    err = 0.0
    for al in [0.0, 0.1, 0.3, 1 / a - 1, (1 / a - 1) + 0.05, 1 / (2 * a), 0.9 / a, 1 / a, 1.3 / a]:
        ser = (1 - a) - 2 * a * np.sum(c * np.cos(2 * np.pi * al * k * a))
        tri = 1 - sum(max(0.0, 1 - abs(al - m / a)) for m in range(-3, 4))
        err = max(err, abs(ser - tri))
    print(f'F3 a={a:.6f}: max |series - triangle form| over 9 points = {err:.1e} (tail bound 2/(pi^2 a N) = {2/(math.pi**2*a*N):.1e}); S_a = |alpha| on [0, {1/a-1:.6f}]')
# F6: extreme entries in Arb (256 bits): P(0), P(theta) - theta (the binding left branch at the jump), P(1/(2h))
for case in ['B', 'B1']:
    name, hq, thq, lamq, k0, etaq, K, M2, M3 = iv.setup(case)
    P0 = iv.derivs(hq, etaq, Fr(0))[0]; Pt = iv.derivs(hq, etaq, thq)[0] - iv.A(thq); Pm = iv.derivs(hq, etaq, Fr(1, 2) / hq)[0]
    print(f'F6 {case}: P(0) = {P0.str(12)}; P(theta) - theta = {Pt.str(12)}; P(1/(2h)) = {Pm.str(12)}')
print(f'maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB')
