"""explore53: Arb certificate of per-gap validity (Lemma 2, REPORT.md) for the J = 4 witness at kappa = 3.34 (degree 64).
Exact rationals: a0 is raised by DELTA (monotone: d rho/d a0 = beta^2/(4p^2) >= 0), phi32 := -a3 exactly.
Checks (python-flint, 128 bits):
 (1) p(x) = sum a_j x^j > 0 on [1, inf): p(1) > 0 and no root of p with real part >= 1 on the real axis;
 (2) n3 > 0 and X_N = largest real root of N (N = 4p*gam - x q^2, cubic); rho(x) > 0 for x > max(X_N, 1);
 (3) for g in (0, g*], g* = pi/sqrt(X_N): every x_k = (k pi/g)^2 > X_N, all terms positive -> valid;
 (4) for g in [g*, pi]: cover by intervals G; S_odd(G) >= sum_{k <= K odd} rho(x_k) > 0 and likewise S_even, with
     K^2 > X_N so the dropped terms are positive.
usage: cert4.py [ncells] [K]"""
import sys, json
from fractions import Fraction as Fr
from flint import arb, acb, fmpq, fmpq_poly, ctx
ctx.prec = 128

W = json.load(open('logs/scanJ_r2_d64.json'))[0]           # kappa = 3.34
assert abs(W['kappa'] - 3.34) < 1e-12
DELTA = Fr(2, 100000)
a = [Fr(x).limit_denominator(10 ** 12) for x in W['a4']]
ph12 = Fr(W['phi4']['(1, 2)']).limit_denominator(10 ** 12)
# phi12 telescopes (global value unchanged); +2e-6 makes the small-gap coefficient n3 > 0 (the SDP left it at -8e-10)
EPS12 = Fr(2, 10 ** 6)
ph12 += EPS12
a[0] += DELTA
ph32 = -a[3]
c = a[3] + ph32; e = 2 * a[3] + ph32
n = [4 * a[0] * c,
     4 * (a[0] * a[4] + a[1] * c) - ph12 ** 2,
     4 * (a[1] * a[4] + a[2] * c) + 2 * ph12 * e,
     4 * (a[2] * a[4] + a[3] * c) - (e * e - 4 * ph12 * a[4])]
q = lambda x: fmpq(x.numerator, x.denominator)
P = fmpq_poly([q(x) for x in a]); N = fmpq_poly([q(x) for x in n])
print('rational form: a =', [str(x) for x in a], 'phi12 =', ph12, 'phi32 =', ph32)
# (1)
p1 = sum(a)
roots = P.complex_roots()
real_ge1 = [r for (r, m) in roots if r.imag.contains(0) and not (r.real < 1)]
print('(1) p(1) =', float(p1), '> 0:', p1 > 0, '; roots of p:', [(str(r.real.mid())[:10], str(r.imag.mid())[:10]) for r, m in roots], '; real roots >= 1:', len(real_ge1))
ok1 = p1 > 0 and not real_ge1
# (2)
rN = [r for (r, m) in N.complex_roots() if r.imag.contains(0)]
XN = max([r.real.upper() for r in rN] + [arb(1)])
print('(2) n3 =', float(n[3]), '> 0:', n[3] > 0, '; real roots of N:', [str(r.real.mid())[:12] for r in rN], '; X_N upper =', XN)
ok2 = n[3] > 0
# (3)-(4)
pi = arb.pi()
gstar = pi / XN.sqrt()
ncell = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
K = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
K = max(K, int(float(XN.sqrt().upper())) + 2)
aa = [arb(q(x)) for x in a]; nn = [arb(q(x)) for x in n]
def rho(x):
    p = aa[0] + x * (aa[1] + x * (aa[2] + x * (aa[3] + x * aa[4])))
    return (nn[0] + x * (nn[1] + x * (nn[2] + x * nn[3]))) / (4 * p)
lo = gstar.lower(); hi = pi.upper()
def sums(G):
    so = arb(0); se = arb(0)
    for k in range(1, K + 1):          # k > K dropped: x_k >= 4 > X_N there, so those terms are positive
        x = (k * pi / G) ** 2
        if k % 2: so += rho(x)
        else: se += rho(x)
    return so, se
sp = sums(arb(pi))
print('    point values at g = pi: S_odd =', sp[0].str(5), ' S_even =', sp[1].str(5))
worst = [None, None]
ok4 = True
stack = [(lo + (hi - lo) * i / ncell, lo + (hi - lo) * (i + 1) / ncell, 0) for i in range(ncell)]
ncells_used = 0
while stack:
    g0, g1, dep = stack.pop()
    G = arb((g0 + g1) / 2, (g1 - g0) / 2)
    so, se = sums(G)
    if so > 0 and se > 0:
        ncells_used += 1
        for j, s in enumerate((so, se)):
            if worst[j] is None or s.lower() < worst[j][0]:
                worst[j] = (s.lower(), float(G.mid()))
        continue
    if dep >= 30:
        ok4 = False
        print('  FAIL cell', G, so, se); break
    m = (g0 + g1) / 2
    stack.append((g0, m, dep + 1)); stack.append((m, g1, dep + 1))
print('(3) g* =', float(gstar.mid()), '; (4) cells used', ncells_used, 'K', K, ': min lower S_odd', worst[0], ' S_even', worst[1])
print('CERTIFIED per-gap validity:', ok1 and ok2 and ok4)
json.dump({'a': [str(x) for x in a], 'phi12': str(ph12), 'phi32': str(ph32), 'delta': str(DELTA), 'eps12': str(EPS12), 'ok': ok1 and ok2 and ok4,
           'XN_upper': str(XN), 'cells': ncell, 'K': K, 'min_Sodd': str(worst[0]), 'min_Seven': str(worst[1])},
          open('logs/cert4_k3.34.json', 'w'), indent=1)
