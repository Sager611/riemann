
# check_good2.py: kappa^l J_l(N) by quadrature (python-flint, exact polynomials, high precision) vs the
# rigorous Good lower bound exp(nu Phi + gf) used in assembly.py  (must have log J >= log LB)
import sys, math
import flint
from flint import arb, acb, fmpq, acb_poly
from math import comb, factorial
import assembly as A
m = int(sys.argv[1]); cases = [tuple(int(v) for v in c.split(':')) for c in sys.argv[2].split(',')]
flint.ctx.prec = int(sys.argv[3]) if len(sys.argv) > 3 else 1600
nu0 = arb('203.5'); A.load(nu0, 'tab203.json', 'lg203.pkl')
nu = arb(m) + arb(1) / 2
P = acb_poly([acb(arb(fmpq(comb(m, j) * factorial(2*m - j), factorial(2*m)))) for j in range(m + 1)])
EJ = 1300
Epoly = acb_poly([acb(arb(fmpq(factorial(m) * factorial(m + j), factorial(j) * factorial(2*m) * factorial(2*m + j + 1)))) for j in range(EJ)])
th0 = arb(factorial(2*m)) / (factorial(m) * 2**m)
logLam = (th0 * nu.sqrt() * nu**(-nu)).log()
yJ = arb('0.85')
def logF(l, n, zeta):
    # log of (e^z X^l / p) zeta^{-n-1} / Lam   (complex log, principal parts only matter for exp)
    z = 2 * nu * zeta
    p = P(z); E = Epoly(z)
    X = z**(2*m + 1) * E / p
    return z - p.log() + l * X.log() - (n + 1) * zeta.log() - logLam
def simpson_sum(vals, h):
    # vals: list of complex logs of integrand values at equally spaced nodes (odd count)
    mx = max(float(v.real.mid()) for v in vals)
    s = acb(0)
    N = len(vals) - 1
    for k, v in enumerate(vals):
        w = 1 if k in (0, N) else (4 if k % 2 == 1 else 2)
        s += w * (v - mx).exp()
    return s * h / 3, mx
for (n, l) in cases:
    a = 2 * l + 1; t = arb(n) / nu - 2 * l
    r = A.saddle(a, t); rf = float(r.mid())
    gf = A.good_factor(l, rf * (1 - 1e-9), rf * (1 + 1e-9), nu0)
    if gf is None: print('n=%d l=%d not good (r*=%.4f)' % (n, l, rf), flush=True); continue
    logLB = float((nu * A.Phi(a, t) + gf).mid())
    rr = arb(rf)
    NA_ = 2400
    if rf <= 0.85 * 0.999:
        phic = arb.pi() / 2
    else:
        phic = (yJ / rr).asin()
    h = 2 * phic / NA_
    vals = []
    for k in range(NA_ + 1):
        ph = -phic + h * k
        zeta = acb(rr * ph.cos(), rr * ph.sin())
        vals.append(logF(l, n, zeta) + zeta.log())     # integrand of (1/2pi) int zeta g dphi
    sA, mA = simpson_sum(vals, h)
    arc = sA / (2 * arb.pi())
    NS = 600
    if rf <= 0.85 * 0.999:
        # upper segment y in [r, yJ]: (1/2pi i) int g d(i y) = (1/2pi) int g dy ; lower = conj
        hs = (yJ - rr) / NS
        vs = [logF(l, n, acb(0, rr + hs * k)) for k in range(NS + 1)]
        sS, mS = simpson_sum(vs, hs)
        seg = 2 * (sS / (2 * arb.pi())).real
    else:
        xc = (rr * rr - yJ * yJ).sqrt(); hs = xc / NS
        vs = [logF(l, n, acb(xc - hs * k, yJ)) for k in range(NS + 1)]   # from xc+i yJ to i yJ: dzeta = -dx
        sS, mS = simpson_sum(vs, hs)
        seg = 2 * (-sS / (2 * arb.pi() * acb(0, 1))).real * 1
        # (1/2pi i) int_{xc}^{0} g(x+iyJ) dx = -(1/2pi i) int_0^{xc}; we summed from xc down: same orientation
        seg = 2 * (sS / (2 * arb.pi() * acb(0, 1))).real
    # combine with scaling
    totarc = arc.real * math.exp(0)  # scaled by e^{mA}
    J = float(arc.real.mid()) * math.exp(0)
    # bring both to common scale: log value
    import mpmath as mp
    vA = mp.mpf(float(arc.real.mid())) * mp.exp(mA)
    vS = mp.mpf(float(seg.mid())) * mp.exp(mS)
    Jtot = vA + vS
    if Jtot > 0:
        lj = float(mp.log(Jtot))
        print('n=%d l=%d r*=%.4f: log(kappa^l J) = %.4f  log(LB) = %.4f  margin %.4f  seg/arc=%.1e %s' % (n, l, rf, lj, logLB, lj - logLB, float(abs(vS/vA)), 'OK' if lj >= logLB else 'VIOLATION'), flush=True)
    else:
        print('n=%d l=%d: J <= 0 ?! %s' % (n, l, mp.nstr(Jtot, 5)), flush=True)
