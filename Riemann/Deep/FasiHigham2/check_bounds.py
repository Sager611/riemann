
# check_bounds.py: numerically integrate |g_l| (exact functions, Arb at high precision so that the
# E-series has no cancellation problem) along the crude contour and compare with the rigorous crude
# bound used in assembly.py (the bound must be >= the numeric value). Cases n:l with N >= 0.
import sys, math
import flint
from flint import arb, acb, fmpq, acb_poly
from math import comb, factorial
import assembly as A
m = int(sys.argv[1]); cases = [tuple(int(v) for v in c.split(':')) for c in sys.argv[2].split(',')]
flint.ctx.prec = int(sys.argv[3]) if len(sys.argv) > 3 else 1800
nu0 = arb('203.5'); A.load(nu0, 'tab203.json', 'lg203.pkl')
nu = arb(m) + arb(1) / 2
P = acb_poly([acb(arb(fmpq(comb(m, j) * factorial(2*m - j), factorial(2*m)))) for j in range(m + 1)])
EJ = 2000
Epoly = acb_poly([acb(arb(fmpq(factorial(m) * factorial(m + j), factorial(j) * factorial(2*m) * factorial(2*m + j + 1)))) for j in range(EJ)])
th0 = arb(factorial(2*m)) / (factorial(m) * 2**m)
logLam = (th0 * nu.sqrt() * nu**(-nu)).log()
logkap = ((2 * arb.pi() * nu).sqrt() * nu**nu * (-nu).exp() / (nu + 1).gamma()).log()
YJ = 0.85
def loggabs(l, n, zeta):
    # log |(e^z X^l / p) zeta^{-n-1}| / (Lam kappa^l)  -- the modulus of the J_l integrand (times zeta)
    z = 2 * nu * zeta
    p = P(z); E = Epoly(z)
    X = z**(2*m + 1) * E / p
    v = (z.real - abs(p).log() + l * abs(X).log() - (n + 1) * abs(zeta).log()) - logLam - l * logkap
    return float(v.mid())
def contour(rc):
    pts = []
    if rc <= YJ * 0.999:
        phic = math.pi / 2; K = 400
        for k in range(K):
            ph = phic * (k + 0.5) / K; pts.append((acb(rc * math.cos(ph), rc * math.sin(ph)), rc * phic / K))
        K2 = 200
        for k in range(K2):
            y = rc + (YJ - rc) * (k + 0.5) / K2; pts.append((acb(0, y), (YJ - rc) / K2))
    else:
        phic = math.asin(YJ / rc); K = 400
        for k in range(K):
            ph = phic * (k + 0.5) / K; pts.append((acb(rc * math.cos(ph), rc * math.sin(ph)), rc * phic / K))
        xc = math.sqrt(rc**2 - YJ**2); K2 = 300
        for k in range(K2):
            x = xc * (k + 0.5) / K2; pts.append((acb(x, YJ), xc / K2))
    return pts
for (n, l) in cases:
    a = 2*l + 1; t = n / float(nu.mid()) - 2*l
    if t < 0:
        print('n=%d l=%d: N < 0, not crude-bounded in the proof (Neg sum, see check_neg_IL.py)' % (n, l)); continue
    r = float(A.saddle(a, arb(t)).mid()) if t > 0 else 0.0
    rc = A.RMAX if r > A.RMAX else max(r, 0.01)
    pcs = A.crude_pieces(l, t, t, rc, nu0)
    lb = max(float(nu.mid()) * float(R.mid()) + float(c.mid()) + math.log(float(ln.mid())) for (R, c, ln) in pcs)
    logs = [loggabs(l, n, z) + math.log(w) for (z, w) in contour(rc)]
    mx = max(logs)
    tot = mx + math.log(sum(math.exp(v - mx) for v in logs)) - math.log(2 * math.pi) + math.log(2)   # both halves
    print('n=%d l=%d t=%.4f r*=%.4f rc=%.4f : log numeric (1/2pi)int|g| = %.3f ; log(crude bound, max piece) = %.3f  %s' % (n, l, t, r, rc, tot, lb, 'OK' if lb >= tot - 1e-6 else 'VIOLATION'), flush=True)
