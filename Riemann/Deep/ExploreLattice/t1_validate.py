# t1_validate.py (explore48): check Chowla-Selberg E* against L-function products for D=-20 (x^2+5y^2) and D=-15 (x^2+xy+4y^2),
# and the known off-line zeros: McPhedran's Potter-Titchmarsh zero 0.9329+15.6682i (x^2+5y^2), explore21's 0.80001+12.03860i (x^2+xy+4y^2).
import math, cmath, time, sys
sys.path.insert(0, '.')
from flint import acb, arb, dirichlet_char, ctx
from epstein_cs import Estar, F_line
from t0_timing import real_char
def Lam_from_L(s, D, combo):
    s = acb(s)
    old = ctx.prec; ctx.prec = 120
    Z = acb(0)
    for coef, ds in combo:
        term = acb(coef)
        for d in ds:
            term *= s.zeta() if d == 1 else s.dirichlet_l(real_char(d))
        Z += term
    L = (arb(abs(D)).sqrt() / (2 * arb.pi())) ** s * s.gamma() * Z
    G = (acb(0, -1) * arb.pi() * s / 2).exp() * L / 2
    ctx.prec = old
    return complex(float(G.real.mid()), float(G.imag.mid()))
cases = [(-20, 0.0, math.sqrt(5.0), [(1, (1, -20)), (1, (-4, 5))]),
         (-15, -0.5, math.sqrt(15.0) / 2, [(1, (1, -15)), (1, (-3, 5))])]
for D, x, y, combo in cases:
    for s in (complex(0.7, 10.0), complex(0.5, 20.0), complex(1.3, 33.0)):
        t0 = time.time(); a = Estar(s, x, y); t1 = time.time(); b = Lam_from_L(s, D, combo)
        print('D=%d s=%s  CS=%.12g%+.12gi  L=%.12g%+.12gi  rel=%.2e  (CS %.3fs)' % (D, s, a.real, a.imag, b.real, b.imag, abs(a-b)/abs(b), t1-t0))
def newton(f, s0, it=40):
    s = s0; h = 1e-6
    for _ in range(it):
        fs = f(s); d = (f(s + h) - fs) / h
        ds = fs / d; s -= ds
        if abs(ds) < 1e-12: break
    return s, abs(f(s))
for (x, y, s0, label) in ((0.0, math.sqrt(5.0), complex(0.9329, 15.6682), 'x^2+5y^2 PT zero'), (-0.5, math.sqrt(15.0)/2, complex(0.80001, 12.0386), 'x^2+xy+4y^2 explore21')):
    s, r = newton(lambda s: Estar(s, x, y), s0)
    print(label, 'root', s, 'resid', r, ' t/(2 pi y) = %.4f' % (s.imag / (2 * math.pi * y)))
# line function check: F(t) sign vs Re of normalized E* on the line
for t in (5.0, 14.0, 30.0):
    Fv, c, b, th = F_line(t, 0.0, math.sqrt(5.0), return_parts=True)
    print('F_line t=%g F=%.6f cos=%.6f beta=%.6f' % (t, Fv, c, b))
