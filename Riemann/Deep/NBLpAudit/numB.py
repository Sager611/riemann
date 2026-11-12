# nblpaudit numerics: A4 residue identity J, B(eps) and the A4 bound chain; Lemma A shape (|zeta|^{+-1} on lines).
# usage: python3 numB.py T [fine]  (T = truncation height; 'fine' = panels 0.01 on [0,1], 0.05 on [1,5]: the pole
# s=1 sits 1-sigma' from the line, which spoils 8-point Gauss on width-0.5 panels near t=0 when 1-sigma' ~ 0.14)
import sys, math, time
import numpy as np
from flint import acb, arb, ctx
ctx.prec = 64
LOGPI = math.log(math.pi)
T = float(sys.argv[1])
def zc(x, t):
    z = acb(x, t).zeta(); return complex(z.real.mid(), z.imag.mid())
def zr(x): z = acb(x).zeta(); return float(z.real.mid())
def gr2(x, t, e):
    s = acb(x, t); s2 = acb(x + e, t)
    lg = (s2 / 2).lgamma() - (s / 2).lgamma()
    a = abs(complex((s2 * (s2 - 1)).real.mid(), (s2 * (s2 - 1)).imag.mid())) / abs(complex((s * (s - 1)).real.mid(), (s * (s - 1)).imag.mid()))
    return (math.exp(float(lg.real.mid()) - e / 2 * LOGPI) * a) ** 2
xg, wg = np.polynomial.legendre.leggauss(8)
CG = 5993.0534
t0 = time.time()
for p, e in [(1.3, 0.2), (1.1, 0.1), (8 / 7, 0.05)]:
    sp = 1 / p - e / 2
    Jf = 2 * math.pi * (zr(2 * sp) / zr(2 * sp + e) / (2 * sp) - 1 / ((2 * sp - 1) * zr(1 + e)))
    IJ = IB = IC = 0.0; marks = {}; last = []
    FINE = len(sys.argv) > 2
    pans = ([(a, 0.01) for a in np.arange(0, 1, 0.01)] + [(a, 0.05) for a in np.arange(1, 5, 0.05)] + [(a, 0.5) for a in np.arange(5, T, 0.5)]) if FINE else [(a, 0.5) for a in np.arange(0.0, T, 0.5)]
    for a, hw in pans:
        for xn, wn in zip(xg, wg):
            t = a + hw / 2 * (xn + 1); w = hw / 2 * wn / (sp * sp + t * t)
            R = zc(sp, t) / zc(sp + e, t)
            IJ += 2 * w * (R.real - 1); IB += 2 * w * abs(1 - R) ** 2; IC += 2 * w * (gr2(sp, t, e) - 1)
            if t > T / 2: last.append(abs(1 - R) ** 2)
        if abs(a + hw - T / 4) < 1e-9 or abs(a + hw - T / 2) < 1e-9: marks[round(a + hw)] = math.pi / sp + IJ
    J = math.pi / sp + IJ
    tailB = 2 * np.mean(last) / T
    twopiB = IB + tailB
    rhs1 = math.pi / sp - 2 * J + (math.pi / sp + IC)
    rhs2 = (2 * math.pi / sp) * (1 - zr(2 * sp) / zr(2 * sp + e)) + 4 * math.pi / ((2 * sp - 1) * zr(1 + e)) + CG * e
    Kp = 4 * p / (2 - p) + CG / (2 * math.pi)
    print("p=%.4f eps=%.2f sigma'=%.4f  J(T/4)=%.8f J(T/2)=%.8f J(T)=%.8f  J_formula=%.8f  diff=%.1e" % (p, e, sp, marks.get(round(T / 4), float('nan')), marks.get(round(T / 2), float('nan')), J, Jf, J - Jf))
    print("   2piB num=%.6f (tail est %.1e)  <= pi/s'-2J+I_c = %.6f  <= final RHS %.4f ;  B=%.3e  K_p*eps=%.3f  B/eps^2=%.3f" % (twopiB, tailB, rhs1, rhs2, twopiB / (2 * math.pi), Kp * e, twopiB / (2 * math.pi) / e ** 2))
    sys.stdout.flush()
print("A4 time %.1f" % (time.time() - t0)); sys.stdout.flush()
# Lemma A shape: max of log|zeta|^{+-1} / log(2+t) over t-ranges
for x in (0.6, 0.75, 0.9):
    row = []
    for lo, hi in ((6, 60), (60, 600), (600, 6000)):
        tt = np.linspace(lo, hi, 6000 if hi > 600 else 3000)
        L = np.array([math.log(abs(zc(x, t))) for t in tt])
        row.append("[%g,%g]: max log|z|=%.3f max -log|z|=%.3f max(log|z|^±)/log(2+t)=%.3f" % (lo, hi, L.max(), -L.min(), np.max(np.abs(L) / np.log(2 + tt))))
    print("sigma=%.2f " % x + " | ".join(row)); sys.stdout.flush()
print("total time %.1f" % (time.time() - t0))
