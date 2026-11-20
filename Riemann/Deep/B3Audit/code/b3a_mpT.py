# b3audit F6: one dual term T(n) entirely in mpmath (dps 25): A_{-1} by trapezoid in y (step and range as b3a_lemB.py, which mpmath.quad
# confirmed to 1e-13), T_pm by trapezoid in delta, Hankel functions mp.hankel1/2. T = 1e4, theta = 0.40, shifts 0; compared with the float route.
import math, mpmath as mp, numpy as np, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from b3a_lemB import run
mp.mp.dps = 25
T, th = 1e4, 0.40; rs = [0.003, 0.3]
Af, dlf, rows = run(T, th, 0.0, 0.0, rs)          # float route (prints its rows)
H = mp.mpf(T)**mp.mpf(th); Tm = mp.mpf(T); K = Tm/H; tau = Tm/H**2
nd = len(dlf); dl = [mp.mpf(float(x)) for x in dlf]; hdl = dl[1] - dl[0]
Y = mp.mpf('6.5'); dmax = max(abs(dl[0]), abs(dl[-1]))
wmax = 2*K*max(abs(1 - (1 + s)*mp.exp(2*t/H)) for s in (-dmax, dmax) for t in (-Y, Y))
h = 2*mp.pi/(wmax + 30); ny = int(2*Y/h) + 2; ys = [-Y + 2*Y*i/(ny - 1) for i in range(ny)]; hy = ys[1] - ys[0]
em1 = [mp.expm1(2*y/H) for y in ys]; c0 = [2*y/H - e for y, e in zip(ys, em1)]; amp = [mp.exp(-y*y + y/H)*hy/mp.sqrt(mp.pi) for y in ys]
A = []
for d in dl:
    A.append(mp.fsum(a*mp.expj(Tm*(c - d*e)) for a, c, e in zip(amp, c0, em1)))
print("max |A_mp - A_float| / max|A| =", float(max(abs(a - complex(b)) for a, b in zip(A, Af))/max(abs(a) for a in A)), flush=True)
xi = [Tm*(1 + d)/(2*mp.pi) for d in dl]; dxi = Tm*hdl/(2*mp.pi)
for r, row in zip(rs, rows):
    w = r*tau; X0 = 2*mp.sqrt(2*mp.pi)*mp.sqrt(w*Tm)
    Tp = -mp.pi*mp.fsum(a*(-0.5j)*mp.hankel1(0, 4*mp.pi*mp.sqrt(w*x)) for a, x in zip(A, xi))*dxi
    Tq = -mp.pi*mp.fsum(a*(0.5j)*mp.hankel2(0, 4*mp.pi*mp.sqrt(w*x)) for a, x in zip(A, xi))*dxi
    sc = H/mp.sqrt(X0)
    print(f"mp: r={r} |T+|/HX0^-1/2={float(abs(Tp)/sc):.8f} |T-|={float(abs(Tq)/sc):.8f} |T|={float(abs(Tp+Tq)/sc):.8f}  float: {row[4]:.8f} {row[5]:.8f} {row[6]:.8f}", flush=True)
