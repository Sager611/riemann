# hcaudit independent checks of cert768 (EPS = 1/1000 as in the certificate): root isolation, Descartes, M3 sanity,
# fine-lattice and local-refinement sup of Phi, R^ = FT of R. Usage: EPS=1/1000 python3 hc_indep.py cert768.json
import sys, json, os
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpz_poly, fmpq, arb, arb_poly, ctx
import numpy as np
ctx.prec = 300
d = json.load(open(sys.argv[1]))
lam = Fr(d['lam']); th = Fr(d['theta']); D = d['D']; s = Fr(d['s']); U = Fr(d['U']); dg = Fr(d['dg'])
B = [Fr(b) for b in d['B']]; B[0] += Fr(os.environ.get('EPS', '0')); A = [Fr(a) for a in d['A']]; beta = Fr(d['beta'])
q2 = lambda x: fmpq(x.numerator, x.denominator)
H = [fmpq_poly([1]), fmpq_poly([0, 2])]
for n in range(1, 2*D): H.append(fmpq_poly([0, 2])*H[n] - 2*n*H[n-1])
P = sum((q2(B[n])*H[2*n] for n in range(D+1)), fmpq_poly([0])); Q = sum((q2(B[n]*(-1)**n)*H[2*n] for n in range(D+1)), fmpq_poly([0]))
p = fmpq_poly([P[2*k] for k in range(D+1)]); q = fmpq_poly([Q[2*k] for k in range(D+1)])
assert all(P[2*k+1] == 0 for k in range(D)) and all(Q[2*k+1] == 0 for k in range(D))
# (a) independent real-root isolation (Arb complex_roots: real roots have exactly zero imaginary part)
def realroots(poly):
    den = poly.denom(); z = fmpz_poly([int((poly[i]*den).p) for i in range(poly.degree()+1)])
    rr = [r for r, m in z.complex_roots() if r.imag == 0]
    return [r.real for r in rr]
rp = realroots(p); rq = realroots(q)
w0 = 2*Fr(314159265, 100000000)*s*s*U*U
print('(a) deg p = %d, real roots of p: %s' % (p.degree(), [x.str(8) for x in rp]))
print('    deg q = %d, real roots of q: %s ; w0 = %.9f, 2 pi s^2 U^2 = %s' % (q.degree(), [x.str(8) for x in rq], float(w0), (2*arb.pi()*arb(q2(s*s*U*U))).str(12)))
print('    p > 0 on [0,inf): %s ; q < 0 on [2 pi s^2U^2, inf): %s' % (all(x < 0 for x in rp) and p(fmpq(0)) > 0,
      all(x < arb(q2(w0)) for x in rq) and q[q.degree()] < 0))
sv = lambda c: sum(1 for i in range(len(c)-1) if c[i]*c[i+1] < 0)
cp = [p[i] for i in range(p.degree()+1) if p[i] != 0]
print('    Descartes sign variations of p (bound on #positive roots): %d' % sv(cp))
# f^ and derivatives in Arb (variable v = sqrt(2 pi) s u)
pi = arb.pi(); sA = arb(q2(s)); k_ = (2*pi).sqrt()*sA
QA = arb_poly([arb(Q[i]) for i in range(Q.degree()+1)]); Q1 = QA.derivative(); Q2_ = Q1.derivative(); Q3 = Q2_.derivative()
def fh(u): v = k_*u; return sA*QA(v)*(-(v*v)/2).exp()
def fh3(u):
    v = k_*u; return sA*k_**3*(Q3(v) - 3*v*Q2_(v) + (3*v*v - 3)*Q1(v) + (3*v - v**3)*QA(v))*(-(v*v)/2).exp()
m3 = max(abs(float(fh3(arb(i)/1000).mid())) for i in range(0, 8001))
print('(b) numerical max|fhat\'\'\'| on [0,8] (step 1e-3) = %.3f  (checker bound M3 = 3536.5)' % m3)

# (c) sup of Phi on [lam,U]^2: fine lattice dg/40 (Arb values, float max), then local refinement around the top points
nodes = [lam + k*dg for k in range(len(A))]
def Aval(x):   # exact piecewise linear, constant tail
    if x >= U: return A[-1]
    k = int((x - lam)/dg); k = min(k, len(A)-2); t = (x - nodes[k])/dg
    return A[k] + (A[k+1] - A[k])*t
mf = 40; de = dg/mf; Nf = 77*mf
up = lambda z: float(z.upper())
F1 = np.array([up(fh(arb(q2(lam + i*de)))) for i in range(Nf+1)])
F2 = np.array([up(fh(arb(q2(2*lam + k*de)))) for k in range(2*Nf+1)])
Al = np.array([float(Aval(lam + i*de)) for i in range(Nf+1)]); xs = np.array([float(lam + i*de) for i in range(Nf+1)])
bf = float(beta); best = []; 
for i in range(Nf+1):
    row = F1[i] + F2[i:i+Nf+1] + Al[i] - Al + bf*(1 - xs[i]); j = int(row.argmax()); best.append((row[j], i, j))
best.sort(reverse=True)
print('(c) fine lattice (delta = dg/40) max Phi = %.10f at (x,y) = (%.5f, %.5f); checker corner max (dg/5) = -0.3478888283' % (best[0][0], xs[best[0][1]], xs[best[0][2]]))
seen = []
for v, i, j in best:
    if all(abs(i - a) + abs(j - b) > 2*mf for a, b in seen): seen.append((i, j))
    if len(seen) == 6: break
glob = -9
for i, j in seen:
    x0, y0 = lam + i*de, lam + j*de; loc = -9; arg = None
    for a in range(-20, 21):
        for b in range(-20, 21):
            x = x0 + a*de/10; y = y0 + b*de/10
            if x < lam or y < lam or x > U or y > U: continue
            val = fh(arb(q2(x))) + fh(arb(q2(x + y))) + arb(q2(Aval(x) - Aval(y) + beta*(1 - x)))
            if up(val) > loc: loc = up(val); arg = (float(x), float(y))
    glob = max(glob, loc)
    print('    local max near (%.4f, %.4f): %.10f at (%.6f, %.6f)' % (float(x0), float(y0), loc, arg[0], arg[1]))
print('    best refined sup estimate %.10f ; c_cert = -0.3478616532 ; slack %.3e' % (glob, -0.3478616532 - glob))
# Q2 (y >= U) and Q3 (x >= U) actual sups (no bound): f^(x)+f^(x+y)+A(x)-A(y)+beta(1-x)
yy = [float(U) + t for t in np.linspace(0, 6, 601)]
q2max = max(F1[i] + up(fh(arb(xs[i] + yv))) + Al[i] - float(A[-1]) + bf*(1 - xs[i]) for i in range(0, Nf+1, 8) for yv in yy[::20])
xx = [float(U) + t for t in np.linspace(0, 6, 301)]
Ag = np.array([float(a) for a in A]); 
q3max = max(up(fh(arb(xv))) + up(fh(arb(xv + float(t)))) + float(A[-1]) - float(Aval(Fr(t).limit_denominator(10**9))) + bf*(1 - xv)
            for xv in xx[::10] for t in list(np.linspace(float(lam), float(U), 78)) + [float(U) + 1, float(U) + 10])
print('    sampled sup on Q2 = %.6f (checker R2 = -0.5640124229), on Q3 = %.6f (checker R3 bound = -0.3481091233)' % (q2max, q3max))
# (d) Lamzouri's R^ is the FT of R (numerical)
import mpmath as mp
mp.mp.dps = 25
Rf = lambda x: (mp.sin(mp.pi*x)/(mp.pi*x))**2/(1 - x*x)
Rh = lambda u: 1 - abs(u) + mp.sin(2*mp.pi*abs(u))/(2*mp.pi)
for x in [0.3, 1.5, 2.7]:
    num = 2*mp.quad(lambda u: Rh(u)*mp.cos(2*mp.pi*u*x), [0, 1])
    print('(d) x = %.1f: FT of R^ = %s, R(x) = %s' % (x, mp.nstr(num, 12), mp.nstr(Rf(x), 12)))
