# consts.py -- lamz2audit: constants of arXiv 2609.02882 (Lamzouri) recomputed; Lemma 3.2 quote in 2610.26841 tested.
# Run: uv run --no-project --with python-flint --with mpmath python3 consts.py > logs/consts.log
import sys, time
from flint import arb, acb, ctx
import mpmath as mp
ctx.prec = 256
t0 = time.time()
s2 = arb(2).sqrt()
def cot(x): return x.cos()/x.sin()
out = []
def P(*a):
    print(*a); sys.stdout.flush()

P("== 1. Closed forms (Arb, 256 bits)")
CMT = arb(1)/2 + cot(1/s2)/s2
C0 = arb(3)/2 - cot(1/s2)/s2
C1 = (C0+1)/2
C2 = (1 + 2*s2 + 2*C0)/(3 + 2*s2)
P("C_MT =", CMT.str(30))
P("C_0 = 2 - C_MT =", C0.str(30), " check 2-C_MT-C_0 contains 0:", (2-CMT-C0).contains(0))
P("C_1 = (C_0+1)/2 =", C1.str(30))
P("C_2 = (1+2sqrt2+2C_0)/(3+2sqrt2) =", C2.str(30))
P("C_2 via (5+2sqrt2-2C_MT)/(3+2sqrt2) =", ((5+2*s2-2*CMT)/(3+2*s2)).str(30))
P("(2+C_0)/3 (BGST25 Thm3(iii) b->0) =", ((2+C0)/3).str(20))
# paper digits: v2 l.588 '1.3274992963205...', v1 '1.3274992963206...'
for lab, s in [("v2 1.3274992963205", "1.3274992963205"), ("v1 1.3274992963206", "1.3274992963206")]:
    lo = arb(s); hi = arb(s) + arb("1e-13")
    P("paper digits", lab, ": C_MT in [digits, digits+1e-13) ?", bool(CMT > lo and CMT < hi))
P("abstract/thm: 67.25%, 83.62%, 88.76%: C_0>0.6725", bool(C0 > arb("0.6725")), " C_1>0.8362", bool(C1 > arb("0.8362")), " C_2>0.8876", bool(C2 > arb("0.8876")))

P("== 2. C(f_0) = int f0^2 + iint |u-v| f0 f0, three routes")
a = s2
A = 1/(s2*(1/s2).sin())          # f_0 = A cos(sqrt2 x) on [-1/2,1/2]
P("int f_0 = 1 (closed form):", (A*2*(a/2).sin()/a).str(25))
P("f_0 > 0 on I: cos(sqrt2/2) =", (a/2).cos().str(15), "> 0; sqrt2/2 < pi/2:", bool(a/2 < arb.pi()/2))
# route (a): Euler-Lagrange: C = 1/2 + f_0(1/2)
Ca = arb(1)/2 + A*(a/2).cos()
P("(a) 1/2 + f_0(1/2) =", Ca.str(30), " equals C_MT:", (Ca-CMT).contains(0))
# route (b): int f0^2 closed + Arb integral of f0(u)*I(u), I(u)=int |u-v| cos(a v) dv closed form
If2 = A*A*(arb(1)/2 + a.sin()/(2*a))
def integrand(u, analytic):
    au = a*u
    Iu = (a/2).sin()/a + 2*((a/2).cos() - au.cos())/(a*a)
    return A*A*au.cos()*Iu
J = acb.integral(integrand, -0.5, 0.5)
Cb = If2 + J.real
P("(b) int f0^2 =", If2.str(20), " iint|u-v|f0f0 =", J.real.str(20), " sum =", Cb.str(25), " - C_MT contains 0:", (Cb-CMT).contains(0))
# route (c): the paper's form Q_0(0) + 2 int_0^1 alpha Q_0(alpha) d alpha, Q_0 = f0*f0 by mpmath double quadrature
mp.mp.dps = 30
Am = 1/(mp.sqrt(2)*mp.sin(1/mp.sqrt(2)))
f0 = lambda x: Am*mp.cos(mp.sqrt(2)*x) if abs(x) <= 0.5 else mp.mpf(0)
def Q0(al):
    lo = max(-0.5, al-0.5); hi = min(0.5, al+0.5)
    if lo >= hi: return mp.mpf(0)
    return mp.quad(lambda u: f0(u)*f0(al-u), [lo, hi])
Cc = Q0(0) + 2*mp.quad(lambda al: al*Q0(al), [0, 1])
P("(c) Q_0(0) + 2 int alpha Q_0 =", mp.nstr(Cc, 22), " diff to C_MT:", mp.nstr(Cc - mp.mpf(CMT.str(40, radius=False)), 3))
P("Q_0(0) = int f0^2 check:", mp.nstr(Q0(0), 20), " vs", If2.str(20))

P("== 3. The quote of Lemma 3.2 in 2610.26841 (l.622-700): delta = eps, any admissible psi_eps, |C - C_MT| < eps ?")
P("Indicator cutoff f = 1_{|x|<=a} cos(sqrt2 x)/norm, a = 1/2 - d (L^2-limit of admissible smooth psi_d):")
P("  E-L gives C = a + cot(sqrt2 a)/sqrt2;  D(d) := C - C_MT;  D(d) - d has derivative -2 + csc^2(sqrt2(1/2-d)) >= -2 + csc^2(1/sqrt2) > 0, value 0 at d=0.")
csc2 = 1/((1/s2).sin()**2)
P("  csc^2(1/sqrt2) - 2 =", (csc2-2).str(20), " > 0:", bool(csc2 - 2 > 0))
for ds in ["0.1", "0.05", "0.01", "0.001", "0.000001"]:
    d = arb(ds); aa = arb(1)/2 - d
    C = aa + cot(s2*aa)/s2
    D = C - CMT
    P("  d =", ds, " D(d) =", D.str(12), " D/d =", (D/d).str(8), " D > d:", bool(D > d))
P("  limit D/d -> csc^2(1/sqrt2) - 1 =", (csc2-1).str(10))
P("t =", round(time.time()-t0, 1), "s")
