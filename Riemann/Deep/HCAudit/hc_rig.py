# hcaudit: rigorous re-certification of cert768 by a second route (EPS = 1/1000): M2 by Arb ball enclosures of f^'' over
# [lam, 2U] (no Cramer/M3), corner maxima on a 8x finer lattice (dg/40) fully in Arb, J by rigorous acb integration.
import sys, json, os
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq, arb, acb, arb_poly, acb_poly, ctx
d = json.load(open(sys.argv[1])); ctx.prec = 128
lam = Fr(d['lam']); th = Fr(d['theta']); D = d['D']; s = Fr(d['s']); U = Fr(d['U']); dg = Fr(d['dg'])
B = [Fr(b) for b in d['B']]; B[0] += Fr(os.environ.get('EPS', '0')); A = [Fr(a) for a in d['A']]; beta = Fr(d['beta'])
q2 = lambda x: fmpq(x.numerator, x.denominator)
H = [fmpq_poly([1]), fmpq_poly([0, 2])]
for n in range(1, 2*D): H.append(fmpq_poly([0, 2])*H[n] - 2*n*H[n-1])
P = sum((q2(B[n])*H[2*n] for n in range(D+1)), fmpq_poly([0])); Q = sum((q2(B[n]*(-1)**n)*H[2*n] for n in range(D+1)), fmpq_poly([0]))
pi = arb.pi(); sA = arb(q2(s)); k_ = (2*pi).sqrt()*sA
QA = arb_poly([arb(Q[i]) for i in range(Q.degree()+1)]); Q1 = QA.derivative(); Q2 = Q1.derivative()
def fh(u): v = k_*u; return sA*QA(v)*(-(v*v)/2).exp()
def fh2(u): v = k_*u; return sA*k_*k_*(Q2(v) - 2*v*Q1(v) + (v*v - 1)*QA(v))*(-(v*v)/2).exp()
# M2 by ball cover of [lam, 2U]
r = Fr(1, 40000); n = int((2*U - lam)/(2*r)) + 1; M2 = arb(0)
for i in range(n):
    ball = arb(q2(lam + (2*i + 1)*r), float(r)*1.000001)
    v = abs(fh2(ball)).upper()
    if v > M2: M2 = v
print('(r1) rigorous M2 = sup_{[lam,2U]}|f^\'\'| <= %s (%d balls of radius %s); checker used 4.553' % (M2.str(8), n, r))
# corner maxima on lattice dg/40 entirely in Arb
m = 40; de = dg/m; N = 77*m
nodes = [lam + j*dg for j in range(len(A))]
Al = [arb(q2(A[i//m] + (A[min(i//m+1, len(A)-1)] - A[i//m])*Fr(i % m, m))) for i in range(N+1)]
F1 = [fh(arb(q2(lam + i*de))) for i in range(N+1)]; F2 = [fh(arb(q2(2*lam + j*de))) for j in range(2*N+1)]
lin = [Al[i] + arb(q2(beta))*(1 - arb(q2(lam + i*de))) for i in range(N+1)]
best = None
for i in range(N+1):
    base = F1[i] + lin[i]
    for j in range(N+1):
        v = (base + F2[i+j] - Al[j]).upper()
        if best is None or v > best: best = v
rem = 3*M2*arb(q2(de))**2/8
R1 = best + rem
R2 = max(((F1[i] + lin[i] - Al[-1]).upper() for i in range(N+1)), key=lambda z: float(z)) + M2*arb(q2(de))**2/8
R3 = arb(q2(A[-1] - min(A) + beta*(1 - U)))
print('(r2) dg/40 lattice, all Arb: corner max %s + rem %s = R1 %s ; R2 %s ; R3 %s' % (best.str(12), rem.str(3), R1.str(12), R2.str(12), R3.str(12)))
c2 = max([R1, R2, R3], key=lambda z: float(z.upper())).upper()
# J by rigorous integration in acb
PA = acb_poly([acb(arb(P[i])) for i in range(P.degree()+1)]); sq = (2*pi).sqrt()/sA
def integrand(t, analytic):
    y = sq*t; return PA(y)*(-(y*y)/2).exp()*t
J = PA(acb(0)).real + 2*acb.integral(integrand, 0, acb(arb(q2(th)))).real
fh0 = sA*arb(Q(fmpq(0)))
print('(r3) J (acb.integral) = %s ; f^(0) = %s' % (J.str(15), fh0.str(12)))
mg = J - fh0 - 2*c2
print('(r4) second-route c = %s ; margin J - f^(0) - 2c = %s -> %s' % (c2.str(12), mg.str(10), 'PASS' if mg > 0 else 'FAIL'))
