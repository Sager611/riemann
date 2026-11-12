# Numerical check of the exact formula (10) of the paper (Gaussian-weighted twisted second moment).
# LHS = int zeta(1/2+a+it) zeta(1/2+b-it) (h/k)^{it} e^{-(t-T)^2/H^2} dt
# RHS = 2 pi sqrt(k/h) sum_n sigma_{a,-b}(n) K_b(n k/h) - 2 pi zeta(a+b) G(1-a) (h/k)^{1/2-a}
# K_b(x) = (2 x^{-b}/sqrt(pi)) int e^{-y^2+(1-2b+2iT)y/H} cos(2 pi x e^{2y/H}) dy, G(s) = exp((s-1/2-iT)^2/H^2)
import sys
from mpmath import mp, mpf, mpc, zeta, quad, exp, cos, sqrt, pi, linspace, log
mp.dps = 18
T, H = mpf(40), mpf(6)
a, b = mpc('0.03', '0.02'), mpc('0.05', '-0.01')
h, k = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1, 1)
r = mpf(h) / k
lhs = quad(lambda t: zeta(mpf(1)/2 + a + 1j*t) * zeta(mpf(1)/2 + b - 1j*t) * exp(1j*t*log(r)) * exp(-(t-T)**2/H**2),
           linspace(T - 6.5*H, T + 6.5*H, 40))
def sig(n):
    return sum(mpf(d)**(-a) * mpf(n//d)**(b) for d in range(1, n+1) if n % d == 0)   # sigma_{a,-b}(n) = sum_{uv=n} u^{-a} v^{b}
def Kb(x):
    f = lambda y: exp(-y**2 + (1 - 2*b + 2j*T)*y/H) * cos(2*pi*x*exp(2*y/H))
    return 2 * x**(-b) / sqrt(pi) * quad(f, linspace(-7, 7, 140))
G = lambda s: exp((s - mpf(1)/2 - 1j*T)**2 / H**2)
S = mpc(0); nmax = int(5 * T * h / (2*pi*k)) + 5
for n in range(1, nmax + 1):
    S += sig(n) * Kb(mpf(n) * k / h)
rhs = 2*pi*sqrt(mpf(k)/h) * S - 2*pi*zeta(a + b) * G(1 - a) * r**(mpf(1)/2 - a)
print("h,k =", h, k, " nmax =", nmax)
print("LHS =", lhs); print("RHS =", rhs); print("|LHS-RHS|/|LHS| =", abs(lhs - rhs)/abs(lhs))
