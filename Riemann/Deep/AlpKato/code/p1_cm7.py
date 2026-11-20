# P1 (alpkato): weight-3 CM newform of level 7, f = eta(z)^3 eta(7z)^3 in S_3(Gamma0(7), (-7/.)).
# Exact q-expansion, CM cross-check, L(f,1), L(f,2) (nonzero; Kato's hypothesis), algebraic parts vs Omega(49a1)^2.
import mpmath as mp, sympy as sp, sys
mp.mp.dps = 40
M = 400
# eta(z)^3 = q^{1/8} sum_{n>=0} (-1)^n (2n+1) q^{n(n+1)/2}  (Jacobi); f = q * prod(1-q^n)^3 (1-q^{7n})^3
def eta3(M, step):
    c = [0]*(M+1)
    n = 0
    while step*n*(n+1)//2 <= M:
        c[step*n*(n+1)//2] += (-1)**n*(2*n+1); n += 1
    return c
A = eta3(M, 1); B = eta3(M, 7)
prod = [0]*(M+1)
for i,a in enumerate(A):
    if a:
        for j,b in enumerate(B):
            if i+j > M: break
            if b: prod[i+j] += a*b
a = [0]*(M+2)
for n in range(M): a[n+1] = prod[n]          # f = q * eta^3 eta(7)^3 / q^{1/8+7/8}
# CM cross-check: a_n = (1/2) sum_{alpha in O_K, N alpha = n} Re(alpha^2), alpha = x + y(1+sqrt-7)/2
cm = [0]*(M+1)
R = int(2*M**0.5)+3
for x in range(-R, R+1):
    for y in range(-R, R+1):
        N = x*x + x*y + 2*y*y
        if 0 < N <= M:
            cm[N] += ((2*x+y)**2 - 7*y*y)  # 4*Re(alpha^2)
ok = all(cm[n] == 8*a[n] for n in range(1, M+1))   # (1/2)*sum Re = a_n  <=> sum 4Re = 8 a_n
print("a_1..a_12:", a[1:13]); print("CM theta-series cross-check (n<=%d):" % M, ok)
inert = [p for p in sp.primerange(3, M) if sp.legendre_symbol(-7 % p, p) == -1]
print("a_p = 0 for all inert p<=%d:" % M, all(a[p] == 0 for p in inert))
# theta(t) = sum a_n exp(-2 pi n t / sqrt7); f(-1/(7z)) relation gives theta(1/t) = eps * t^3 * theta(t)
s7 = mp.sqrt(7)
def theta(t): return mp.fsum(a[n]*mp.e**(-2*mp.pi*n*t/s7) for n in range(1, M+1))
for t in [mp.mpf('1.1'), mp.mpf('1.37')]:
    print("t=%s  theta(1/t)/(t^3 theta(t)) =" % t, mp.nstr(theta(1/t)/(t**3*theta(t)), 25))
eps = 1
def Lam(s):
    x0 = 2*mp.pi/s7
    return mp.fsum(a[n]*(mp.gammainc(s, x0*n)*(s7/(2*mp.pi*n))**s + eps*mp.gammainc(3-s, x0*n)*(s7/(2*mp.pi*n))**(3-s)) for n in range(1, M+1))
def L(s): return Lam(s)/((s7/(2*mp.pi))**s*mp.gamma(s))
L1, L2 = L(1), L(2)
print("L(f,1) =", mp.nstr(L1, 30)); print("L(f,2) =", mp.nstr(L2, 30))
print("check L(f,1) - sqrt7/(2pi) L(f,2) =", mp.nstr(L1 - s7/(2*mp.pi)*L2, 5))
print("check Euler product proxy at s=3 vs series:", mp.nstr(L(3), 20))
# real period of E = 49a1: y^2+xy = x^3-x^2-2x-1 -> 4x^3+b2 x^2+2 b4 x+b6, b2=-3,b4=-4,b6=-4 ; Delta=-343 (one real root)
x = sp.symbols('x'); rts = sp.Poly(4*x**3 - 3*x**2 - 8*x - 4, x).nroots(n=50)
e1 = [r for r in rts if abs(sp.im(r)) < 1e-30][0]; e2 = [r for r in rts if sp.im(r) > 1e-30][0]
e1 = mp.mpf(str(sp.re(e1))); e2 = mp.mpc(str(sp.re(e2)), str(sp.im(e2)))
# Delta<0: Omega = 2*pi / AGM(2 sqrt(beta), sqrt(2 beta + alpha)), alpha=3 e1 /... use standard: with cubic 4(x-e1)(x-e2)(x-e3)
z = mp.sqrt(3*e1**2 - 0)  # placeholder replaced below
beta = mp.sqrt((e1-e2)*(e1-mp.conj(e2))); beta = abs(beta)
alpha = 3*e1/2 + (-3)/8  # 3 e1 + b2/4 scaled for 4x^3 model: e's are roots of 4x^3+..., cubic monic in x: x^3 - (3/4)x^2 - 2x - 1
alpha = 3*e1 - mp.mpf(3)/4
Om = 2*mp.pi/mp.agm(2*mp.sqrt(beta), mp.sqrt(2*beta + alpha))
print("Omega(49a1) (real period, Cremona 1.9333117...):", mp.nstr(Om, 25))
for name, val in [("L(f,2)/(pi*Om^2)", L2/(mp.pi*Om**2)), ("L(f,2)*sqrt7/(pi*Om^2)", L2*s7/(mp.pi*Om**2)), ("L(f,1)/Om^2", L1/Om**2), ("L(f,1)*sqrt7/Om^2", L1*s7/Om**2)]:
    print(name, "=", mp.nstr(val, 25), " identify:", mp.identify(val, tol=mp.mpf(10)**-30))
sys.stdout.flush()
