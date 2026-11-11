# Synthetic checks of the deterministic steps of Lamzouri 2611.14486 with OFF-LINE zeros (impossible with real data).
# N1: identity (IntegralExpressionGeneral, l.302-305; BGST Lemma 3) for invariant multisets; failure without invariance
#     or with the RH weight 4/(4+(gamma-gamma')^2).
# N5: inequalities of Prop SmallSetZeros (l.360-389): |w|<=4/3, grouping (factor 8), separation, both Hoelder steps;
#     exponent arithmetic (exact rationals).
import math, random
from fractions import Fraction as Fr
import numpy as np

def w(u): return 4/(4 - u*u)

def lhs(B, T, a, weight=w):
    return sum(T**((r - s)*a)*weight(r - s) for r in B for s in B)

def rhs(B, T, a, A=4000.0, n=1_600_001):
    t = np.linspace(-A, A, n); S = np.zeros_like(t, dtype=complex)
    for r in B:
        S += T**(a*(r - 0.5))/(1 - (r - (0.5 + 1j*t))**2)
    y = np.abs(S)**2; hstep = t[1] - t[0]
    simpson = hstep/3*(y[0] + y[-1] + 4*y[1:-1:2].sum() + 2*y[2:-1:2].sum())
    return 2/math.pi*simpson

def N1():
    B = [0.5+3.1j, 0.5+7.4j, 0.8+5.0j, 0.2+5.0j, 0.6+5.6j, 0.4+5.6j, 0.5+9.0j, 0.5+9.0j, 0.95+12.0j, 0.05+12.0j]
    T = 1000.0; worst = 0.0
    for a in (0.37, -0.81, 1.3):
        L = lhs(B, T, a); Rv = rhs(B, T, a)
        worst = max(worst, abs(L - Rv)/abs(Rv))
        print(f'N1 invariant B (3 off-line pairs, 1 double, 2 simple on-line), T={T:g}, alpha={a}: LHS={L.real:.10f}{L.imag:+.1e}i  RHS={Rv:.10f}  rel.err={abs(L-Rv)/abs(Rv):.1e}')
    Bn = [r for r in B if r != 0.2+5.0j]  # drop one partner: not invariant
    a = 0.37; L = lhs(Bn, T, a); Rv = rhs(Bn, T, a)
    print(f'N1 non-invariant B (partner of 0.8+5i dropped), alpha={a}: LHS={L.real:.6f}{L.imag:+.6f}i  RHS={Rv:.6f}  -> identity fails: {abs(L-Rv)/abs(Rv) > 1e-3}')
    Lw = lhs(B, T, a, weight=lambda u: 4/(4 + (u.imag)**2))
    print(f'N1 invariant B with RH weight 4/(4+(g-g\')^2): LHS={Lw.real:.6f}  RHS={rhs(B, T, a):.6f}  -> fails: {abs(Lw-rhs(B,T,a))/rhs(B,T,a) > 1e-3}')
    return worst

def N5(seed=1, th=0.55):
    random.seed(seed); T = 1e6; LT = math.log(T)
    E = []
    for _ in range(25):
        g = random.uniform(0, 30)
        if random.random() < 0.5:
            b = random.uniform(0.001, 0.49); E += [complex(0.5 + b, g), complex(0.5 - b, g)]
        else:
            E += [complex(0.5, g)]*random.choice((1, 1, 2))
    Z = E + [complex(0.5, random.uniform(0, 30)) for _ in range(200)]   # Z contains E (extra zeros on the line)
    K = lambda r, s: 1/(1 + ((r.imag - s.imag)*LT)**2)
    ok = True
    wmax = max(abs(w(r - s)) for r in E for s in E)
    wb = all(abs(w(r - s)) <= 4/(4 + (r.imag - s.imag)**2 - (r.real - s.real)**2) + 1e-12 for r in E for s in E)
    ok &= wmax <= 4/3 and wb
    A = sum(T**(th*abs(r.real - s.real))*K(r, s) for r in E for s in E)
    Eg = [r for r in E if r.real >= 0.5]
    Bv = 2*sum((T**(th*(r.real - s.real)) + T**(th*(1 - r.real - s.real)) + T**(th*(r.real + s.real - 1)) + T**(th*(s.real - r.real)))*K(r, s) for r in Eg for s in Eg)
    C8 = sum(T**(th*(r.real + s.real - 1))*K(r, s) for r in Eg for s in Eg)
    D = sum(T**(th*(2*r.real - 1))*sum(K(r, s) for s in Eg) for r in Eg)
    eps = (1 - 7*th/4)/2; q = math.floor(1/eps) + 1; p = q/(q - 1)
    H1 = sum(T**(th*p*(2*r.real - 1)) for r in Eg)**(1/p)*sum(sum(K(r, s) for s in Z if s.real >= 0.5)**q for r in Z if r.real >= 0.5)**(1/q)
    H2 = len(E)**(1/q)*sum(T**(th*p*p*(2*r.real - 1)) for r in Z if r.real >= 0.5)**(1/p)
    inner = sum(T**(th*p*(2*r.real - 1)) for r in Eg)
    ok &= A <= Bv*(1 + 1e-12) and Bv <= 8*C8*(1 + 1e-12) and C8 <= D*(1 + 1e-12) and D <= H1*(1 + 1e-12) and inner <= H2*(1 + 1e-12)
    print(f'N5 seed={seed} theta={th} (#E={len(E)}, T=1e6, q={q}): |w|max={wmax:.4f}<=4/3 & bound: {wb};  A={A:.4e} <= B={Bv:.4e} <= 8C={8*C8:.4e};  C={C8:.4e} <= D={D:.4e} <= Hoelder1={H1:.4e};  inner={inner:.4e} <= Hoelder2={H2:.4e}  -> {ok}')
    return ok

def exponents():
    for th in (Fr(11, 20), Fr(4, 7) - Fr(1, 10**5)):
        e2 = (Fr(4, 7) - th)/2; t2 = th + e2                      # Thm Main2: g supported in [-t2, t2]
        eps = (1 - Fr(7, 4)*t2)/2; q = math.floor(1/eps) + 1; p = Fr(q, q - 1)
        c_prop = Fr(1, 8) - Fr(7, 32)*t2; c_main = Fr(1, 16) - Fr(7, 64)*th
        print(f'exponents theta={float(th):.8f}: t2={float(t2):.8f} eps={float(eps):.3e} q={q} 2*t2*p^2={float(2*t2*p*p):.8f} (<8/7={8/7:.8f}: {2*t2*p*p < Fr(8,7)})  1/(pq)={float(1/(p*q)):.3e} >= eps/2={float(eps/2):.3e}: {1/(p*q) >= eps/2}  delta-exponent 1/8-7t2/32={float(c_prop):.4e} == 1/16-7th/64={float(c_main):.4e}: {c_prop == c_main}')

if __name__ == '__main__':
    worst = N1()
    print(f'N1 worst relative error (invariant): {worst:.1e}')
    print('N5 all ok:', all(N5(seed=s, th=th) for s in (1, 2, 3) for th in (0.52, 0.55, 0.56)))
    exponents()
