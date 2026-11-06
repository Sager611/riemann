
import sympy as sp
y,A,B,C,D,dd = sp.symbols('y A B C D d')
# generic polynomial p of degree 7 with symbolic coefficients, check identities for V_m and Jacobi recurrence
deg = 7
cs = sp.symbols('c0:%d'%(deg+1))
p = sum(cs[i]*y**i for i in range(deg+1))
d = 7
def E(f): return sp.expand(y*sp.diff(f,y))
def T(k): return y**k*sp.diff(p,y,k)   # times p
# V_m * A * p = y^m * (E(E-d)(E+A) p)^{(m)}
EEE = E(E(E(p))) + (A-d)*E(E(p)) - A*d*E(p)   # E(E-d)(E+A) = E^3 + (A-d)E^2 - A d E
ok1 = True
for m in range(0,deg+1):
    lhs = sp.expand(y**m*sp.diff(EEE,y,m))
    beta = A*(2*m+1-d) - d*(2*m+1) + 3*m**2 + 3*m + 1
    gam = m*(m-d)*(m+A)
    rhs = sp.expand(T(m+3) + (A-d+3+3*m)*T(m+2) + beta*T(m+1) + gam*T(m))
    if sp.simplify(lhs-rhs) != 0:
        ok1 = False; print('V_m mismatch at m', m)
print('V_m identity:', ok1)
# Jacobi operator
W = sp.expand(y*(1-y/A)*sp.diff(p,y,2) + (B - y + (d-1)*y/A)*sp.diff(p,y) + d*p)
ok2 = True
for m in range(0,deg-1):
    lhs = sp.expand(y*y**m*sp.diff(W,y,m))
    rhs = sp.expand((1-y/A)*T(m+2) + (B+m-y+(d-1-2*m)*y/A)*T(m+1) + (d-m)*(1+sp.Rational(m)/A)*y*T(m))
    if sp.simplify(lhs-rhs) != 0:
        ok2 = False; print('Jacobi rec mismatch at m', m, sp.factor(sp.simplify(lhs-rhs)))
print('Jacobi recurrence identity:', ok2)
# perturbed Jacobi factorization for p_F = 3F2(-d,A,C;B,D;D y/(A C)), numeric rational params
import random
for trial in range(3):
    Av, Bv, Cv, Dv = [sp.Rational(random.randint(50,400), random.randint(1,7)) for _ in range(4)]
    dv = random.randint(3,9)
    pF = sum((-1)**j*sp.binomial(dv,j)*sp.rf(Av,j)*sp.rf(Cv,j)*Dv**j/(Av**j*Cv**j*sp.rf(Bv,j)*sp.rf(Dv,j))*y**j for j in range(dv+1))
    Wv = sp.expand(y*(1-y/Av)*sp.diff(pF,y,2) + (Bv - y + (dv-1)*y/Av)*sp.diff(pF,y) + dv*pF)
    lhs = sp.expand(E(Wv) + Dv*Wv)
    epsp = (Cv-Dv)/Cv
    EE3 = E(E(E(pF))) + (Av-dv)*E(E(pF)) - Av*dv*E(pF)
    rhs = sp.expand(-epsp/Av*EE3)
    print('factorization trial', trial, sp.simplify(lhs-rhs)==0)
    # also check the 3F2 ODE
    ode = E(E(E(pF))) + (Bv-1+Dv-1)*E(E(pF)) + (Bv-1)*(Dv-1)*E(pF) - Dv*y/(Av*Cv)*sp.expand( E(E(E(pF))) + (Av+Cv-dv)*E(E(pF)) + (Av*Cv - dv*Av - dv*Cv)*E(pF) - dv*Av*Cv*pF )
    print('  3F2 ODE holds:', sp.simplify(sp.expand(ode))==0)
