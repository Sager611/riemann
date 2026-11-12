
import sympy as sp
N, d, z = sp.symbols('N d z')
rel = d**3 + 3*N*d - 2*N
a = N + d; u = N - d/2; v = N**2 - d*N + d**2 + 3*N; P = N*(N+1)*(N+2)
def red(e):  # remainder of numerator modulo rel
    num, den = sp.fraction(sp.together(sp.expand(e)))
    q, r = sp.div(sp.expand(num), rel, d, N)
    return sp.simplify(r)
print('I1 av-P - rel', sp.expand(a*v - P - rel))
print('I2', sp.expand(2*u*a + v - 3*N*(N+1)), 'I3', sp.expand(a + 2*u - 3*N))
print('I4 S(a)', sp.factor(sp.expand(a**2 - 2*u*a + v)))
print('I5 v-a^2', sp.factor(sp.expand(v - a**2)), ' I6 v-u^2', sp.factor(sp.expand(v - u**2)))
R = v/(3*(d**2+N))
KEYd = sp.expand(v*(v-u**2) - 3*a**2*(d**2+N))
print('KEY', sp.factor(KEYd), '|', sp.expand(KEYd*sp.Rational(4,9)))
# factorization P q(z) = (a - z)(z^2 - 2uz + v)
Pq = P - 3*N*(N+1)*z + 3*N*z**2 - z**3
print('fact', sp.expand((a - z)*(z**2 - 2*u*z + v) - Pq - rel))
# q' at a
Pqd = sp.diff(Pq, z)
print('P q\'(a) + 3(d^2+N)', sp.expand(Pqd.subs(z, a) + 3*(d**2+N)))
# NOTE closed form vs ours
K_over_R2 = 3*a**2*(d**2+N)/(v*(v-u**2))
note = 8*N*(N+d)**2*(1-d)**2/(d**3*((N+d)**2+3*N*(1-d))*(3*N+d**2+d+1))
E = sp.expand((d**2+N)*d**3*(3*N+d**2+d+1) - 2*N*(1-d)**2*(d**2+4*N))
q, r = sp.div(E, rel, N, d)
print('E = Q*rel + r: Q =', sp.factor(q), ' r =', r)
q2, r2 = sp.div(E, rel, d, N)
print('order d: Q =', sp.expand(q2), ' r =', r2)
print('check v form', sp.expand(v - ((N+d)**2 + 3*N*(1-d))))
