# check_entries.py (explore19): DH/chi5 closed-form entries vs x-space quadrature of the archimedean part (independent of the
# digamma closed forms):  A(f1,f2) = (log(q/pi) - gamma) <f1,f2> + int_0^inf [2 e^{-x} <f1,f2> - 2 e^{(1-2 z0) x} C12(x)]/(e^x - e^{-x}) dx
# (derivation: psi(z) = -gamma + int_0^inf (e^{-u} - e^{-zu})/(1 - e^{-u}) du, u = 2x), plus -2 sum wt(n) C12(log n) (direct).
import sys, math
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx, fmpq
import mpmath as mp
import wmat
ctx.prec = 200; mp.mp.dps = 40
c = 31; cq = fmpq(c)
F = wmat.Form(cq, 'dh')
T, idx = F.block(8, 'even')
L = mp.log(c); z0 = mp.mpf(3)/4; q = 5
def w(j): return 2*mp.pi*j/L
def Nj(j): return mp.sqrt(1/L) if j == 0 else mp.sqrt(2/L)
def C(i, j, x):   # autocorrelation of phi_i, phi_j at shift x (0 <= x <= L), by quadrature
    if x >= L: return mp.mpf(0)
    a, b = -L/2 + x, L/2
    return mp.quad(lambda y: Nj(i)*mp.cos(w(i)*y)*Nj(j)*mp.cos(w(j)*(y - x)), [a, b])
def entry(i, j):
    ip = 1 if i == j else 0
    f = lambda x: (2*mp.exp(-x)*ip - 2*mp.exp((1 - 2*z0)*x)*C(i, j, x))/(mp.exp(x) - mp.exp(-x))
    A = (mp.log(q/mp.pi) - mp.euler)*ip + mp.quad(f, [0, 0.5, 1, 2, L, mp.inf])
    a_, b_, k = wmat.dh_coeffs(c + 1)
    P = mp.mpf(0)
    for n in range(2, c):
        P += mp.mpf(b_[n].mid().str(45, radius=False))/mp.sqrt(n)*C(i, j, mp.log(n))
    return A - 2*P
for (i, j) in [(0, 0), (1, 1), (0, 2), (3, 5), (7, 7)]:
    e1 = entry(i, j); e2 = mp.mpf(T[i, j].mid().str(45, radius=False))
    print(i, j, mp.nstr(e1, 20), mp.nstr(e2, 20), 'diff', mp.nstr(e1 - e2, 3), flush=True)
