# B1/B2 numbers. B2: Nevanlinna half-plane characteristic of f(z) = F(1/2 - i z) on Im z > 0 (= right half-plane sigma > 1/2):
#   S(r) = A(r) + B(r) + C(r,inf), A(r) = (1/pi) int_1^r (1/t^2 - 1/r^2) [log+|f(t)| + log+|f(-t)|] dt  -> A(inf) = (2/pi) int_1^inf log+|F(1/2+it)| t^-2 dt,
#   B(r) <= (4/(pi r)) log F(2)-type bound -> 0, no poles with |z| > 1.  Tail t > T0 bounded with |zeta(1/2+it)| <= 0.618 t^(1/6) log t (Hiary 2016, t >= 3).
import numpy as np, math, sys
from fz import Fn, line_eval
T0 = 3000.0; dt = 0.002
for kind in ('zeta', 'dh'):
    f = Fn(kind, T0 + 20); n = int((T0 - 1)/dt) + 1; t = 1 + dt*np.arange(n)
    F, _ = line_eval(f, 0.5, 1.0, n, dt, Bk=50); lp = np.maximum(np.log(np.abs(F)), 0)
    A = (2/math.pi)*np.trapezoid(lp/t**2, t)
    tail = (2/math.pi)*((1/6)*(math.log(T0) + 1) + math.log(math.log(T0)) + 1/math.log(T0) + max(0, math.log(0.618)))/T0
    print(f"{kind}: A(inf) partial (t<= {T0:.0f}) = {A:.5f}; tail <= {tail:.5f} ({'Hiary' if kind=='zeta' else 'Hiary-type, heuristic for DH'})", flush=True)
# B1: sum over zeros 1/gamma^2 vs Hadamard (unconditional) sum_rho 1/(rho(1-rho)) = 2 + gamma - log(4 pi)
import mpmath as mp
mp.mp.dps = 20; g = [float(mp.im(mp.zetazero(k))) for k in range(1, 101)]
s2 = sum(1/x**2 for x in g); s4 = sum(1/x**4 for x in g); had = float(1 + mp.euler/2 - mp.log(4*mp.pi)/2)
tl = lambda T: (math.log(T/(2*math.pi)) + 1)/(math.pi*T)
print(f"sum_{{n<=100}} 1/gamma^2 = {s2:.7f}; tail bound (1/pi)(log(T/2pi)+1)/T at T=gamma_100: {tl(g[-1]):.6f}; Hadamard sum_(gamma>0) 1/(1/4+gamma^2) = {had:.9f} (RH to 3e12)")
print(f"sum 1/gamma^4 (n<=100) = {s4:.3e}; sum_(gamma>3e12) 1/gamma^2 <= {tl(3e12):.2e}")
