# pyright: reportMissingImports=false
"""Second check (Prop. 11.12, eq. (61)): R(100) = sum_s |log Theta_s(100)| / mu_1(100) <= 0.149 claimed.
Telescoping (53): mu_r = mu_1 + sum_{s<r} log Theta_s, so sum_{s<r} |log Theta_s| >= mu_1 - mu_r, and
mu_r(100) -> mu_inf = 2 log(t_101/t_100) (paper Cor. 11.2; unconditional here: first 101 zeros simple, on the line).
Hence R(100) >= 1 - mu_inf/mu_1(100)."""
import sys
sys.argv = ['x', '102']
src = open('/Users/adriansager/playground/riemann/.scratch/gershon/arb_c1.py').read().split('g = [gamma_k')[0]
exec(src)
from flint import acb, arb
g = {k: gamma_k(k) for k in (99, 100, 101)}
mu1 = (g[100]**2/(g[99]*g[101])).log()
print("mu_1(100) =", mu1.str(15), " a(100) = 100*eps =", (100*(1-g[99]*g[101]/g[100]**2)).str(10))
try:
    z = acb.zeta_zeros(100, 2)
    t100, t101 = z[0].imag, z[1].imag
except Exception as e:
    print("zeta_zeros unavailable:", e); t100, t101 = arb("236.524229665816205802880"), arb("237.769820480925204349")
print("t_100 =", t100.str(15), " t_101 =", t101.str(15))
muinf = 2*(t101/t100).log()
print("mu_inf =", muinf.str(12), " mu_inf/mu_1 =", (muinf/mu1).str(10))
print("R(100) >= 1 - mu_inf/mu_1 =", (1 - muinf/mu1).str(10), " paper claims R(100) < 0.149")
print("certified mu_inf < (4/5) mu_1:", bool(muinf < arb(4)/5*mu1))
