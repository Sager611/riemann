# tcheck.py (wsgsaudit): CCM's archimedean term itself, (1/2pi) int phi(t) |F(t)|^2 dt with 2 theta'(t) = phi(t) = -log pi + Re psi(1/4+it/2)
# (CCM eq. (3.10), l. 535-538, 455-461), by t-quadrature for phi_0 = L^{-1/2} on I, vs weil_matrix's closed form arch_diag(0) (Gauss/digamma route).
import sys
sys.path.insert(0, '/tmp/wsgs_copy')
import mpmath as mp
from flint import fmpq, ctx
from weil_matrix import WeilMat
mp.mp.dps = 20; cn = int(sys.argv[1]); L = mp.log(cn); N0 = 1/mp.sqrt(L); T = mp.mpf(4000)
phi = lambda t: -mp.log(mp.pi) + mp.re(mp.digamma(mp.mpf(1)/4 + 0.5j*t))
F2 = lambda t: (N0*L)**2 if t == 0 else (N0*2*mp.sin(t*L/2)/t)**2
pts = [mp.mpf(0)] + [2*mp.pi*m/L for m in range(1, int(T*L/(2*mp.pi)) + 1)]
I = sum(mp.quad(lambda t: phi(t)*F2(t), [pts[i], pts[i+1]]) for i in range(len(pts) - 1))
Te = pts[-1]
tail = 2*N0**2*((mp.log(Te/(2*mp.pi)) + 1)/Te)          # int_Te^inf phi * 4N0^2 sin^2/t^2 ~ 2N0^2 int log(t/2pi)/t^2 (oscillatory part O(1/T^2))
val = (I + tail)/mp.pi                                    # (1/2pi) int_R = (1/pi) int_0^inf (even integrand)
ctx.prec = 128; cf = WeilMat(fmpq(cn)).arch_diag(0, 'even')
print(f"c={cn} t-quadrature {mp.nstr(val, 12)}  closed form {cf.mid().str(12, radius=False)}  diff {mp.nstr(val - mp.mpf(cf.mid().str(30, radius=False)), 3)} (T={mp.nstr(Te, 6)})")
