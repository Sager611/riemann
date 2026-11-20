# latrough.py (wsgsaudit): the lattice identity (Step 2) in the rough regime the L^2 claim needs: f = phi_0 = L^{-1/2} 1_I (jumps at the
# endpoints, F(r_k) ~ 1/k). Lattice sum with my own mu(k) to K plus the smooth tail (4 N0^2 L^2/pi^2)(log(K/2L) + 1)/K vs weil_matrix T[0,0].
import sys
sys.path.insert(0, '/tmp/wsgs_copy')
import mpmath as mp
from flint import fmpq, ctx
from weil_matrix import WeilMat
mp.mp.dps = 20; cn = int(sys.argv[1]); K = int(sys.argv[2]); c = mp.mpf(cn); L = mp.log(c); N0 = 1/mp.sqrt(L)
def vM(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return mp.log(p) if m == 1 else mp.mpf(0)
atoms = [(mp.log(n), vM(n)/mp.sqrt(n)) for n in range(2, cn) if vM(n) != 0]
nb = int(60/float(mp.log10(c))) + 3
def mu(k):
    r = mp.pi*k/L; v = -mp.log(mp.pi) + mp.re(mp.digamma(mp.mpf(1)/4 + 0.5j*r))
    v += 2*(-1)**k*sum(c**(-(2*n + mp.mpf(1)/2))*(2*n + mp.mpf(1)/2)/((2*n + mp.mpf(1)/2)**2 + r*r) for n in range(nb))
    return v - 2*sum(w*mp.cos(r*ln) for ln, w in atoms)
s = mu(0)*(N0*L)**2
for k in range(1, K + 1, 2): s += 2*mu(k)*(2*N0*L/(mp.pi*k))**2
tail = 4*N0**2*L**2/mp.pi**2*(mp.log(K/(2*L)) + 1)/K   # sum over odd k > K of 2 mu(k) F_k^2, mu ~ log(k/2L), density 1/2
lat = 2*(N0*4*mp.sinh(L/4))**2 + (s + tail)/(2*L)
ctx.prec = 128; T, _ = WeilMat(fmpq(cn)).block(0, 'even')
print(f"c={cn} K={K} lattice {mp.nstr(lat, 12)}  matrix {T[0,0].mid().str(12, radius=False)}  diff {mp.nstr(lat - mp.mpf(T[0,0].mid().str(25, radius=False)), 3)}  (tail term {mp.nstr(tail/(2*L), 3)})")
