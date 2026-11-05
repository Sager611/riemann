
# check_neg_IL.py: numeric |integrals| for I_L and the Neg sum vs the rigorous bounds (m = 203)
import sys, math
import mpmath as mp
from math import comb, factorial
import assembly as A
from flint import arb
mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 400
m = int(sys.argv[1]); ns = [int(v) for v in sys.argv[2].split(',')]
nu0 = arb('203.5'); A.load(nu0, 'tab203.json', 'lg203.pkl')
nu = mp.mpf(m) + mp.mpf(1)/2
pc = [mp.mpf(comb(m,j)*factorial(2*m-j))/factorial(2*m) for j in range(m+1)]
qc = [(-1)**j*pc[j] for j in range(m+1)]

def poly(c, z):
    s = mp.mpc(0)
    for a in reversed(c): s = s*z + a
    return s
Ec = [mp.mpf(factorial(m)*factorial(m+j))/(factorial(j)*factorial(2*m)*factorial(2*m+j+1)) for j in range(1100)]
th0 = mp.mpf(factorial(2*m))/(factorial(m)*2**m)
Lam = th0*mp.sqrt(nu)*nu**(-nu)
kap = mp.sqrt(2*mp.pi*nu)*nu**nu*mp.exp(-nu)/mp.gamma(nu+1)
yJ = mp.mpf('0.85')
def lse(v):
    mx = max(v); return mx + mp.log(sum(mp.exp(x-mx) for x in v))
for n in ns:
    tau = n/float(nu)
    # I_L numeric: (1/2pi) int_{pi/2}^{3pi/2} |1/q(2 nu zeta)| / Lam * yJ^{-n} dphi   (normalized like c_n*(2nu)^n/Lam)
    K = 200; vals = []
    for k in range(K):
        ph = mp.pi/2 + mp.pi*(k+0.5)/K; zeta = yJ*mp.expj(ph); z = 2*nu*zeta
        vals.append(-mp.log(abs(poly(qc, z))) - mp.log(Lam) - n*mp.log(yJ) + mp.log(mp.pi/K))
    ILnum = lse(vals) - mp.log(2*mp.pi)
    pcs = A.IL_pieces(tau, nu0)
    ILb = float(lse([float(nu)*float(R.mid()) + float(c.mid()) + math.log(float(ln.mid())) for (R, c, ln) in pcs]))
    # Neg numeric: sum over l >= LN of (1/2pi) int_{-yJ}^{yJ} |g_l(iy)| dy
    LN = int(math.floor(tau/2)) + 1
    negvals = []
    for l in range(LN, LN + 12):
        K2 = 160; v = []
        for k in range(K2):
            y = yJ*(k+0.5)/K2; zeta = 1j*y; z = 2*nu*zeta
            p = poly(pc, z); E = poly(Ec, z); X = z**(2*m+1)*E/p
            v.append(mp.log(abs(mp.exp(z)/p)) + l*mp.log(abs(X)) - (n+1)*mp.log(y) - mp.log(Lam) - l*mp.log(kap) + mp.log(yJ/K2))
        negvals.append(lse(v) - mp.log(2*mp.pi) + mp.log(2))
    NEGnum = lse(negvals)
    npcs = A.neg_pieces(LN, tau, nu0)
    NEGb = float(lse([float(nu)*float(R.mid()) + float(c.mid()) + math.log(float(ln.mid())) for (R, c, ln) in npcs]))
    print('n=%d tau=%.3f: I_L numeric %.3f bound %.3f %s | Neg numeric %.3f bound %.3f %s' % (n, tau, float(ILnum), ILb, 'OK' if ILb >= ILnum else 'VIOL', float(NEGnum), NEGb, 'OK' if NEGb >= NEGnum else 'VIOL'), flush=True)
