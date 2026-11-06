# gw.py (explore17): float64 Galerkin pieces of the Weil form on W = [-L/2, L/2], L = ln c (same basis/conventions as
# Riemann/Deep/WeilSimple/weil_matrix.py: even phi_j = N_j cos(w_j x), j = 0..N; odd psi_j = sqrt(2/L) sin(w_j x), j = 1..N; w_j = 2 pi j/L).
#   Q(f) = V(L)||f||^2 + int_0^L (C_f(0) - C_f(x)) k(x) dx + 2<f,cosh(x/2)>^2 - 2<f,sinh(x/2)>^2 - 2 sum_{n<c} w_n C_f(log n),
#   k(x) = e^{x/2}/sinh x, V(L) = psi(1/4) - log pi + int_L^inf k,  w_n = Lambda(n)/sqrt n,  C_f(x) = int f(y) f(y - x) dy.
# Pieces: LK = int_0^L (I - C(x)) k(x) dx  (Kahane / Levy-Khintchine positive piece, PSD on R);  LOG = log|t| multiplier
#   = int_0^L (I - C(x))/x dx - (gamma + log L) I  (exactly scale invariant);  POLE (rank 2);  PR = prime atoms.
import numpy as np, math
EG = 0.5772156649015329
PSI14 = -EG - math.pi/2 - 3*math.log(2)      # psi(1/4)
PSI34 = -EG + math.pi/2 - 3*math.log(2)      # psi(3/4)

def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0

class Basis:
    def __init__(s, c, N, sector):
        s.c = c; s.L = math.log(c); s.N = N; s.sector = sector
        s.idx = np.arange(0, N + 1) if sector == 'even' else np.arange(1, N + 1)
        s.w = 2*math.pi*s.idx/s.L
        s.Nj = np.where(s.idx == 0, 1/math.sqrt(s.L), math.sqrt(2/s.L))
    def C(s, x):   # matrix C_ij(x), 0 <= x <= L (symmetric)
        L = s.L; w = s.w; n = len(w); x = float(x)
        I, J = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
        wi, wj = w[I], w[J]; sg = np.where((s.idx[I] + s.idx[J]) % 2 == 0, 1.0, -1.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            if s.sector == 'even':
                off = s.Nj[I]*s.Nj[J]*sg*(wj*np.sin(wj*x) - wi*np.sin(wi*x))/(wi*wi - wj*wj)
            else:
                off = (2/L)*sg*(wi*np.sin(wj*x) - wj*np.sin(wi*x))/(wi*wi - wj*wj)
        M = np.where(I == J, 0.0, off)
        d = np.empty(n)
        for a in range(n):
            if s.idx[a] == 0: d[a] = 1 - x/L
            else:
                ww = w[a]; d[a] = (1 - x/L)*math.cos(ww*x) + (-1 if s.sector == 'even' else 1)*math.sin(ww*x)/(ww*L)
        M[np.arange(n), np.arange(n)] = d
        return M
    def jump(s, kfun, nq=800):   # int_0^L (I - C(x)) k(x) dx by Gauss-Legendre
        xg, wg = np.polynomial.legendre.leggauss(nq); x = (xg + 1)*s.L/2; wq = wg*s.L/2
        n = len(s.idx); M = np.zeros((n, n)); I = np.eye(n)
        for xx, ww in zip(x, wq): M += ww*kfun(xx)*(I - s.C(xx))
        return M
    def pole(s):
        sL4 = math.sinh(s.L/4); sg = np.where(s.idx % 2 == 0, 1.0, -1.0)
        if s.sector == 'even':
            v = s.Nj*sg*sL4/(s.w**2 + 0.25); return 2*np.outer(v, v)
        v = -2*s.Nj*sg*s.w*sL4/(s.w**2 + 0.25); return -2*np.outer(v, v)
    def primes(s, weights=None):   # -2 sum_{n<c} w_n C(log n); weights: dict n -> w_n (default Lambda(n)/sqrt n)
        n0 = len(s.idx); M = np.zeros((n0, n0))
        for n in range(2, int(math.ceil(s.c))):
            if n >= s.c: break
            wn = (mangoldt(n)/math.sqrt(n)) if weights is None else weights.get(n, 0.0)
            if wn != 0.0: M += -2*wn*s.C(math.log(n))
        return M

def kz(x):  return math.exp(x/2)/math.sinh(x)      # zeta (even gamma factor) jump density
def kchi_odd(x): return math.exp(-x/2)/math.sinh(x)   # odd character gamma factor
def tail(L, a0):   # int_L^inf of 2 sum_k e^{-(2k+a0)x}: = 2 sum_k e^{-(2k+a0)L}/(2k+a0)
    t = 0.0; k = 0
    while True:
        b = 2*k + a0; term = 2*math.exp(-b*L)/b; t += term; k += 1
        if term < 1e-18: return t
def Vz(L):  return PSI14 - math.log(math.pi) + tail(L, 0.5)
def Vdh(L, q=5): return math.log(q/math.pi) + PSI34 + tail(L, 1.5)

def weil_parts(c, N, sector, nq=800):
    B = Basis(c, N, sector); n = len(B.idx)
    LK = B.jump(kz, nq); LOG = B.jump(lambda x: 1/x, nq) - (EG + math.log(B.L))*np.eye(n)
    return dict(B=B, LK=LK, V=Vz(B.L)*np.eye(n), LOG=LOG, POLE=B.pole(), PR=B.primes())

def dh_coeffs(nmax):
    kap = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
    a = [0.0]*(nmax + 1)
    for n in range(1, nmax + 1): a[n] = {1: 1.0, 2: kap, 3: -kap, 4: -1.0, 0: 0.0}[n % 5]
    b = [0.0]*(nmax + 1)   # -f'/f = sum b_n n^-s : sum_{d|n} a_d b_{n/d} = a_n log n
    for n in range(2, nmax + 1):
        t = a[n]*math.log(n)
        for d in range(2, n + 1):
            if n % d == 0: t -= a[d]*b[n//d]
        b[n] = t
    return kap, a, b
