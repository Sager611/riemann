# semi_pair2.py: copy of ../../ACLGFals/code/semi_pair.py (aclgfals) with (i) a tabulated kernel C_S (cost independent of the number
# of Laurent terms), (ii) the scale t, (iii) a merged eigen-list cut to its complete range, (iv) lower memory (one vector matrix).
# Model (../../ACLGFals/DERIVATION.md sec. 3): Gamma_t = (1 - P_{1/t}) U P_t is unitarily equivalent to the symmetric operator
#   (A G)(w) = 2 t int_0^1 C_S(t^2 u w) G(u) du on L^2(0,1),  C_S(z) = sum_n gamma(n) cos(2 pi n z),
#   n = prod_{p in S} p^{k_p} (k_p >= -1), gamma(n) = prod_p g_p(k_p), g_p(k) = 1 - 1/p (k >= 0), -1/p (k = -1), cut n <= nmax.
# Mellin data in the rescaled frame: |M a_n(i tau)|^2 = |M e~ + conj(phi_t) M g~|^2/(2 + 2 kappa), phi_t = t^{2 i tau} u_K(i tau),
#   u_K = exp(-2 i vartheta) sum_n gamma(n) n^{-1/2 + i tau} (the cut symbol, used consistently: f = U_K^* g).
# m_n(R) = (1/2pi) int chi_R |M a_n|^2 dtau;  T(n0, R) = sum_{n > n0} kappa_n m_n(R)  (ACLG eq. H2ii, l. 1257).
# Usage: semi_pair2.py primes(comma or -) nmax R1,R2,.. [npw] [t]   env: SP_EIGH=k (pairs per end, default 400), SP_F32=1,
#   SP_EXACT=1 (direct cosine sum, as semi_pair.py), SP_PPW (table points per top wavelength, default 32).
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'): os.environ[k] = '1'
import sys, json, math, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hankel import chi_cut
from scipy.special import loggamma
import scipy.linalg as sla
DT = np.float32 if os.environ.get('SP_F32', '0') == '1' else np.float64
EXACT = os.environ.get('SP_EXACT', '0') == '1'
PPW = float(os.environ.get('SP_PPW', '32'))

t0 = time.time()
primes = [] if sys.argv[1] == '-' else [int(q) for q in sys.argv[1].split(',')]
nmax = float(sys.argv[2]); Rs = [float(r) for r in sys.argv[3].split(',')]
npw = float(sys.argv[4]) if len(sys.argv) > 4 else 7.0
tsc = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
terms = {1.0: 1.0}
for p in primes:
    new = {}
    for n, gm in terms.items():
        k = -1
        while True:
            m = n*p**k
            if m > nmax: break
            new[m] = new.get(m, 0.0) + gm*((1 - 1/p) if k >= 0 else -1/p)
            k += 1
    terms = new
freqs = np.array(sorted(terms)); gams = np.array([terms[f] for f in freqs])
fe = freqs*tsc*tsc                      # frequencies of the rescaled kernel cos(2 pi n t^2 u w)
ftop = max(float(fe.max()), 1.0)
# grid: composite GL (q = 20) in s on [0,1], u = s^2, panels <= 3*7/npw wavelengths of cos(2 pi ftop s^2)
q = 20
gx, gw = np.polynomial.legendre.leggauss(q)
edges = [0.0]
while edges[-1] < 1.0:
    a = edges[-1]; wdt = min(0.05, 3.0*7.0/(npw*2*ftop*max(a, 1e-3)))
    edges.append(min(1.0, a + wdt))
e = np.array(edges); npan = e.size - 1
s = ((e[:-1, None] + e[1:, None])/2 + (e[1:, None] - e[:-1, None])/2*gx[None, :]).ravel()
ws = ((e[1:, None] - e[:-1, None])/2*gw[None, :]).ravel()
u = s*s; om = 2*s*ws
n = u.size
sq = np.sqrt(om)
# table of 2 t C_S(t^2 z), z in [0, 1 + pad], uniform step h; 8-point Lagrange interpolation (offsets -3..4)
h = 1.0/(PPW*ftop)
nz = int(math.ceil(1.0/h)) + 8
zt = (np.arange(nz) - 3)*h
tab = np.zeros(nz)
CH = max(64, (1 << 21)//max(1, fe.size))            # chunk rows: <= 2M table entries per temporary (memory cap)
for c0 in range(0, nz, CH):
    c1 = min(nz, c0 + CH)
    tab[c0:c1] = (2*tsc)*(np.cos(2*math.pi*np.outer(zt[c0:c1], fe)) @ gams)
OFF = np.arange(-3, 5)
DEN = np.array([np.prod([j - m for m in OFF if m != j]) for j in OFF], dtype=float)
def kern_block(ub, uall):
    z = np.outer(ub, uall)
    if EXACT:
        out = np.zeros_like(z)
        for f, gm in zip(fe, gams): out += (2*tsc*gm)*np.cos(2*math.pi*f*z)
        return out
    x = z/h; i0 = np.floor(x).astype(np.int64); fr = x - i0; i0 += 3   # table index of z_{i0} (offset 3)
    out = np.zeros_like(z)
    for jj, j in enumerate(OFF):
        L = np.ones_like(z)
        for m in OFF:
            if m != j: L *= (fr - m)
        out += tab[i0 + j]*(L/DEN[jj])
    return out
def build_K():
    K = np.empty((n, n), dtype=DT, order='F')
    for i0 in range(0, n, 128):
        i1 = min(n, i0 + 128)
        K[i0:i1] = kern_block(u[i0:i1], u)*(sq[i0:i1, None]*sq[None, :])
    return K
K = build_K()
tb = round(time.time() - t0, 1)
kE = min(int(os.environ.get('SP_EIGH', '400')), n//2)
mu1, v1 = sla.eigh(K, subset_by_index=[n - kE, n - 1], driver='evr', overwrite_a=True)
del K; K = build_K()
mu2, v2 = sla.eigh(K, subset_by_index=[0, kE - 1], driver='evr', overwrite_a=True)
del K
cut = max(abs(float(mu1[0])), abs(float(mu2[-1])))     # merged |mu| list is complete above this value
mu = np.concatenate([mu1, mu2]); o = np.argsort(-np.abs(mu)); o = o[np.abs(mu[o]) >= cut]
Gm = np.concatenate([v1, v2], axis=1)[:, o].astype(np.float64); del v1, v2
mu = mu[o].astype(np.float64); S = np.abs(mu); sgn = np.sign(mu)
Gm /= sq[:, None]                                       # G_n(u_j); H_n = sgn_n G_n (K symmetric)
res = {'info': {'primes': primes, 'nmax': nmax, 't': tsc, 'npw': npw, 'f32': DT == np.float32, 'exact': EXACT, 'ppw': PPW,
                'nterms': int(freqs.size), 'nodes': int(n), 'kE': kE, 'cut': cut, 'sec_kernel': tb, 'sec_eig': round(time.time() - t0, 1)},
       'kappa': [float(x) for x in S], 'mass': {}, 'check': {}}
for R in Rs:
    tp = max(4, int(math.ceil(2*R/1.5)))
    te = np.linspace(0, 2*R, tp + 1)
    tau = ((te[:-1, None] + te[1:, None])/2 + (te[1:, None] - te[:-1, None])/2*gx[None, :]).ravel()
    tw = ((te[1:, None] - te[:-1, None])/2*gw[None, :]).ravel()
    chi = chi_cut(tau, R)
    th = np.imag(loggamma(0.25 + 0.5j*tau)) - 0.5*tau*math.log(math.pi)
    uK = np.exp(-2j*th)*sum(gm*f**-0.5*np.exp(1j*tau*math.log(f)) for f, gm in zip(freqs, gams))
    ub = np.conj(uK*np.exp(2j*tau*math.log(tsc)))
    W = np.zeros((tau.size, n), dtype=complex)
    g16, w16 = np.polynomial.legendre.leggauss(16)
    for pi_ in range(npan):
        sa, sb = e[pi_], e[pi_ + 1]; nodes = s[pi_*q:(pi_ + 1)*q]
        ya, yb = math.log(max(sa, 1e-7)), math.log(sb)
        npy = max(1, int(math.ceil((yb - ya)*4*R/6.0)))
        ye = np.linspace(ya, yb, npy + 1)
        y = ((ye[:-1, None] + ye[1:, None])/2 + (ye[1:, None] - ye[:-1, None])/2*g16[None, :]).ravel()
        wy = ((ye[1:, None] - ye[:-1, None])/2*w16[None, :]).ravel()
        sf = np.exp(y); wf = wy*sf
        bw = np.array([1.0/np.prod([nodes[j] - nodes[k] for k in range(q) if k != j]) for j in range(q)])
        diff = sf[:, None] - nodes[None, :]
        diff[diff == 0] = 1e-300
        L = (bw[None, :]/diff)/np.sum(bw[None, :]/diff, axis=1, keepdims=True)
        LW = L*(2*wf)[:, None]
        for c0 in range(0, tau.size, 256):
            c1 = min(tau.size, c0 + 256)
            W[c0:c1, pi_*q:(pi_ + 1)*q] = np.exp(2j*np.outer(tau[c0:c1], y)) @ LW
    m = np.zeros(S.size); cross = np.zeros(S.size)
    TB = 200
    for c0 in range(0, tau.size, TB):
        c1 = min(tau.size, c0 + TB)
        Wr = np.ascontiguousarray(W[c0:c1].real); Wi = np.ascontiguousarray(W[c0:c1].imag)
        A1 = Wr @ Gm; A2 = Wi @ Gm; del Wr, Wi
        Me = A1 + 1j*A2; Mg = (A1 - 1j*A2)*sgn[None, :]; del A1, A2
        wt = (tw*chi)[c0:c1]
        m += (2/(2*math.pi))*(wt[None, :] @ (np.abs(Me + ub[c0:c1, None]*Mg)**2)).ravel()
        cross += (2/(2*math.pi))*(wt[None, :] @ np.real(np.conj(Me)*ub[c0:c1, None]*Mg)).ravel()
        del Me, Mg
    del W
    m /= (2 + 2*S)
    res['mass'][str(R)] = [float(x) for x in m]
    res['check'][str(R)] = {'cross_over_kappa_first5': [float(cross[i]/S[i]) for i in range(min(5, S.size))]}
res['info']['sec'] = round(time.time() - t0, 1)
print(json.dumps(res))
