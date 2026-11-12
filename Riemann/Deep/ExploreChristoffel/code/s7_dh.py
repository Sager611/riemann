# explore38 s7: generic Burnol/Christoffel bound (on-line zeros + interior vectors for off-line zeros, poles of V in Re s < 1/2),
# validated on zeta against s3; applied to Davenport-Heilbronn (DH) with its own NB Gram (digamma, copied from ExploreHankel/nbgram.py).
# usage: uv run --no-project --with numpy --with scipy --with mpmath python3 s7_dh.py NMAX_DH K_ON
import sys, math, time, json, resource
import numpy as np
import mpmath as mp
from math import gcd
from scipy.special import digamma
mp.mp.dps = 30
t0 = time.time()
NMAX = int(sys.argv[1]); KON = int(sys.argv[2])
OUT = '/Users/adriansager/playground/riemann/.scratch/explore38/'
EG = 0.57721566490153286061

def bound(N, pts, Vv, dphi_on, online_mask, p, R, That, tau, ret_parts=False):
    """pts: complex w (on-line: Re = 1/2 exactly, limit; interior: Re < 1/2). Vv = V(w). dphi_on = Re V'/V at on-line points.
    p, R: real poles of V in Re s < 1/2 and residues; That = (T, psi_w) values; tau: residues of V*T^ at p."""
    lam = 1.0 / N; L = math.log(N); n = len(pts)
    w = pts; wc = np.conj(w)
    q = R[None, :] / (p[None, :] - w[:, None])          # [a,k] coefficient of t^{-p_k} in Vpsi_a
    H = np.empty((n, n), dtype=np.complex128)
    E3 = lam ** (1 - p[:, None] - p[None, :]) / (1 - p[:, None] - p[None, :])   # [k,l]
    qE = q @ E3                                           # [a,l]
    for r0 in range(0, n, 128):
        r1 = min(n, r0 + 128)
        den = 1 - w[r0:r1, None] - wc[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            blk = (1 - Vv[r0:r1, None] * np.conj(Vv)[None, :] * lam ** den) / den
        for i in range(r0, r1):
            if online_mask[i]:
                blk[i - r0, i] = L + dphi_on[i]
        # pole cross terms
        P1 = lam ** (1 - w[r0:r1, None] - p[None, :]) / (1 - w[r0:r1, None] - p[None, :])          # [a,k]
        blk -= Vv[r0:r1, None] * (P1 @ np.conj(q).T)
        P2 = lam ** (1 - p[:, None] - wc[None, :]) / (1 - p[:, None] - wc[None, :])                 # [k,b]
        blk -= np.conj(Vv)[None, :] * (q[r0:r1] @ P2)
        blk -= qE[r0:r1] @ np.conj(q).T
        H[r0:r1] = np.conj(blk)
    P3 = lam ** (1 - p[None, :] - wc[:, None]) / (1 - p[None, :] - wc[:, None])     # [a,k]
    c = That - (np.conj(Vv)[:, None] * P3) @ tau - np.conj(q) @ (E3 @ tau)
    x = np.linalg.solve(H, c)
    B = float(np.real(np.vdot(c, x)))
    if ret_parts:
        return B, float(np.linalg.eigvalsh(H).min()) if n <= 1500 else float('nan')
    return B

# ---------------- zeta validation (target chi_1, poles -2j) ----------------
zc = np.load(OUT + 'zeros_cache.npz'); gz = zc['g']; thz = zc['th']; th1z = zc['th1']
Kz = 200
gam = np.concatenate([gz[:Kz], -gz[:Kz]]); rho = 0.5 + 1j * gam
Vz = (rho / (1 - rho)) ** 3 * np.exp(2j * np.concatenate([thz[:Kz], -thz[:Kz]]))
dz = np.concatenate([2 * th1z[:Kz]] * 2) + 3 / np.abs(rho) ** 2
jj = np.arange(1, 25)
pz = -2.0 * jj
Rz = np.array([-2 * (-1) ** j * (2 * math.pi) ** (2 * j) * (2 * j / (2 * j + 1)) ** 3 / math.factorial(2 * j) for j in jj])
Tz = (rho - 1) / rho ** 2
tz = Rz * (pz - 1) / pz ** 2
for N in ([10, 100, 1280] if NMAX <= 400 else []):
    print(f'zeta check K={Kz}: N={N} B_generic={bound(N, rho, Vz, dz, np.ones(len(rho), bool), pz, Rz, Tz, tz):.10e}  (s3 K=200: N=10 1.1051872396e-02, N=100 6.6647406667e-03, N=1280 4.6981726073e-03)', flush=True)

# ---------------- DH ----------------
kap = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
S = [0.0, 1.0, 1.0 + kap, 1.0, 0.0]
def gdh_coprime(h, k):
    P = 5 * h * k
    br = np.union1d(np.arange(0, P + 1, h), np.arange(0, P + 1, k)).astype(np.float64)
    a, b = br[:-1], br[1:]
    c = np.array(S)[(a // h).astype(np.int64) % 5] * np.array(S)[(a // k).astype(np.int64) % 5]
    m = c != 0
    return float(np.sum(c[m] / P * (digamma(b[m] / P) - digamma(a[m] / P))))
G = np.zeros((NMAX, NMAX)); cache = {}
for h in range(1, NMAX + 1):
    for k in range(h, NMAX + 1):
        g_ = gcd(h, k); a_, b_ = h // g_, k // g_
        if (a_, b_) not in cache:
            cache[(a_, b_)] = gdh_coprime(a_, b_)
        G[h - 1, k - 1] = G[k - 1, h - 1] = cache[(a_, b_)] / g_
del cache
beta = sum(S[r] * (digamma((r + 1) / 5) - digamma(r / 5)) / 5 for r in range(1, 5))
bvec = np.array([beta / k for k in range(1, NMAX + 1)])
Lc = np.linalg.cholesky(G); y = np.linalg.solve(Lc, bvec); d2 = 1 - np.cumsum(y * y)
del G, Lc
print(f'DH Gram N<={NMAX} built {time.time()-t0:.1f}s; d2(10)={d2[9]:.8f} d2(400)={d2[min(399,NMAX-1)]:.8f} (explore13: d2 logN 0.122 at 10, 0.170 at 400)', flush=True)
on = np.array(json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/ExplorePotential/logs/dh_ft_2500.json'))['online'])
off = np.array(json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/ExplorePotential/logs/dh_off_2500.json'))['offline'])
print(f'DH zeros: {len(on)} on-line < 2500, {len(off)} off-line right members (beta>1/2) < 2500')
g1 = on[:KON]
def Vf(w):
    w = mp.mpc(w)
    return complex((w / (1 - w)) * mp.power(5 / mp.pi, w - 0.5) * mp.gamma((1 + w) / 2) / mp.gamma((2 - w) / 2))
rho_on = np.concatenate([0.5 + 1j * g1, 0.5 - 1j * g1])
V_on = np.array([Vf(r) for r in rho_on])
d_on = np.array([float(1 / (0.25 + t * t) + math.log(5 / math.pi) + mp.re(mp.digamma(mp.mpc(0.75, t / 2)))) for t in np.concatenate([g1, -g1])])
r0f = np.concatenate([off[:, 0] + 1j * off[:, 1], off[:, 0] - 1j * off[:, 1]])     # right members and conjugates
w_off = 1 - np.conj(r0f)
V_off = np.array([Vf(x) for x in w_off])
kk = np.arange(0, 20)
pd = -(2.0 * kk + 1)
Rd = np.array([float((p_ / (1 - p_)) * mp.power(5 / mp.pi, p_ - 0.5) * 2 * (-1) ** int(k_) / (mp.factorial(int(k_)) * mp.gamma((2 - p_) / 2))) for k_, p_ in zip(kk, pd)])
td = Rd / pd
print(f'V_f checks: |V_f(rho_1)| = {abs(V_on[0]):.12f}; residues R_0..2 = {Rd[:3]}', flush=True)
# diag tail over on-line zeros > KON (to 2500), then density log(5t/2pi)/2pi
gt = on[KON:]
def tail(L):
    s = 2 * float(np.sum((1 / (0.25 + gt ** 2)) / (L + np.log(5 * gt / (2 * math.pi)))))
    u = np.linspace(0, 60, 60001); t = on[-1] * np.exp(u); lg = np.log(5 * t / (2 * math.pi))
    return s + float(np.trapezoid(2 * (lg / (2 * math.pi)) / (t * t * (L + lg)) * t, u))
# C_f (all zeros) from explore13: 0.19527; on-line share from zeros we have:
Con = float(np.sum(1 / np.abs(rho_on) ** 2)) + 2 * float(np.sum(1 / (0.25 + gt ** 2)))
rows = []
NS = [n for n in [5, 10, 20, 40, 60, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000] if n <= NMAX]
NS = sorted(set(NS + list(range(5, NMAX + 1, 1 if NMAX <= 400 else 2))))
for N in NS:
    L = math.log(N)
    pts = np.concatenate([rho_on, w_off]); Vv = np.concatenate([V_on, V_off])
    mask = np.concatenate([np.ones(len(rho_on), bool), np.zeros(len(w_off), bool)])
    dph = np.concatenate([d_on, np.zeros(len(w_off))])
    That = 1 / (1 - np.conj(pts))
    Bon = bound(N, rho_on, V_on, d_on, np.ones(len(rho_on), bool), pd, Rd, 1 / (1 - np.conj(rho_on)), td)
    Ball, emin = bound(N, pts, Vv, dph, mask, pd, Rd, That, td, ret_parts=True) if N in (5, 100, NMAX) else (bound(N, pts, Vv, dph, mask, pd, Rd, That, td), float('nan'))
    tl = tail(L)
    rows.append((N, L, d2[N - 1], Bon + tl, Ball + tl))
    if N in (5, 10, 20, 40, 100, 200, 400, 600, 800, 1000) or N == NMAX:
        print(f'DH N={N:5d} d2={d2[N-1]:.6e} d2L={d2[N-1]*L:.5f}  B_on+tail={Bon+tl:.6e} B_all+tail={Ball+tl:.6e} ratio={(Ball+tl)/d2[N-1]:.4f}  L^2(d2-B_all)={L*L*(d2[N-1]-Ball-tl):+.5f}  eigmin={emin:.3f}', flush=True)
np.save(OUT + 'dh_compare.npy', np.array(rows))
print(f'DH on-line C share (zeros < 2500 + tail-free) = {Con:.6f}; done {time.time()-t0:.1f}s maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB')
