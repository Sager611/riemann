# t3_zeros.py — zeros off the critical line of the L-series attached to torus knots (Hikami's strange identities)
# A: trefoil chi12 (Dirichlet L, control); B: T(2,5) knot series Im(psi20); C: scalar-FE eigencombination G (DH-type) of the T(2,5) vector;
# D: T(3,4) series (1/2)[(24/.) - (12/.)]; E: T(2,7) series mod 28.
import numpy as np, mpmath as mp, sys, json, time
T0 = time.time()
TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 200.0
mp.mp.dps = 30

def tk(s, t):  # torus knot periodic function mod 2st (Hikami): +1 at +-(st-s-t), -1 at +-(st+s-t)
    f = 2*s*t; a = np.zeros(f)
    for r, v in [(s*t-s-t, 1), (s*t+s-t, -1), (s*t-s+t, -1), (s*t+s+t, 1)]:
        a[r % f] += v
    return f, a

# psi20: even primitive char mod 20 = psi4 * psi5, psi5(2) = i
def psi20():
    a = np.zeros(20, complex)
    for n in range(20):
        if np.gcd(n, 20) != 1: continue
        p4 = 1 if n % 4 == 1 else -1
        k = [j for j in range(4) if pow(2, j, 5) == n % 5][0]
        a[n] = p4 * (1j)**k
    return a
P = psi20()
tau = sum(P[n]*np.exp(2j*np.pi*n/20) for n in range(20))
eps = tau/np.sqrt(20)
print('psi20(3) =', P[3], ' root number eps =', eps, ' |eps| =', abs(eps), ' tan(arg eps) =', np.tan(np.angle(eps)), ' 1/golden =', (np.sqrt(5)-1)/2)
f25, a25 = tk(2, 5)
print('T(2,5) periodic fn equals Im psi20:', np.allclose(a25, P.imag))
G = 2*np.real(eps**-0.5 * P)          # scalar FE, sign +1
H = 2*np.imag(eps**-0.5 * P)          # scalar FE, sign -1 (H/(2i) up to normalisation)
fam = {}
f12, a12 = tk(2, 3); fam['A_trefoil_chi12'] = (12, a12)
fam['B_T25_Im_psi20'] = (20, a25)
fam['C_T25_G_plus'] = (20, G)
fam['C2_T25_H_minus'] = (20, H)
f34, a34 = tk(3, 4); fam['D_T34'] = (24, a34)
f27, a27 = tk(2, 7); fam['E_T27'] = (28, a27)

# Bernoulli numbers for Euler–Maclaurin
B2 = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330]
from math import factorial
NEM = 150
def hurwitz_np(s, x):
    n = np.arange(NEM)[None, :] + x
    S = np.exp(-s[:, None]*np.log(n)).sum(1)
    X = NEM + x; lX = np.log(X)
    S += np.exp((1-s)*lX)/(s-1) + 0.5*np.exp(-s*lX)
    poch = s.copy()
    for j in range(1, len(B2)+1):
        S += B2[j-1]/factorial(2*j)*poch*np.exp((-s-2*j+1)*lX)
        poch = poch*(s+2*j-1)*(s+2*j)
    return S
def L_np(s, f, a):
    out = np.zeros_like(s)
    for r in range(1, f+1):
        c = a[r % f]
        if c != 0: out += c*hurwitz_np(s, r/f)
    return out*np.exp(-s*np.log(f))
def L_mp(s, f, a):
    return mp.power(f, -s)*mp.fsum(complex(a[r % f])*mp.zeta(s, mp.mpf(r)/f) for r in range(1, f+1) if a[r % f] != 0)

# sanity: compare numpy vs mpmath at a test point
s0 = 0.7+123.4j
print('EM check:', abs(L_np(np.array([s0]), 20, G)[0] - complex(L_mp(mp.mpc(s0), 20, G))))

res = {}
sig = np.arange(0.53, 1.80, 0.02)
tt = np.arange(0.5, TMAX, 0.04)
for name, (f, a) in fam.items():
    Z = np.empty((len(sig), len(tt)))
    for i, sg in enumerate(sig):
        for j0 in range(0, len(tt), 4000):
            s = sg + 1j*tt[j0:j0+4000]
            Z[i, j0:j0+4000] = np.abs(L_np(s, f, a))
    # local minima in 2D
    cand = []
    for i in range(1, len(sig)-1):
        row = Z[i]
        m = (row[1:-1] < row[:-2]) & (row[1:-1] < row[2:]) & (row[1:-1] < Z[i-1, 1:-1]) & (row[1:-1] < Z[i+1, 1:-1])
        for j in np.nonzero(m)[0] + 1:
            if Z[i, j] < 0.15: cand.append((sig[i], tt[j]))
    zeros = []
    for (sg, t) in cand:
        try:
            z = mp.findroot(lambda s: L_mp(s, f, a), mp.mpc(sg, t), tol=1e-20, maxsteps=50)
        except Exception:
            continue
        if abs(L_mp(z, f, a)) < 1e-15 and z.real > 0.5 + 1e-8 and 0 < z.imag < TMAX:
            zz = (round(float(z.real), 10), round(float(z.imag), 8))
            if all(abs(zz[0]-w[0]) + abs(zz[1]-w[1]) > 1e-6 for w in zeros):
                zeros.append(zz)
    zeros.sort(key=lambda w: w[1])
    res[name] = zeros
    print(f'{name}: {len(zeros)} zeros with Re s > 1/2, 0 < t < {TMAX:.0f}; first:', zeros[:6], f' [{time.time()-T0:.0f}s]', flush=True)
json.dump(res, open('t3_zeros_out.json', 'w'), indent=1)
