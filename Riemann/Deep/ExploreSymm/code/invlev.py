# Quantized inverse Levinson problem: P-periodic f(z) = p(w), w = exp(-z/Pz) (Pz = zeta zeros per period),
# p(w) = prod_j (1 - w/w_j), c_0 = 1. Constraint on every measured line x_l: mean|f|^2 = S(r_l) <= C_l, r_l = exp(-x_l/Pz).
# Maximise M = #roots with |w_j| <= 1 (zeros with Re z >= 0). LP (smeared) bound: M <= Pz * min_l h_l/|x_l|.
# usage: python3 invlev.py profile.npz Pz dens M nstart seed
import numpy as np, sys, math, time
from scipy.optimize import minimize
prof = np.load(sys.argv[1]); xs = prof['xs']; hl = prof['h']
Pz = float(sys.argv[2]); dens = float(sys.argv[3]); M = int(sys.argv[4]); nstart = int(sys.argv[5]); seed = int(sys.argv[6])
K = int(round(dens*Pz)); Nf = 2*K + 2
r = np.exp(-xs/Pz); logC = 2*hl
ang = 2*np.pi*np.arange(Nf)/Nf
W = r[:, None]*np.exp(1j*ang)[None, :]          # nl x Nf
lp_bound = Pz*min(hl[i]/abs(xs[i]) for i in range(len(xs)) if xs[i] < 0)
tau = 0.01
sgn = np.where(np.arange(K) < M, -1.0, 1.0)     # u_j = sgn_j * a_j^2
def unpack(p):
    a = p[:K]; v = p[K:]; u = sgn*a*a; return a, u, v
def F(p, grad=True):
    a, u, v = unpack(p); wj = np.exp(u + 1j*v)
    nl = W.shape[0]; logS = np.zeros(nl); oms = []
    CH = 8
    for c0 in range(0, nl, CH):
        Wc = W[c0:c0+CH]
        D = wj[None, None, :] - Wc[:, :, None]
        lp2 = np.sum(np.log(np.abs(D)**2), axis=2) - 2*np.sum(u)
        mx = lp2.max(axis=1, keepdims=True); e = np.exp(lp2 - mx)
        logS[c0:c0+CH] = np.log(e.mean(axis=1)) + mx[:, 0]
        oms.append(e/e.sum(axis=1, keepdims=True)); del D
    g = (logS - logC)/tau; gm = g.max(); pi = np.exp(g - gm); Z = pi.sum(); pi /= Z
    val = tau*(np.log(Z) + gm)
    if not grad: return val, logS
    du = np.zeros(K); dv = np.zeros(K)
    for ci, c0 in enumerate(range(0, nl, CH)):
        Wc = W[c0:c0+CH]; G = wj[None, None, :]/(wj[None, None, :] - Wc[:, :, None])
        wts = pi[c0:c0+CH, None]*oms[ci]
        du += np.einsum('lk,lkj->j', wts, 2*G.real); dv += np.einsum('lk,lkj->j', wts, -2*G.imag)
    du -= 2
    da = du*sgn*2*a
    return val, np.concatenate([da, dv])
rng = np.random.default_rng(seed)
def init(kind):
    a = np.zeros(K); v = np.zeros(K)
    # inside roots: just inside the unit circle, equally spaced
    a[:M] = math.sqrt(0.2/Pz); v[:M] = 2*np.pi*(np.arange(M) + rng.uniform(0, 1, M)*0.3)/max(M, 1)
    nout = K - M
    if kind == 0:   # outside x-positions spread like the measured zero density of V*Psi (roughly uniform in [-6, -0.2])
        xo = -np.sort(rng.uniform(0.2, 6.0, nout))
    elif kind == 1:
        xo = -rng.exponential(2.0, nout) - 0.05
    else:
        xo = -np.linspace(0.3, 8.0, nout)
    a[M:] = np.sqrt(-xo/Pz)
    v[M:] = rng.uniform(0, 2*np.pi, nout) if kind != 2 else 2*np.pi*np.arange(nout)*0.6180339887 % (2*np.pi)
    return np.concatenate([a, v])
best = None; t0 = time.time()
for st in range(nstart):
    p0 = init(st % 3)
    res = minimize(F, p0, jac=True, method='L-BFGS-B', options={'maxiter': 4000, 'maxfun': 8000})
    val, logS = F(res.x, grad=False)
    viol = np.max(logS - logC)
    if best is None or viol < best[0]: best = (viol, res.x.copy(), logS)
    print(f"start {st} kind {st%3}: max violation (log S - log C) = {viol:+.5f}  smooth {val:+.5f}  it {res.nit}  {time.time()-t0:.0f}s", flush=True)
viol, x, logS = best; a, u, v = unpack(x)
nin = int(np.sum(u <= 1e-12))
print(f"Pz={Pz} K={K} M={M} LP bound M<={lp_bound:.3f} (nu<={lp_bound/Pz:.4f}); best max violation {viol:+.5f}; roots in closed disc {nin}")
print("line x, h_target, h_competitor:", ' '.join(f"{xs[i]:.2f}:{hl[i]:.4f}/{logS[i]/2:.4f}" for i in range(len(xs))))
xr = -u*Pz   # zero real parts in x units (x = -Pz*log|w|)
print("inside roots x:", np.round(np.sort(xr[:M]), 3))
print("outside roots x quantiles:", np.round(np.quantile(xr[M:], [0, .1, .25, .5, .75, .9, 1]), 3) if K > M else '')
np.savez(f"inv_Pz{int(Pz)}_M{M}_s{seed}.npz", x=x, u=u, v=v, logS=logS, xs=xs, hl=hl, viol=viol, Pz=Pz, K=K, M=M)
