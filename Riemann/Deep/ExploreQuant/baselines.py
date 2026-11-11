# explore29: GUE (sine-kernel) and CUE / Poisson baselines for the quantizer statistics; rigorous TV(log p) lower bound.
# usage: uv run --no-project --with numpy --with scipy python3 baselines.py <dl dir>
import sys, gzip, math, time
import numpy as np

DL = sys.argv[1]
rng = np.random.default_rng(20261111)
t0 = time.time()

# ---- exact sine-kernel gap probability E(0;s) by Bornemann's Nystrom method
xg, wg = np.polynomial.legendre.leggauss(48)
def E0(s):
    if s == 0: return 1.0
    x = (xg + 1)*s/2; w = wg*s/2
    D = x[:, None] - x[None, :]
    with np.errstate(invalid='ignore', divide='ignore'):
        K = np.where(np.abs(D) < 1e-15, 1.0, np.sin(math.pi*D)/(math.pi*D))
    A = np.sqrt(w)[:, None]*K*np.sqrt(w)[None, :]
    return float(np.linalg.det(np.eye(len(x)) - A))
sg, swg = np.polynomial.legendre.leggauss(160)
S = (sg + 1)*4.0; SW = swg*4.0          # [0, 8]
Ev = np.array([E0(s) for s in S])
m2 = 2*np.sum(SW*Ev); m3 = 6*np.sum(SW*S*Ev); m4 = 12*np.sum(SW*S**2*Ev)
print(f"sine kernel (GUE, exact): E s^2 = {m2:.6f} (var {m2-1:.6f})  E s^3 = {m3:.6f} ({10*math.log10(m3):.3f} dB)  E s^4 = {m4:.6f}  [E0(8)={E0(8.0):.2e}] [{time.time()-t0:.1f}s]")
# E log s: p(s) = E0''(s); E[log s] = int log s p ds; use p from finite differences on a grid
h = 2e-3; grid = np.arange(h, 6.0, h)
Eg = np.array([E0(s) for s in np.concatenate(([0.0], grid, [6.0+h]))])
p = (Eg[2:] - 2*Eg[1:-1] + Eg[:-2])/h**2
print(f"  check int p = {np.sum(p)*h:.6f}, int s p = {np.sum(grid*p)*h:.6f}, E log s = {np.sum(np.log(grid)*p)*h:+.6f} [{time.time()-t0:.1f}s]")

# ---- the rigorous oscillation bound: osc_V log p >= eps(o), o = |g_n - g_{n-1}|/(g_n + g_{n-1})
# max centroid offset on [-1,1] for p in [1, e^eps] is max_u K u(2-u)/(2(2+K u)), K = e^eps - 1
uu = np.linspace(1e-6, 2, 20001)
def maxoff(K):
    return np.max(K*uu*(2-uu)/(2*(2 + K*uu)))
Ks = np.exp(np.linspace(-12, 12, 2401))
offs = np.array([maxoff(K) for K in Ks])   # increasing in K
def eps_of(o):
    o = np.asarray(o)
    return np.log1p(np.interp(o, offs, Ks, left=0.0, right=np.inf))
print(f"  oscillation table: o=0.1 -> eps {float(eps_of(0.1)):.4f}, 0.3 -> {float(eps_of(0.3)):.4f}, 0.5 -> {float(eps_of(0.5)):.4f}, 0.8 -> {float(eps_of(0.8)):.4f}")

def report(name, d):
    d = d/np.mean(d)
    pk = np.mean((d[1:-1] > d[:-2]) & (d[1:-1] > d[2:]))
    c1 = np.corrcoef(d[:-1], d[1:])[0, 1]
    o = np.abs(d[1:] - d[:-1])/(d[1:] + d[:-1])
    ev = eps_of(o)
    print(f"{name}: n={len(d)}  E d^2={np.mean(d**2):.5f}  E d^3={np.mean(d**3):.5f}  E d^4={np.mean(d**4):.5f}  corr1={c1:+.4f}  peakfrac={pk:.4f}  "
          f"E log d={np.mean(np.log(d)):+.5f}  osc bound: mean eps={np.mean(ev):.4f}, median {np.median(ev):.4f} per cell")

# ---- CUE
sp = []
Nm = 200
for k in range(1000):
    Z = (rng.standard_normal((Nm, Nm)) + 1j*rng.standard_normal((Nm, Nm)))/math.sqrt(2)
    Q, R = np.linalg.qr(Z)
    Q = Q*(np.diag(R)/np.abs(np.diag(R)))
    th = np.sort(np.angle(np.linalg.eigvals(Q)))
    dd = np.diff(np.concatenate((th, [th[0] + 2*math.pi])))*Nm/(2*math.pi)
    sp.append(dd)
report(f"CUE N={Nm} x 1000 (per-matrix sequences concatenated)", np.concatenate(sp))
print(f"  [{time.time()-t0:.1f}s]")
report("Poisson (iid exponential)", rng.exponential(size=200000))
report("crystal + 0.1 jitter", 1 + 0.1*(rng.standard_normal(200001)[1:] - rng.standard_normal(200001)[:-1]))

with gzip.open(DL + '/zeros1.gz', 'rt') as f:
    g = np.array([float(x) for x in f.read().split()])
th = g/2*np.log(g/(2*math.pi)) - g/2 - math.pi/8 + 1/(48*g)
d = np.diff(th/math.pi)
report("zeta zeros 10001..100000", d[10000:])
for nm, base in [('zeros3', 267653395647.0), ('zeros5', 1.3709199099319953e21)]:
    v = []
    for line in open(DL + '/' + nm):
        try: v.append(float(line.strip()))
        except ValueError: pass
    report(f"zeta {nm}", np.diff(np.array(v)))
print(f"done [{time.time()-t0:.1f}s]")
