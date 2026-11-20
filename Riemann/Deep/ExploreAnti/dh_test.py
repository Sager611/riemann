# explore55 adversarial: the same 1D/2D statistic for the Davenport-Heilbronn function (on-line zeros, f'(rho) = -i e^{-i vartheta} Z_DH'),
# "mollifier" = truncated Dirichlet inverse of DH's coefficients (no Euler product); compared with zeta at the same height.
import sys, time, resource, numpy as np
RSS = lambda tag: print('#', tag, round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20, 1), 'MB', flush=True)
sys.path.insert(0, '../ExploreConcentration/code')
import dhlib
from moll import bound_complex, bound_real, mobius
from rs import zprime, theta as zth
TMAX = float(sys.argv[1]); h = 1e-5
def Zc(t):
    t = np.atleast_1d(t); return np.concatenate([dhlib.Zdh(t[i:i+512]) for i in range(0, len(t), 512)])
ts = np.arange(5.0, TMAX, 0.02); z = Zc(ts)
RSS('grid'); idx = np.nonzero(np.sign(z[:-1])*np.sign(z[1:]) < 0)[0]; lo, hi, flo = ts[idx], ts[idx+1], z[idx]
for it in range(40):
    mid = (lo+hi)/2; fm = Zc(mid); sm = np.sign(fm) == np.sign(flo)
    lo = np.where(sm, mid, lo); flo = np.where(sm, fm, flo); hi = np.where(sm, hi, mid)
RSS('bisect'); g = (lo+hi)/2; zpd = (Zc(g+h) - Zc(g-h))/(2*h)
Dd = -1j*np.exp(-1j*dhlib.vartheta(g))*zpd
a = np.array([0.0, 1.0, dhlib.k, -dhlib.k, -1.0])        # a_n by n mod 5
def inv(y):
    b = np.zeros(y+1); b[1] = 1
    for n in range(2, y+1):
        s = 0.0
        for d in range(1, n):
            if n % d == 0: s += b[d]*a[(n//d) % 5]
        b[n] = -s
    return b
def stats(gg, D, coef, y, tag):
    n = np.arange(1, int(y)+1); x = np.log(y/n)/np.log(y); w = coef[n]*n**-0.5*x
    M = np.exp(-1j*np.outer(gg, np.log(n))) @ w; Y = D*M; L = np.log(gg.max()/(2*np.pi))
    zz = Y/(Y.mean()/abs(Y.mean())); u = abs(Y.mean())/L; A = np.mean(abs(zz)**2)/L**2; B = np.mean(zz**2)/L**2
    print(f"{tag}: N={len(gg)} L={L:.2f} y={y:.1f} CGG={bound_complex(Y[:, None]):.4f} 2D={bound_real(Y[:, None]):.4f} u={u:.4f} A={A:.4f} B={B.real:.4f}{B.imag:+.4f}i B-u^2={B.real-u*u:.4f} A-u^2={A-u*u:.4f}", flush=True)
RSS('deriv'); print('DH on-line zeros', len(g), 'in [5,', TMAX, '] ; |Z_DH(g)| max', float(np.max(abs(Zc(g)))), flush=True)
gz = np.loadtxt('../../../.scratch/explore55/zeros1.txt'); gz = gz[gz < TMAX]; Dz = -1j*np.exp(-1j*zth(gz))*zprime(gz)
for th in [0.25, 0.5]:
    yd = TMAX**th; stats(g, Dd, inv(int(yd)), yd, f'DH   theta={th}')
    mu = mobius(int(yd)).astype(float); stats(gz, Dz, mu, yd, f'zeta theta={th}')
for th in [0.25, 0.5]:   # DH with the zeta mollifier (Mobius): wrong inverse
    yd = TMAX**th; stats(g, Dd, mobius(int(yd)).astype(float), yd, f'DH(mu-mollifier) theta={th}')
