# explore55: empirical 1D (CGG) vs 2D bounds with the asymptotically optimal polynomials, vs the recipe predictions.
import numpy as np, sys, json, re
from moll import mobius, bound_complex, bound_real
from rs import zprime, theta as rstheta
def asym(th, a):                      # P = a x^2 + (1-a) x
    I = 0.5 - a/6; Jp = 1 + a*a/3; u = 0.5 + th*I; A = u*u + 1/12 + Jp/(12*th); B = u*u - 1/12
    return dict(u=u, A=A, B=B, k1D=u*u/A, k2D=2*u*u/(A+B))
def best(th, key):
    aa = np.linspace(-2, 1, 30001); v = np.array([asym(th, a)[key] for a in aa]); i = v.argmax(); return aa[i], v[i]
def mol(g, y, polys, chunk=1000):
    mu = mobius(int(y)); n = np.nonzero(mu[1:int(y)+1])[0] + 1; x = np.log(y/n)/np.log(y); w = mu[n]*n**-0.5
    Wt = np.stack([w*p(x) for p in polys], 1); ln = np.log(n.astype(float)); M = np.empty((len(g), len(polys)), complex)
    for a in range(0, len(g), chunk):
        M[a:a+chunk] = np.exp(-1j*np.outer(g[a:a+chunk], ln)) @ Wt
    return M
def run(g, D, tag):
    T = g.max(); L = np.log(T/(2*np.pi)); out = []
    for th in [0.25, 0.5]:
        a1, k1 = best(th, 'k1D'); a2, k2 = best(th, 'k2D')
        polys = [lambda x, a=a1: a*x*x+(1-a)*x, lambda x, a=a2: a*x*x+(1-a)*x, lambda x: x, lambda x: x*x, lambda x: x**3]
        Y = D[:, None]*mol(g, T**th, polys)
        z = Y[:, 2]; ph = z.mean()/abs(z.mean()); z = z/ph
        r = dict(tag=tag, L=round(L, 3), N=len(g), theta=th, a_1D=round(a1, 4), a_2D=round(a2, 4),
                 pred_1D=round(k1, 4), pred_2D=round(k2, 4), CGG_P1D=bound_complex(Y[:, [0]]), real2D_P2D=bound_real(Y[:, [1]]),
                 CGG_J3=bound_complex(Y[:, 2:5]), real2D_J3=bound_real(Y[:, 2:5]),
                 u_x=abs(z.mean())/L, A_x=np.mean(abs(z)**2)/L**2, B_x=np.mean(z**2).real/L**2, ImB_x=np.mean(z**2).imag/L**2,
                 pred_x=asym(th, 0.0))
        out.append(r); print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
    return out
if __name__ == '__main__':
    src = sys.argv[1]
    if src.endswith('zeros3.txt'):
        import mpmath as mp; mp.mp.dps = 40; t0 = 267653395647
        d = np.array([float(l) for l in open(src) if re.fullmatch(r'\s*[0-9]+\.[0-9]+\s*', l)])
        zp, _ = np.load(src.replace('.txt', '_zp10000.npy'))
        th = np.array([float(mp.fmod(mp.siegeltheta(mp.mpf(t0)+mp.mpf(str(x))), 2*mp.pi)) for x in d])
        print('high-height mode needs moll_hi phases; skipped', file=sys.stderr); sys.exit(0)
    g = np.loadtxt(src); a, b = int(sys.argv[2]), int(sys.argv[3]); g = g[a:b]
    D = -1j*np.exp(-1j*rstheta(g))*zprime(g); run(g, D, src.split('/')[-1]+f'[{a}:{b}]')
