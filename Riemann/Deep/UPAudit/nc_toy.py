# N-C: Lemma Z against planted off-line zeros in the explicit-formula toy model (no primes needed, so L can be large).
# Toy "prime" measure nu_s(u) = -phi(u) sum_{rho in Z} e^{(rho-s)u} (Z1 without pole/trivial terms);
# M_toy(s) = sum_{j<=k} (-1)^j/j! int f(v/L) nu_s^{*j}(v) dv  (continuous analogue of Z2: no collisions).
# Background: real zeta zeros near 15000 rescaled to density L/2pi; planted pairs rho = 1/2 +- d + i gamma_p.
# Ratios: R+ = |M(1/2+y+it)|/((1+W_{-y})(1+On)^{k-1}), R- = |M(1/2-y+it)|/(X^y (1+W)(1+On)^{k-1}), A = 12.
# Naive route (DPC l.665-670) would allow |M| ~ X^{k d}; Lemma Z says X^{d} up to (1+On)^{k-1} factors.
import numpy as np, sys, json, math
ZEROS = np.load('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/zeta_zeros.npy')
def S(x):
    x = np.clip(np.asarray(x, float), 0.0, 1.0)
    g = lambda y: np.where(y > 0, np.exp(-1/np.maximum(y, 1e-300)), 0.0)
    return g(x)/(g(x)+g(1-x))
class Toy:
    def __init__(m, L, th=0.56, eta=0.1, b=0.66, du=0.01, A=12):
        m.L, m.th, m.eta, m.b, m.du, m.A = L, th, eta, b, du, A; m.X = math.exp(th*L); m.k = int(th/eta)
        m.u = np.arange(0, b*L+du, du); v = m.u/L
        rho = (1-S((v-th)/(b-th)))*S(v/eta-1); m.ph = rho/np.maximum(m.u, 1e-300)
    def f(m, x):
        ins = np.abs(x) < m.th; z = np.where(ins, 1-(x/m.th)**2, 1.0)
        return np.where(ins, np.exp(1-1/np.maximum(z, 1e-300)), 0.0)*(1+0.5j*x)
    def M(m, zs, s):
        E = np.zeros(len(m.u), complex)
        for i in range(0, len(zs), 100): E += np.exp(np.outer(m.u, zs[i:i+100]-s)).sum(axis=1)
        g = -m.ph*E; tot = complex(m.f(np.array(0.0))); pw = g.copy()
        for j in range(1, m.k+1):
            v = np.arange(len(pw))*m.du
            tot += (-1)**j/math.factorial(j)*np.sum(m.f(v/m.L)*pw)*m.du
            if j < m.k: pw = np.convolve(pw, g)*m.du
        return abs(tot)
    def kap(m, x): return np.minimum(1.0, np.maximum(x, 1e-300)**(-m.A))
    def bounds(m, zs, t, y):
        kk = m.kap(m.L*np.abs(zs.imag-t)); d = zs.real-0.5
        On = kk.sum(); W = (m.X**np.maximum(d, 0)*kk).sum(); Wy = (m.X**np.maximum(d-y, 0)*kk).sum()
        return (1+Wy)*(1+On)**(m.k-1), m.X**y*(1+W)*(1+On)**(m.k-1), W, On
def background(L, t0, Gam):
    z = ZEROS[(ZEROS > 15000-400) & (ZEROS < 15000+400)]; dens = np.log(15000/(2*np.pi))/(2*np.pi)
    g = (z-15000)*dens/(L/(2*np.pi)) + t0; g = g[np.abs(g-t0) <= Gam]
    return 0.5+1j*g
def planted(cfg, L, t0, a):
    d = a/L; h = 2*np.pi*0.64/L
    P = {'single': [(0.0, d)], 'single_off': [(0.7/L, d)],
         'cluster3': [(0.0, d), (0.3/L, 0.8*d), (-0.4/L, 0.6*d)],
         'coincident3': [(0.0, d), (0.0, d), (0.0, d)],
         'cluster5': [(0.0, d), (0.2/L, d), (-0.2/L, 0.9*d), (0.5/L, 0.7*d), (-0.6/L, d)],
         'hcoincident': [(0.0, d), (h, d)]}[cfg]
    out = []
    for (dg, dd) in P: out += [0.5+dd+1j*(t0+dg), 0.5-dd+1j*(t0+dg)]
    return np.array(out)
def bounds_A(m, zs, t, y, A):
    kk = np.minimum(1.0, np.maximum(m.L*np.abs(zs.imag-t), 1e-300)**(-float(A))); d = zs.real-0.5
    On = kk.sum(); W = (m.X**np.maximum(d, 0)*kk).sum(); Wy = (m.X**np.maximum(d-y, 0)*kk).sum()
    return (1+Wy)*(1+On)**(m.k-1), m.X**y*(1+W)*(1+On)**(m.k-1)
if __name__ == '__main__':
    tiny = len(sys.argv) > 1
    Ls = [20.0] if tiny else [20.0, 30.0, 44.0]
    cfgs = ['single', 'cluster3'] if tiny else ['single', 'single_off', 'cluster3', 'coincident3', 'cluster5', 'hcoincident']
    res = []
    for L in Ls:
        m = Toy(L, du=0.01 if L < 40 else 0.012); t0 = 1000.0; Gam = 30.0
        bg = background(L, t0, Gam)
        for a in ([2, 8] if tiny else [2, 5, 8, 12, round(0.48*L, 2)]):
            d = a/L
            for cfg in cfgs:
                zs = np.concatenate([bg, planted(cfg, L, t0, a)])
                worst = {f'R{sg}_A{A}': 0.0 for sg in '+-' for A in (2, 4, 12)}; Mc = {}
                for dt in np.linspace(-3/L, 3/L, 13):
                    t = t0+dt
                    for y in (0.0, 0.5*d, d, 0.5):
                        Mp = m.M(zs, 0.5+y+1j*t); Mm = m.M(zs, 0.5-y+1j*t)
                        for A in (2, 4, 12):
                            Bp, Bm = bounds_A(m, zs, t, y, A)
                            worst[f'R+_A{A}'] = max(worst[f'R+_A{A}'], Mp/Bp); worst[f'R-_A{A}'] = max(worst[f'R-_A{A}'], Mm/Bm)
                        if abs(dt) < 1e-12: Mc[(round(y, 6), '+')] = Mp; Mc[(round(y, 6), '-')] = Mm
                M0 = Mc[(0.0, '+')]
                ex = {f'e+({c})': math.log(Mc[(round(c*d, 6), '+')]/M0)/(m.th*L) for c in (0.5, 1.0)}
                ex.update({f'e-({c})': math.log(Mc[(round(c*d, 6), '-')]/M0)/(m.th*L) for c in (0.5, 1.0)})
                row = dict(L=L, cfg=cfg, a=a, d=round(d, 4), k=m.k, logXd=round(m.th*a, 3), logM0=math.log(M0), **worst, **ex)
                res.append(row); print(json.dumps({kk: (round(v, 4) if isinstance(v, float) else v) for kk, v in row.items()}), flush=True)
    json.dump(res, open('/Users/adriansager/playground/riemann/Riemann/Deep/UPAudit/logs/nc_toy%s.json' % ('_tiny' if tiny else ''), 'w'), indent=1)
