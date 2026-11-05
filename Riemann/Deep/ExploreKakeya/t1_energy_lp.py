
# t1_energy_lp.py -- what would a better additive-energy bound for large-value sets buy inside the Guth-Maynard chain?
# Exponents in units of log N: |W| = N^r, T0 = N^t0 (interval length), E(W) = N^e, V = N^s.
# GM chain (arXiv 2405.20552): |W|^3 N^{6s-3} <~ N^3 + S2 + S3  (Prop 4.6 + (5.5));
#   S2 <~ N^2|W|^2 + T N |W|^{2-1/k} + N^2 |W|^2 (T^{1/2}/|W|^{3/4})^{1/k}  (Prop 6.1, any k)
#   S3 <~ T^2 |W|^{3/2} + T N |W|^{1/2} E^{1/2}                           (Prop 10.1)
#   E  <~ energy hypothesis H (GM Prop 11.1 = H0)
# LV(s,tau): split [0,N^tau] into intervals of length N^t0 (t0 in [1,4/3]; Prop 11.1 needs N >= T^{3/4});
# also MVT and Huxley. Zero density: GM sec. 13 / Ingham zero detector: lengths N^k, N = T^nu, nu in (0,1/2].
import numpy as np, sys
R = np.linspace(0.0, 3.0, 30001)
T0GRID = np.linspace(1.0, 4/3, 201)

def e_upper(H, r, s, t0):
    if H.endswith('+C73'): H = H[:-4]
    if H == 'H1C': H = 'H1'
    if H == 'H0':   # GM Prop 11.1
        e = np.maximum.reduce([r + 4 - 4*s, 21*r/8 + t0/4 + 1 - 2*s, 3*r + 1 - 2*s])
    elif H == 'H1': # Heath-Brown third term removed (Conj 7.3): drop the |W|^{21/8} T^{1/4} N^{1-2s} term
        e = np.maximum(r + 4 - 4*s, 3*r + 1 - 2*s)
    elif H == 'H2': # 'union of short APs' structure theorem: E <~ |W|^2 N^{2-2s} + |W|^4/T + |W|^2
        e = np.maximum.reduce([2*r + 2 - 2*s, 4*r - t0, 2*r])
    elif H == 'H3': # random energy
        e = np.maximum(2*r, 4*r - t0)
    elif H.startswith('C'): # GM bottleneck term improved by N^{-c}
        c = float(H[1:])
        e = np.maximum.reduce([r + 4 - 4*s - c, 21*r/8 + t0/4 + 1 - 2*s, 3*r + 1 - 2*s])
    return np.minimum(e, 3*r)

def rho0(s, t0, H, kmax=12):
    lhs = 3*R + 6*s - 3
    eu = e_upper(H, R, s, t0)
    feas_E = eu >= 2*R - 1e-12          # E >= |W|^2 always
    best = 9.0
    for k in range(1, kmax+1):
        s2c = 2 + 2*R + (t0/2 - 0.75*R)/k        # Heath-Brown third term inside S2 (Prop 6.1)
        if H.endswith('+C73'): s2c = np.full_like(R, -9.0)   # Conj 7.3: no third term
        rhs = np.maximum.reduce([np.full_like(R, 3.0), 2 + 2*R, t0 + 1 + (2 - 1.0/k)*R,
                                 s2c, 2*t0 + 1.5*R, t0 + 1 + R/2 + eu/2])
        ok = (lhs <= rhs + 1e-12) & feas_E
        rmax = R[ok].max() if ok.any() else 0.0
        best = min(best, rmax)
    return best

_tab = {}
def table(s, H):
    key = (round(s, 6), H)
    if key not in _tab:
        _tab[key] = np.array([rho0(s, t0, H) for t0 in T0GRID])
    return _tab[key]

def LV(s, tau, H):
    if tau <= 1: return 2 - 2*s
    best = min(max(2 - 2*s, tau + 1 - 2*s), max(2 - 2*s, tau + 4 - 6*s))   # MVT, Huxley
    if H is None or H == 'none': return best
    tab = table(s, H)
    m = T0GRID <= tau + 1e-12
    best = min(best, float(np.min((tau - T0GRID[m]) + tab[m])))
    return best

def density(s, H, nus=None):
    if nus is None: nus = np.linspace(0.01, 0.5, 491)
    worst = 2*(1 - s)
    for nu in nus:
        vals = []
        k = 1
        while k*nu <= 3.0:
            x = k*nu
            vals.append(x*LV(s, 1/x, H))
            k += 1
        worst = max(worst, min(vals))
    return worst/(1 - s)

def mvt_only_density(s, nus=None):
    if nus is None: nus = np.linspace(0.01, 0.5, 491)
    worst = 2*(1 - s)
    for nu in nus:
        vals = []
        k = 1
        while k*nu <= 3.0:
            x = k*nu; tau = 1/x
            lv = 2 - 2*s if tau <= 1 else max(2 - 2*s, tau + 1 - 2*s)
            vals.append(x*lv); k += 1
        worst = max(worst, min(vals))
    return worst/(1 - s)

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'check':
        for s in [0.70, 0.75, 0.80]:
            print('s=%.2f rho0(H0,6/5)=%.4f GM(18-20s)/5=%.4f' % (s, rho0(s, 1.2, 'H0'), (18-20*s)/5))
        for s in [0.6, 0.7, 0.75, 0.8]:
            print('s=%.2f A(MVT only)=%.4f Ingham 3/(2-s)=%.4f | A(MVT+Huxley)=%.4f Huxley 3/(3s-1)=%.4f' % (
                s, mvt_only_density(s), 3/(2-s), density(s, 'none'), 3/(3*s-1)), flush=True)
    elif mode == 'dens':
        for H in sys.argv[2].split(','):
            for s in [float(x) for x in sys.argv[3].split(',')]:
                print('%s s=%.3f A=%.4f' % (H, s, density(s, H)), flush=True)
