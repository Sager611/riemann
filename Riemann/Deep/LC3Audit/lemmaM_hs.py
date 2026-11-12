# Lemma M numerical falsification test (lc3audit): merged multiset Z' vs Z on unfolded zeta zeros with planted clusters.
# HS norms are computed exactly (up to quadrature) on the Fourier side: sum_ij w_i w_j g(p_i-p_j) = int R(a)|sum_i w_i e(-a p_i)|^2 da,
# g = K^2 = R^, R = f*f supported in [-1,1]. Bounds use tabulated |g''''|.
# Run: bash .scratch/memcap.sh 290 900 uv run --no-project --with numpy python3 Riemann/Deep/LC3Audit/lemmaM_hs.py
import numpy as np, time
t0 = time.time(); rng = np.random.default_rng(20261112)
Z = np.load('Riemann/Deep/ExploreHelio/results/zeta_zeros.npy')
gam = Z[20000:20800]                                   # 800 zeros, t ~ 18000-18600
y0 = gam/(2*np.pi)*np.log(gam/(2*np.pi*np.e)) + 7/8    # unfolded, mean spacing 1
ap = np.array([122290,-45012,16241,1887,2711,3731,-5307,2982,-5153,5358])*1e-6
def S(x):
    x = np.clip(x, 0, 1); a = np.exp(-1/np.maximum(x, 1e-300))*(x > 0); b = np.exp(-1/np.maximum(1-x, 1e-300))*(x < 1)
    return a/(a+b)
def fwin(kind, v, eps=0.05):
    v = np.abs(v); ins = v < 0.5
    if kind == 'second': f = (1 + sum(ap[j]*np.cos(2*np.pi*(j+1)*v) for j in range(10)))*S((0.5-v)/eps)**2
    else: f = np.exp(-2/np.maximum(1-4*v*v, 1e-300))
    return f*ins
gl_x, gl_w = np.polynomial.legendre.leggauss(1200)
def norm_const(kind):
    u = 0.5*gl_x; return np.sum(0.5*gl_w*fwin(kind, u))
def Rfun(kind, a):                                     # R(a) = int f(u) f(a-u) du (f normalised)
    out = np.zeros(a.size); c = norm_const(kind)
    for i, al in enumerate(a):
        lo, hi = max(-0.5, al-0.5), min(0.5, al+0.5)
        if hi <= lo: continue
        uu = 0.5*(hi-lo)*gl_x + 0.5*(hi+lo); out[i] = 0.5*(hi-lo)*np.sum(gl_w*fwin(kind, uu)*fwin(kind, al-uu))
    return out/c**2
_x, _w = np.polynomial.legendre.leggauss(64); _e = np.linspace(-1, 1, 65)      # composite GL: 64 panels x 64 nodes on [-1, 1]
al_x = np.concatenate([0.5*(b-a)*_x + 0.5*(a+b) for a, b in zip(_e[:-1], _e[1:])]); al_w = np.concatenate([0.5*(b-a)*_w for a, b in zip(_e[:-1], _e[1:])])
# (leggauss(4096) would build a 4096^2 companion matrix: > 290 MB)
def hs2(kind, Rv, p, wts):                             # sum_ij w_i w_j g(p_i - p_j), chunked over points
    acc = np.zeros(al_x.size, complex)
    for s in range(0, p.size, 400):
        acc += np.exp(-2j*np.pi*np.outer(al_x, p[s:s+400]))@wts[s:s+400]
    return float(np.sum(al_w*Rv*np.abs(acc)**2))
def g4tab(kind, xmax, h=0.005):                        # |g''''| on [0, xmax], g = K^2, K = int f cos(2 pi u x)
    u = 0.25*(gl_x+1); w = 0.25*gl_w; c = norm_const(kind); f = 2*w*fwin(kind, u)/c; tu = 2*np.pi*u
    x = np.arange(0, xmax+h, h); out = np.zeros(x.size); out2 = np.zeros(x.size)
    for s in range(0, x.size, 1000):
        ph = np.outer(x[s:s+1000], tu); co = np.cos(ph); si = np.sin(ph)
        k0 = co@f; k1 = -(si@(f*tu)); k2 = -(co@(f*tu**2)); k3 = si@(f*tu**3); k4 = co@(f*tu**4)
        out[s:s+1000] = np.abs(2*k0*k4 + 8*k1*k3 + 6*k2**2); out2[s:s+1000] = np.abs(2*(k1**2 + k0*k2))
    return x, out, out2

def gx(kind, x):                                    # g = K^2 by direct quadrature in u (x-space check)
    u = 0.25*(gl_x+1); w = 0.25*gl_w; c = norm_const(kind); f = 2*w*fwin(kind, u)/c; tu = 2*np.pi*u; x = np.atleast_1d(x); out = np.zeros(x.size)
    for s in range(0, x.size, 1000): out[s:s+1000] = (np.cos(np.outer(x[s:s+1000], tu))@f)**2
    return out

def plant(y, kap, frac=0.06):
    kp = np.sqrt(kap); pts = [y.copy()]; idx = rng.choice(y.size, size=int(4*frac*y.size), replace=False); q = np.array_split(idx, 4)
    pts.append(y[q[0]] + rng.uniform(0.2, 0.95, q[0].size)*kap)                    # T1 tight pairs
    m = max(1, int(np.floor(kp/kap)) - 1); d = 0.97*min(kap, kp/(m+0.5))            # T2 spread clusters, diam <= kp
    for i in q[1]: pts.append(y[i] + d*np.arange(1, m+1))
    pts.append(y[q[2]]); pts.append(y[q[2]][::2] + 0.5*kap)                          # T3 double points (+ neighbour)
    n = int(np.ceil(2*kp/(0.9*kap)))                                                 # T4 long chains, extent ~ 2 kp
    for i in q[3]: pts.append(y[i] + 0.9*kap*np.arange(1, n+1))
    return np.sort(np.concatenate(pts))

def components(p, kap):
    y, m = np.unique(p, return_counts=True); br = np.where(np.diff(y) >= kap)[0] + 1
    st = np.concatenate([[0], br]); en = np.concatenate([br, [y.size]]); return y, m.astype(float), st, en

TABS = {}
def slide_max(tab, r):
    out = tab.copy()
    for k in range(1, r+1): out[:-k] = np.maximum(out[:-k], tab[k:]); out[k:] = np.maximum(out[k:], tab[:-k])
    return out

def run(kind, kap, planted=True, log=print):
    kp = np.sqrt(kap); p = plant(y0, kap) if planted else y0.copy()
    y, m, st, en = components(p, kap); N = m.sum(); diam = y[en-1] - y[st]; size = en - st
    mer = np.where((size >= 2) & (diam <= kp))[0]; lng = np.where(diam > kp)[0]
    ystar = [st[c] for c in range(st.size) if size[c] == 1 and m[st[c]] == 1]
    Lpts = [i for c in lng for i in range(st[c], en[c]) if m[i] == 1]
    # Z' and its simple points (combinatorial check of M4(ii))
    newy, newm = [], []
    for c in range(st.size):
        sl = slice(st[c], en[c])
        if c in set(mer.tolist()): newy.append(np.sum(m[sl]*y[sl])/m[sl].sum()); newm.append(m[sl].sum())
        else: newy += list(y[sl]); newm += list(m[sl])
    newy = np.array(newy); newm = np.array(newm); simple_new = set(np.round(newy[newm == 1], 12).tolist())
    simple_pred = set(np.round(y[ystar + Lpts], 12).tolist()); ok_comb = (simple_new == simple_pred) and np.unique(np.round(newy, 12)).size == newy.size
    Rv = Rfun(kind, al_x)
    # per-cluster Fourier transforms of mu_C (centroid) and mu_C^end (endpoint)
    Mh = np.zeros((al_x.size, mer.size), complex); Me = np.zeros_like(Mh); cen = np.zeros(mer.size); mC = np.zeros(mer.size)
    for k, c in enumerate(mer):
        ys, ms = y[st[c]:en[c]], m[st[c]:en[c]]; yb = np.sum(ms*ys)/ms.sum(); cen[k] = yb; mC[k] = ms.sum()
        th = 2*np.pi*np.outer(al_x, ys - yb); Mh[:, k] = np.exp(-2j*np.pi*al_x*yb)*((2j*np.sin(th/2)*np.exp(-1j*th/2))@ms)
        th = 2*np.pi*np.outer(al_x, ys - ys[0]); Me[:, k] = np.exp(-2j*np.pi*al_x*ys[0])*((2j*np.sin(th/2)*np.exp(-1j*th/2))@ms)
    E2 = float(np.sum(al_w*Rv*np.abs(Mh.sum(1))**2)) if mer.size else 0.0
    E2e = float(np.sum(al_w*Rv*np.abs(Me.sum(1))**2)) if mer.size else 0.0
    # bounds: termwise |tr D_C D_C'| <= 2 psi(Delta) m_C m_C' kp^4, B2 = sum of these, B2p = 2 kp^4 sum_{a,b merged} m_a m_b psit(y_a-y_b)
    span = y[-1] - y[0] + 2
    if TABS.get((kind, 'span'), 0) < span: TABS[kind] = g4tab(kind, span + 5); TABS[(kind, 'span')] = span + 5
    xg, tab, tab2 = TABS[kind]; h = xg[1] - xg[0]; psi2 = slide_max(tab2, int(np.ceil(2*kp/h)) + 1)
    psi = slide_max(tab, int(np.ceil(2*kp/h)) + 1); psit = slide_max(tab, int(np.ceil(4*kp/h)) + 1); Psi1 = slide_max(tab, int(np.ceil(1/h)) + 1)
    look = lambda arr, d: arr[np.minimum((np.abs(d)/h).astype(int), arr.size-1)]
    worst = 0.0; B2 = 0.0
    if mer.size:
        Tm = np.real(Mh.conj().T @ (Mh*(al_w*Rv)[:, None])); Dl = cen[:, None] - cen[None, :]
        Bm = 2*look(psi, Dl)*np.outer(mC, mC)*kp**4; B2 = float(Bm.sum())
        Be = float(np.sum(look(psi2, Dl)*np.outer(mC, mC)))*kp**2     # endpoint version: |M| <= |s||t| sup|g''| <= kp^2 sup|g''|
        log(f"   endpoint merge: exact {E2e:.4e} <= first-order bound kp^2 sum m m sup|g''| = {Be:.4e}: {E2e <= Be}")
        dg = np.sqrt(np.abs(np.diag(Tm))); tol = 1e-11*np.outer(dg, dg) + 1e-300
        worst = float(np.max((np.abs(Tm) - tol)/np.maximum(Bm, 1e-300))); worst_diag = float(np.max(np.diag(Tm)/np.diag(Bm)))
        # x-space cross-check of the Fourier quadrature on the first 25 merged clusters
        sub = list(range(min(25, mer.size))); Ex = 0.0
        for k1 in sub:
            for k2 in sub:
                c1, c2 = mer[k1], mer[k2]; ya, ma = y[st[c1]:en[c1]], m[st[c1]:en[c1]]; yb, mb = y[st[c2]:en[c2]], m[st[c2]:en[c2]]
                sa, tb, D = ya - cen[k1], yb - cen[k2], cen[k1] - cen[k2]
                S_, T_ = np.meshgrid(sa, tb, indexing='ij'); args = np.concatenate([[D], (D - T_).ravel(), (D + S_).ravel(), (D + S_ - T_).ravel()])
                gv = gx(kind, args); n = S_.size; Mv = gv[0] - gv[1:1+n] - gv[1+n:1+2*n] + gv[1+2*n:]
                Ex += float(np.sum(np.outer(ma, mb).ravel()*Mv))
        Ef = float(np.sum(Tm[np.ix_(sub, sub)]))
        log(f"   x-space check (25 clusters): {Ex:.10e} vs Fourier {Ef:.10e} | diagonal max T_CC/bound={worst_diag:.3e}")
    mp = np.concatenate([np.arange(st[c], en[c]) for c in mer]) if mer.size else np.array([], int)
    B2p = 0.0
    for s in range(0, mp.size, 200):
        a = mp[s:s+200]; B2p += float(np.sum(np.outer(m[a], m[mp])*look(psit, y[a][:, None] - y[mp][None, :])))
    B2p *= 2*kp**4
    # c_eps-type constant over all real points, pair counts, |L| bound, triangle inequality
    SPsi = 0.0
    for s in range(0, y.size, 150):
        a = np.arange(s, min(s+150, y.size)); SPsi += float(np.sum(np.outer(m[a], m)*look(Psi1, y[a][:, None] - y[None, :])))
    cum = np.concatenate([[0], np.cumsum(m)]); pair1 = float(np.sum(m*(cum[np.searchsorted(y, y+1, 'right')] - cum[np.searchsorted(y, y-1, 'left')])))
    trA = hs2(kind, Rv, y, m); trA1 = hs2(kind, Rv, newy, newm)
    tri = abs(np.sqrt(trA1) - np.sqrt(trA)) <= np.sqrt(E2)*(1 + 1e-9) + 1e-9
    log(f"{kind:6s} kap={kap:.4g} planted={planted} N={N:.0f} comps={st.size} merged={mer.size} |Y*|={len(ystar)} |L|={len(Lpts)} comb_ok={ok_comb}")
    log(f"   HS^2: exact(centroid)={E2:.4e} exact(endpoint)={E2e:.4e} | termwise max |T|/bound={worst:.3e} | E2<=B2: {E2<=B2} B2={B2:.4e} B2'={B2p:.4e}")
    log(f"   ||A'-A||/(kap sqrt N)={np.sqrt(E2)/(kap*np.sqrt(N)):.4e}  bound sqrt(2 sum Psi1)/sqrt N={np.sqrt(2*SPsi/N):.4e} | endpoint ||A'-A||/(sqrt(kap) sqrt N)={np.sqrt(E2e)/(np.sqrt(kap)*np.sqrt(N)):.4e}")
    log(f"   trA^2/N={trA/N:.6f} trA'^2/N={trA1/N:.6f} triangle_ok={tri} | |L|={len(Lpts)} <= 4 sqrt(kap) pairs(d<=1)={4*np.sqrt(kap)*pair1:.1f}: {len(Lpts) <= 4*np.sqrt(kap)*pair1} | pairs(d<=1)/N={pair1/N:.3f}")
    return dict(E2=E2, B2=B2, worst=worst)

if __name__ == '__main__':
    for kind in ['bump', 'second']:
        for kap in [1/16, 0.02, 0.005, 0.001]:
            run(kind, kap, planted=True, log=lambda s: print(s, flush=True))
        for kap in [1/16, 0.02, 0.2, 0.3]:
            run(kind, kap, planted=False, log=lambda s: print(s, flush=True))
    print(f"time {time.time()-t0:.1f}s")
