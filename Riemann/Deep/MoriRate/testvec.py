# testvec.py (morirate audit): Mori 24047156 (paper III) endpoint-smoothed prolate vector k~ = E(h~)|[-a,a], rebuilt from its definition
# (sec. 1.3), and W(k~) computed from the zero side under RH (verified at these heights): W(u) = sum_rho |u^(gamma)|^2 over both signs,
# u^(z) = int u(y) e^{izy} dy. Reproduces Table 3 (mu = 5, 7: W(k)/B = 1.6, 1.8; W(k~)/B~ = 2.1, 2.3; ||k||^2/H^2 = 0.222685, 0.221439).
# Usage: testvec.py mu T   (zeros up to height T; tail beyond T from the mean of t^2 |u^(t)|^2 over the top zeros)
import sys, os, math, json, time
import mpmath as mp

def legendre_vals(u, kmax):
    v = [mp.mpf(1), u]
    for k in range(1, kmax):
        v.append(((2*k+1)*u*v[k] - k*v[k-1])/(k+1))
    return v

def prolate(c, idx, digits):
    import numpy as np
    from scipy.linalg import eigh_tridiagonal
    cf = float(c)
    kmax = int(1.6*cf + 2.5*digits) + 20
    if kmax % 2: kmax += 1
    ks = list(range(0, kmax+1, 2)); m = len(ks)
    D = [k*(k+1) + c**2*mp.mpf(2*k*(k+1) - 1)/((2*k+3)*(2*k-1)) for k in ks]
    E = [c**2*mp.mpf((k+1)*(k+2))/((2*k+3)*mp.sqrt(mp.mpf((2*k+1)*(2*k+5)))) for k in ks[:-1]]
    w, V = eigh_tridiagonal(np.array([float(x) for x in D]), np.array([float(x) for x in E]), select='i', select_range=(0, 2))
    chi = mp.mpf(w[idx]); v = [mp.mpf(x) for x in V[:, idx]]
    for it in range(10):
        a = [D[i] - chi for i in range(m)]
        cp = [mp.mpf(0)]*m; dp = [mp.mpf(0)]*m
        cp[0] = E[0]/a[0]; dp[0] = v[0]/a[0]
        for i in range(1, m):
            den = a[i] - E[i-1]*cp[i-1]
            cp[i] = E[i]/den if i < m-1 else mp.mpf(0)
            dp[i] = (v[i] - E[i-1]*dp[i-1])/den
        y = [mp.mpf(0)]*m; y[-1] = dp[-1]
        for i in range(m-2, -1, -1): y[i] = dp[i] - cp[i]*y[i+1]
        nr = mp.sqrt(mp.fsum(t*t for t in y)); v = [t/nr for t in y]
        Tv = [D[i]*v[i] + (E[i]*v[i+1] if i < m-1 else 0) + (E[i-1]*v[i-1] if i > 0 else 0) for i in range(m)]
        chi = mp.fsum(v[i]*Tv[i] for i in range(m))
    if v[0] < 0: v = [-t for t in v]
    coef = [v[j]*mp.sqrt(ks[j] + mp.mpf(1)/2) for j in range(m)]
    return ks, coef, kmax

def gl_nodes(n):
    xs, ws = [], []
    for i in range(1, n+1):
        x = mp.cos(mp.pi*(i - mp.mpf(1)/4)/(n + mp.mpf(1)/2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for k in range(1, n): p0, p1 = p1, ((2*k+1)*x*p1 - k*p0)/(k+1)
            dp = n*(x*p1 - p0)/(x*x - 1)
            dx = p1/dp; x -= dx
            if abs(dx) < mp.mpf(10)**(-mp.mp.dps + 5): break
        p0, p1 = mp.mpf(1), x
        for k in range(1, n): p0, p1 = p1, ((2*k+1)*x*p1 - k*p0)/(k+1)
        dp = n*(x*p1 - p0)/(x*x - 1)
        xs.append(x); ws.append(2/((1 - x*x)*dp*dp))
    return xs, ws

def main():
    mu = int(sys.argv[1]); T = float(sys.argv[2])
    c0 = 2*math.pi*mu
    d4est = 1210*c0**4.5*math.exp(-2*c0)
    digits = max(40, int(-math.log10(d4est)) + 22)
    mp.mp.dps = digits
    t0 = time.time()
    lam = mp.sqrt(mu); a = mp.log(lam); c = 2*mp.pi*mu
    P = {}
    for n, idx in ((0, 0), (4, 2)):
        ks, coef, kmax = prolate(c, idx, digits)
        P[n] = (ks, coef, kmax)
    def psi(n, u):
        ks, coef, kmax = P[n]
        vals = legendre_vals(u, kmax)
        return mp.fsum(coef[j]*vals[ks[j]] for j in range(len(ks)))
    d = {}; I = {}
    for n in (0, 4):
        ks, coef, kmax = P[n]
        beta0 = coef[0]/mp.sqrt(mp.mpf(1)/2)
        psi0 = psi(n, mp.mpf(0))
        lamn = c*beta0**2/(mp.pi*psi0**2)
        d[n] = 1 - lamn
        I[n] = mp.sqrt(lam)*mp.sqrt(2)*beta0
    EL = lambda u: mp.exp(-c**2*(1 - u*u)/4)
    Lf = lambda u: (1 - EL(u))**2
    DL = lambda u: 2*EL(u) - EL(u)**2
    J = {}
    for n in (0, 4):
        f = lambda u: psi(n, u)*DL(u)
        J[n] = mp.sqrt(lam)*2*mp.quad(f, [0, 1 - 200/c**2, 1 - 20/c**2, 1 - 2/c**2, 1])
    It = {n: I[n] - J[n] for n in (0, 4)}
    B = I[4]**2*d[0] + I[0]**2*d[4]
    Bt = It[4]**2*d[0] + It[0]**2*d[4]
    H2 = I[4]**2 + I[0]**2
    print(json.dumps({'mu': mu, 'dps': digits, 'c': float(c), 'd0': mp.nstr(d[0], 6), 'd4': mp.nstr(d[4], 6), 'I0': mp.nstr(I[0], 10), 'I4': mp.nstr(I[4], 10),
                      'J0': mp.nstr(J[0], 6), 'J4': mp.nstr(J[4], 6), 'B': mp.nstr(B, 8), 'Bt': mp.nstr(Bt, 8), 'H2': mp.nstr(H2, 10), 't': round(time.time()-t0, 1)}), flush=True)
    # h (unsmoothed) and h~ on |x| < lam ; u = x/lam
    def hval(x, smooth):
        u = x/lam
        if abs(u) >= 1: return mp.mpf(0)
        if smooth: return Lf(u)*(It[4]*psi(0, u) - It[0]*psi(4, u))/mp.sqrt(lam)
        return (I[4]*psi(0, u) - I[0]*psi(4, u))/mp.sqrt(lam)
    # quadrature nodes on [-a, a] with breakpoints a - log m and geometric refinement towards each right end
    bps = sorted(set([-a, a] + [a - mp.log(m) for m in range(2, mu + 1) if a - mp.log(m) > -a]))
    gx, gw = gl_nodes(24)
    hmax = mp.mpf(min(0.01, 15.0/T))
    nodes = []
    for p, q in zip(bps[:-1], bps[1:]):
        cuts = [p]
        for e in (1, 2, 3, 4, 5, 6, 7):
            r = q - mp.mpf(10)**(-e)
            if r > cuts[-1]: cuts.append(r)
        cuts.append(q)
        for s0, s1 in zip(cuts[:-1], cuts[1:]):
            npiece = int(mp.ceil((s1 - s0)/hmax))
            for j in range(npiece):
                lo = s0 + (s1 - s0)*j/npiece; hi = s0 + (s1 - s0)*(j+1)/npiece
                for xg, wg in zip(gx, gw):
                    nodes.append(((lo + hi)/2 + (hi - lo)/2*xg, (hi - lo)/2*wg))
    print(json.dumps({'nodes': len(nodes), 't': round(time.time()-t0, 1)}), flush=True)
    out = {}
    for smooth in (False, True):
        kv = []
        for (y, w) in nodes:
            ey = mp.exp(y); s = mp.mpf(0); m = 1
            while m*ey < lam:
                s += hval(m*ey, smooth); m += 1
            kv.append((y, w, mp.exp(y/2)*s))
        nk2 = mp.fsum(w*k*k for (y, w, k) in kv)
        # h norm
        if smooth:
            nh2 = mp.quad(lambda u: (Lf(u)*(It[4]*psi(0, u) - It[0]*psi(4, u)))**2, [-1, -1 + 2/c**2, -1 + 20/c**2, 0, 1 - 20/c**2, 1 - 2/c**2, 1])
        else:
            nh2 = H2
        out[smooth] = (kv, nk2, nh2)
        print(json.dumps({'smooth': smooth, 'normk2': mp.nstr(nk2, 10), 'normh2': mp.nstr(nh2, 10), 'ratio': mp.nstr(nk2/nh2, 8), 't': round(time.time()-t0, 1)}), flush=True)
    # zeros
    zf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeros_tv_dps%d.txt' % digits)
    zs = [l.strip() for l in open(zf)] if os.path.exists(zf) else []
    while not zs or mp.mpf(zs[-1]) < T:
        zs.append(mp.nstr(mp.zetazero(len(zs) + 1).imag, digits))
        if len(zs) % 100 == 0:
            with open(zf, 'w') as f: f.write('\n'.join(zs) + '\n')
    with open(zf, 'w') as f: f.write('\n'.join(zs) + '\n')
    zs = [mp.mpf(z) for z in zs if mp.mpf(z) <= T]
    print(json.dumps({'nzeros': len(zs), 'T': T, 't': round(time.time()-t0, 1)}), flush=True)
    for smooth in (False, True):
        kv, nk2, nh2 = out[smooth]
        terms = []
        for g in zs:
            re = mp.fsum(w*k*mp.cos(g*y) for (y, w, k) in kv)
            im = mp.fsum(w*k*mp.sin(g*y) for (y, w, k) in kv)
            terms.append(2*(re*re + im*im))
        S = mp.fsum(terms)
        top = zs[-min(200, len(zs)//3):]; tt = terms[-len(top):]
        A = mp.fsum(t*g*g/2 for t, g in zip(tt, top))/len(top)
        tail = 2*A*(mp.log(T/(2*mp.pi)) + 1)/(2*mp.pi*T)
        BB = Bt if smooth else B
        W = S + tail
        part = {}
        for lim in (100, 300, 1000):
            part[lim] = float(mp.fsum(t for t, g in zip(terms, zs) if g <= lim)/BB)
        print(json.dumps({'smooth': smooth, 'W': mp.nstr(W, 8), 'W_over_B': float(W/BB), 'tail_over_B': float(tail/BB), 'partial_over_B': part,
                          'Rayleigh': mp.nstr(W/nk2, 6), 't': round(time.time()-t0, 1)}), flush=True)

if __name__ == '__main__':
    main()
