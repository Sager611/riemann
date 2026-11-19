# FodgeAudit falsify.py: tests the INTERMEDIATE inequalities of Fodge's Lemma C.3 proof (Zenodo 24325329 v40, p.30-31) on exact data.
# Phase: g(x) = (t/2pi)(log(N+x+m) - log(N+x)); psi_j = g - j x. Taylor coefficients of g at the midpoint x0 in mpmath (dps 40),
# g(x0) mod 1 and g'(x0)-j in mp, the rest in double; I_j by composite 20-point Gauss-Legendre with <= 1 cycle per panel.
# Checks (ratio = LHS/RHS, a failure is a ratio > 1):
#  E1  j in E: |I_j| <= 4 (pi lam)^(-1/2)
#  N1  nonstationary: |R_j| <= (g''(a)/|al-j|^3 + g''(b)/|be-j|^3 + TV(g''/psi'^3))/(4pi^2)  (TV from a grid: lower estimate of the RHS)
#  N2  nonstationary: |R_j| <= (h lam/|al-j|^3 + h lam/|be-j|^3 + eps3/(2d^2) + h lam/d^3)/(4pi^2)
#  S1  stationary: |R_j| <= a0/(2pi^2)(u_b^-3 + |u_a|^-3) + TV(b1)/(2pi), R_j = I_j - a0 e(psi(x_j)+1/8) - B_j(b) + B_j(a)
#  S2  stationary: |R_j| <= h^1.5 lam/(2pi^2)((be-j)^-3 + (j-al)^-3) + h^1.5 Xi (b-a)/(4 pi lam)  (per-j form of the proof's sums)
#  S3  u_b >= (be-j)/sqrt(h lam), |u_a| >= (j-al)/sqrt(h lam), a0 <= lam^-1/2 (ratio = RHS/LHS of each)
#  B1  |sum_{j not in E} B_j(b)| <= 1/2 and same at a (closed form, checked against 2e5-term symmetric partial sums once per case)
#  M1  sum_stat g''(x_j)^-1/2 <= mu' L lam^1/2 + (sqrt h - 1)/(2 sqrt lam)
#  D1  decomposition: |S - sum_E I_j - main - boundary - sum_window R_j| <= tail bound of N2 outside the window (identity check)
#  C3  final: |S| <= mu' L lam^1/2 + 4.524 lam^-1/2 + 1.0004 + 2.24 lam
import sys, math, json
import numpy as np
import mpmath as mp
mp.mp.dps = 40
GLx, GLw = np.polynomial.legendre.leggauss(20)
Z3 = 1.2020569031595942

def case(N, L, Lp, m, lam0, mode=None, delta=1e-6, D=20, max_stat=40):
    assert 1 <= m < Lp <= L and N >= 10**4 and L <= N/75
    a, b = -0.5, Lp - m - 0.5
    Nm = mp.mpf(N); Np = N - 0.5
    t = mp.pi*(Nm+L)**3*lam0/m
    gp = lambda x, t=None: -(tt/(2*mp.pi))*m/((Nm+x)*(Nm+x+m))
    tt = t
    if mode in ('beta', 'alpha'):      # put be = k + 1/2 + delta (stationary j = k with be - j = 1/2 + delta) or al = k - 1/2 - delta
        x_e = b if mode == 'beta' else a
        kap = (1/(2*mp.pi))*m/((Nm+x_e)*(Nm+x_e+m))
        v = -t*kap
        k = mp.floor(v) if mode == 'beta' else mp.ceil(v)
        target = k + mp.mpf(1)/2 + delta if mode == 'beta' else k - mp.mpf(1)/2 - delta
        tt = -target/kap
    elif mode == 'eps_half':           # be = k + 1/2 - delta: eps -> 1/2 (largest closed-form term pi - 2)
        kap = (1/(2*mp.pi))*m/((Nm+b)*(Nm+b+m)); v = -t*kap; k = mp.floor(v)
        tt = -(k + mp.mpf(1)/2 - delta)/kap
    t = tt
    W = mp.pi*(Nm+L)**3/t; lam = float(m/W); h = float(((Nm+L)/(Nm-mp.mpf(1)/2))**3); e3 = 3/Np; e4 = 12/Np**2
    g = lambda x: (t/(2*mp.pi))*(mp.log(Nm+x+m) - mp.log(Nm+x))
    al = gp(mp.mpf(a)); be = gp(mp.mpf(b)); alf, bef = float(al), float(be)
    ja = int(mp.floor(al + mp.mpf(1)/2)); jb = int(mp.floor(be + mp.mpf(1)/2)); E = {ja, jb}
    x0 = (a+b)/2; X0 = Nm + mp.mpf(x0)
    K = 18
    ck = [None]*(K+1)
    for k in range(1, K+1):
        ck[k] = (t/(2*mp.pi))*(-1)**(k-1)/k*((X0+m)**(-k) - X0**(-k))
    c0 = g(mp.mpf(x0)); c0m = float(c0 - mp.floor(c0)); c1 = ck[1]
    cq = np.array([float(ck[k]) for k in range(2, K+1)])
    def Q(y):   # sum_{k>=2} c_k y^k
        acc = np.zeros_like(y)
        for c in cq[::-1]: acc = acc*y + c
        return acc*y*y
    def Qp(y):  # derivative
        acc = np.zeros_like(y)
        for k in range(K, 1, -1): acc = acc*y + k*cq[k-2]
        return acc*y
    for yy in (a-x0, b-x0, 0.37*(b-x0)):   # truncation check against mp
        ex = g(mp.mpf(x0+yy)) - c0 - c1*yy
        assert abs(float(ex) - Q(np.array([yy]))[0]) < 1e-8*(1+abs(float(ex))), 'Taylor'
    tf2 = float(t/(2*mp.pi))
    g2 = lambda x: tf2*m*(2*(N+x)+m)/((N+x)**2*(N+x+m)**2)    # g'' (works on floats and arrays)
    # exact sum
    S = mp.fsum(mp.expjpi(2*g(mp.mpf(n))) for n in range(0, Lp-m)); S = complex(S)
    # Gauss-Legendre nodes adapted to the widest tested frequency
    jlo = int(math.floor(alf)) - D; jhi = int(math.ceil(bef)) + D
    maxd = max(abs(alf-jlo), abs(bef-jhi), abs(bef-alf)) + 1
    npan = int(math.ceil((b-a)*maxd/2)) + int(b-a) + 4
    edges = np.linspace(a, b, npan+1); hw = (edges[1]-edges[0])/2
    xs = (edges[:-1,None] + hw*(GLx[None,:]+1)).ravel(); ws = np.tile(GLw*hw, npan)
    ys = xs - x0; base = c0m + Q(ys)
    c1f = lambda j: float(c1 - j)
    def Ij(j):
        ph = base + c1f(j)*ys - j*x0          # psi_j = g(x0) - j x0 + (c1 - j) y + Q(y)
        return complex(np.sum(ws*np.exp(2j*np.pi*(ph - np.floor(ph)))))
    def psi(j, x):   # psi_j(x) mod-free (c0m + (c1-j)(x-x0) + Q)
        y = np.asarray(x, dtype=float) - x0
        return c0m - j*x0 + c1f(j)*y + Q(y)
    def psip(j, x):
        y = np.asarray(x, dtype=float) - x0
        return c1f(j) + Qp(y)
    gb = g(mp.mpf(b)); ga = g(mp.mpf(a)); egb = complex(mp.expjpi(2*(gb-mp.floor(gb)))); ega = complex(mp.expjpi(2*(ga-mp.floor(ga))))
    def B(j, x, ex):   # B_j(x) = e(g(x)) e(-jx)/(2 pi i (g'(x)-j)), x half-integer: e(-jx) = (-1)^j
        return ex*(-1)**(j % 2)/(2j*math.pi*float((al if x == a else be) - j))
    R = {}; out = {k: (0.0, None) for k in ('E1','N1','N2','S1','S2','S3','B1','M1','D1','C3')}
    def rec(k, ratio, info):
        if ratio > out[k][0]: out[k] = (ratio, info)
    xgrid = np.linspace(a, b, 4001); gg_grid = g2(xgrid)
    sum_E = 0; sum_main = 0; sum_R = 0; mainF = 0.0; tailN2 = 0.0
    Xi = e3**2*(h/6 + h**1.5/2)**2 + h**0.5*(e3**2*(7*h*h/36 + h**2.5/2) + e4*(h/12 + h**1.5/3))
    js = list(range(jlo, jhi+1))
    stat = [j for j in js if alf < j < bef and j not in E]; subset = len(stat) > max_stat
    if subset:
        keep = set(stat[:15] + stat[-15:] + stat[15:-15:max(1, (len(stat)-30)//10)])
        for j in stat:   # main-term sum over ALL stationary j (no integrals needed)
            Qv = -t*m/(2*mp.pi*j); xj = float((-m + mp.sqrt(m*m + 4*Qv))/2 - Nm); mainF += g2(xj)**-0.5
        js = [j for j in js if not (alf < j < bef and j not in E) or j in keep]
    for j in js:
        if j in E:
            I = Ij(j); sum_E += I
            rec('E1', abs(I)/(4/math.sqrt(math.pi*lam)), j); continue
        if j < alf or j > bef:          # nonstationary
            I = Ij(j); Rj = I - B(j, b, egb) + B(j, a, ega); sum_R += Rj
            d = min(abs(alf-j), abs(bef-j))
            pp = psip(j, xgrid); fv = gg_grid/pp**3
            rhs1 = (g2(a)/abs(alf-j)**3 + g2(b)/abs(bef-j)**3 + np.sum(np.abs(np.diff(fv))))/(4*math.pi**2)
            rhs2 = (h*lam/abs(alf-j)**3 + h*lam/abs(bef-j)**3 + e3/(2*d*d) + h*lam/d**3)/(4*math.pi**2)
            rec('N1', abs(Rj)/rhs1, (j, d)); rec('N2', abs(Rj)/rhs2, (j, d))
        else:                           # stationary
            I = Ij(j)
            # x_j: (N+x)(N+x+m) = -t m/(2 pi j)
            Qv = -t*m/(2*mp.pi*j); xj = float((-m + mp.sqrt(m*m + 4*Qv))/2 - Nm)
            assert a < xj < b
            a0 = g2(xj)**-0.5; pj = float(psi(j, xj))
            main = a0*complex(mp.expjpi(2*(mp.mpf(pj) - mp.floor(pj) + mp.mpf(1)/8)))
            Rj = I - main - B(j, b, egb) + B(j, a, ega)
            sum_main += main; sum_R += Rj; mainF += 0 if subset else a0
            ub = math.sqrt(2*float(psi(j, b) - psi(j, xj))); ua = -math.sqrt(2*float(psi(j, a) - psi(j, xj)))
            # TV of b1(u) = 1/psi'(x) - a0/u(x) on |x - x_j| >= 0.05 (lower estimate)
            xg = xgrid[np.abs(xgrid - xj) >= 0.05]
            ug = np.sign(xg - xj)*np.sqrt(np.maximum(2*(psi(j, xg) - pj), 0)); b1 = 1/psip(j, xg) - a0/ug
            left = xg < xj; tv = np.sum(np.abs(np.diff(b1[left]))) + np.sum(np.abs(np.diff(b1[~left])))
            rhs1 = a0/(2*math.pi**2)*(ub**-3 + abs(ua)**-3) + tv/(2*math.pi)
            rhs2 = h**1.5*lam/(2*math.pi**2)*((bef-j)**-3 + (j-alf)**-3) + h**1.5*Xi*(b-a)/(4*math.pi*lam)
            rec('S1', abs(Rj)/rhs1, (j, round(bef-j, 6), round(j-alf, 6))); rec('S2', abs(Rj)/rhs2, (j, round(bef-j, 6)))
            r3 = max(((bef-j)/math.sqrt(h*lam))/ub, ((j-alf)/math.sqrt(h*lam))/abs(ua), a0*math.sqrt(lam))
            rec('S3', r3, j)
    # boundary closed form
    def bsum(x, ex, v):  # sum_{j not in E} (-1)^j/(v - j)
        eps = float(v - jb) if x == b else float(v - ja); jc = jb if x == b else ja; jo = ja if x == b else jb
        cf = (-1)**(jc % 2)*((math.pi/math.sin(math.pi*eps) - 1/eps) if eps != 0 else 0.0)
        if jo != jc: cf -= (-1)**(jo % 2)/float(v - jo)
        return cf
    Bb = egb*bsum(b, egb, be)/(2j*math.pi); Ba = ega*bsum(a, ega, al)/(2j*math.pi)
    # check closed form vs symmetric partial sums (centre 0, J = 2e5)
    jj = np.arange(-200000, 200001); jj = jj[(jj != ja) & (jj != jb)]
    ps = np.sum(((-1.0)**(jj % 2))/(bef - jj)); cfb = bsum(b, egb, be)
    assert abs(ps - cfb) < 1e-4, ('closed form', ps, cfb)
    rec('B1', max(abs(Bb), abs(Ba))/0.5, None)
    # main terms
    mu1 = 2*(N+L)**1.5/L*(1/math.sqrt(Np) - 1/math.sqrt(Np+L))
    rec('M1', mainF/(mu1*L*math.sqrt(lam) + (math.sqrt(h)-1)/(2*math.sqrt(lam))), None)
    # decomposition: tail of N2 outside window
    tail = 0.0
    for side_d0 in (abs(alf - (jlo-1)), abs((jhi+1) - bef)):
        dd = side_d0 + np.arange(0, 200000)
        tail += np.sum((3*h*lam/dd**3 + e3/(2*dd*dd)))/(4*math.pi**2)
        tail += (3*h*lam/(2*dd[-1]**2) + e3/(2*dd[-1]))/(4*math.pi**2)
    if not subset:
        resid = abs(S - sum_E - sum_main - (Bb - Ba) - sum_R)
        rec('D1', resid/(tail + 1e-9), (resid, tail))
    rhsC3 = mu1*L*math.sqrt(lam) + 4.524/math.sqrt(lam) + 1.0004 + 2.24*lam
    rec('C3', abs(S)/rhsC3, (abs(S), rhsC3))
    return dict(N=N, L=L, Lp=Lp, m=m, lam=lam, mode=mode, alpha=alf, beta=bef, nstat=max(0, int(math.floor(bef-0.5)) - int(math.ceil(alf+0.5)) + 1),
                out={k: [v[0], v[1]] for k, v in out.items()})

if __name__ == '__main__':
    cases = json.loads(sys.argv[1]) if len(sys.argv) > 1 else [[10**4, 133, 133, 1, 0.01, None]]
    for cs in cases:
        r = case(*cs)
        o = r['out']
        print(json.dumps({'N': r['N'], 'L': r['L'], 'Lp': r['Lp'], 'm': r['m'], 'lam': round(r['lam'], 6), 'mode': r['mode'], 'nstat': r['nstat'],
                          **{k: (round(v[0], 5), v[1]) for k, v in o.items()}}), flush=True)
