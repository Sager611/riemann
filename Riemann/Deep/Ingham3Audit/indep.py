
# indep.py -- pbpaudit: independent float64 recomputation of the Ingham3 bound R for x-cells [La, Lb] (k-th powers).
# Own kernel envelope G, own lambda-quadrature of S and F(sigma) from their Stieltjes definitions (no moment tables, no p-grid,
# no sigma-cell freezing of exponents, no case (i)/(ii) majorants), Gauss-Legendre in sigma, best admissible KLN row at each sigma.
# KLN rows: c0, c1, q0, q1, s_min taken from cert.py's KP objects (validated there by reproducing KLN Table 1).
# Not a proof (floats); it is a check that cert.py's certified cell values are >= a direct evaluation and close to it.
import sys, math, json, time
import numpy as np
sys.path.insert(0, 'Riemann/Deep/Ingham3')
import cert
t00 = time.time()
P = [cert.KP(r) for r in cert.KLN_ROWS]
rows = [(float(p.c0.mid()), float(p.c1.mid()), float(p.q0.mid()), float(p.q1.mid()), float(p.smin.mid()), p.name) for p in P]
for (c0, c1, q0, q1, sm, nm) in rows:
    assert q1 > 0 and q0 - q1*0.8733819666132076 > 0      # B_P > 0 and decreasing on [0.515, B]
B = 0.8733819666132076; S0 = 0.515; T0 = 3e12; l0 = math.log(T0); s1 = 0.2173
def Genv(lam, dmax):
    y = lam/4.0
    with np.errstate(divide='ignore', invalid='ignore'):
        sinc = np.where(y > 0, np.sin(y)/np.where(y > 0, y, 1), 1.0)
        S = np.where(y <= np.pi, np.maximum(sinc, s1), np.minimum(s1, 1/np.where(y > 0, y, 1)))
        corr = 2*dmax*np.minimum(0.5, 1/np.where(lam > 0, lam, 1e-300))
    return (S + corr)**2
y1 = 2.0
for _ in range(200):   # sin y/y = s1 on (0, pi)
    f = math.sin(y1)/y1 - s1; y1 -= f/((y1*math.cos(y1) - math.sin(y1))/y1**2)
BPS = sorted([2.0, 4*y1, 4*math.pi, 4/s1])
def ugrid(ua, ub, n=1500):
    pts = [ua] + [math.log(b) for b in BPS if ua < math.log(b) < ub] + [ub]
    us, ws = [], []
    for a, b in zip(pts[:-1], pts[1:]):
        x, w = np.polynomial.legendre.leggauss(40)
        m = max(1, int(n*(b - a)/(ub - ua)) + 1); h = (b - a)/m
        for j in range(m):
            us.append(a + h*j + h*(x + 1)/2); ws.append(w*h/2)
    return np.concatenate(us), np.concatenate(ws)
def R_dir(k, La, Lb, dmax_list=(1e-9, 0.025, 0.15)):
    da = k*math.log1p(math.exp(-La/k)); db = k*math.log1p(math.exp(-Lb/k))
    dmax = [d for d in dmax_list if da < d][0]
    D = -math.log(db); lam0 = db*T0
    # S: 14 log 14 G(14 delta) + (1/delta) int_{14 delta}^oo G(lam)(1 + D + log lam) dlam
    u, w = ugrid(max(math.log(14*db), -60.0), math.log(1e13))
    lam = np.exp(u); G = Genv(lam, dmax)
    Sd = 14*math.log(14)*Genv(np.array([14*db]), dmax)[0] + np.sum(w*G*(1 + D + u)*lam)/db + (1 + D)*(1+dmax)**2*math.exp(-60.0)/db
    # F(sigma): int_{T0}^oo N_K(sigma,t) d(-G(delta t)) with N_K = A (log t)^{5-2s} t^{a} + B log^2 t
    u, w = ugrid(math.log(lam0), math.log(1e13)); lam = np.exp(u); G = Genv(lam, dmax); Lt = D + u   # log t
    G0 = Genv(np.array([lam0]), dmax)[0]
    def logcontrib(s):
        a = 8*(1 - s)/3; e = 5 - 2*s
        best = None
        IA = np.sum(w*G*lam**a*Lt**(e - 1)*(e + a*Lt))              # int G lam^{a-1}(..) dlam, dlam = lam du (row-independent)
        IB = np.sum(w*G*Lt)                                          # int G (D + log lam)/lam dlam
        for (c0, c1, q0, q1, sm, nm) in rows:
            if s <= sm: continue
            lA = c0 + c1*s; Bc = q0 - q1*s
            # log of A delta^{-a} IA + N_K(T0) G0 + 2 B IB
            tA1 = lA + a*D + math.log(IA)                            # A delta^{-a} I_A
            tA2 = lA + e*math.log(l0) + a*l0 + math.log(G0)          # A l0^e T0^a G(lam0)  (boundary term N_K(T0) G(delta T0))
            mx = max(tA1, tA2); tN = mx + math.log(math.exp(tA1 - mx) + math.exp(tA2 - mx))
            v = tN + math.log1p((Bc*l0**2*G0 + 2*Bc*IB)*math.exp(-tN))   # + B (l0^2 G0 + 2 I_B)
            if best is None or v < best: best = v
        return best + math.log(La) + (s - 1)*La
    x, wq = np.polynomial.legendre.leggauss(6)
    edges = [S0 + 0.0025*j for j in range(0, 144)] + [B]
    edges = [e_ for e_ in edges if e_ <= B]
    if edges[-1] != B: edges.append(B)
    tot = 0.0
    for a_, b_ in zip(edges[:-1], edges[1:]):
        ss = a_ + (b_ - a_)*(x + 1)/2
        tot += sum(wq_i*(b_ - a_)/2*math.exp(logcontrib(s_)) for s_, wq_i in zip(ss, wq))
    Z = 2*(math.exp((S0 - 1)*La)*Sd + tot)
    J = math.floor((Lb + da)/math.log(2))
    PP = 2*(Lb + da)*(math.exp(da/2 - La/2)/4 + J*math.exp(da - 2*La/3)/9 + math.log(J)*math.exp(-La - math.log(db)))
    return Z + PP + 2*math.exp(-3*La), dict(Z=Z, PP=PP, Spart=2*math.exp((S0 - 1)*La)*Sd, dmax=dmax, lam0=lam0)
if __name__ == '__main__':
    cells = json.load(open(sys.argv[1]))
    out = []
    for c in cells:
        t = time.time()
        r, det = R_dir(c['k'], c['La'], c['Lb'])
        rp, _ = R_dir(c['k'], c['La'], c['La'])
        out.append(dict(c, R_dir_cell=r, R_dir_point=rp, Z=det['Z'], PP=det['PP'], Spart=det['Spart']))
        print('%-10s k=%d [%.6f, %.6f]  R_dir(cell)=%.6g  R_dir(point La)=%.6g  Spart=%.3g PP=%.3g  %.1fs' % (c.get('tag',''), c['k'], c['La'], c['Lb'], r, rp, det['Spart'], det['PP'], time.time() - t), flush=True)
    json.dump(out, open(sys.argv[2], 'w'), indent=1)
    print('total %.1fs' % (time.time() - t00))
