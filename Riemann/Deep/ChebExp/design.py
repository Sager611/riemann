# design (non-rigorous): chain of Lambda-intervals with kernels; output JSON for cert.py
# kernel: m=3, F(u)=e^{ss u} B_3(u;h), Phi(s)=((1-e^{-h(s-ss)})/(s-ss))^3 ; Hadamard point s0
import numpy as np, math, sys, json
from mpmath import digamma as mpdig
LP2 = math.log(math.pi / 2)  # log Delta_L >= LP2 * n_L  (Minkowski, L != Q)
LOG3 = math.log(3.0)
GAM = np.concatenate([np.linspace(0, 2, 4001), np.linspace(2.0005, 40, 6000), np.logspace(math.log10(40.001), 5, 400)])

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
PR = primes_upto(3_000_000)
XEXACT = 2.0e12   # exact prime-power sum when x <= XEXACT

def g_gamma(s0):
    a = (math.log(math.pi) - float(mpdig(s0 / 2))) / 2
    b = (math.log(2 * math.pi) - float(mpdig(s0))) / 2
    return min(a, b)

def Phi(s, h, ss):
    z = s - ss; w = h * z / 2
    with np.errstate(all='ignore'):
        sh = np.where(np.abs(w) < 1e-12, 1.0 + 0 * w, np.sinh(w) / np.where(w == 0, 1, w))
    return (h * np.exp(-w) * sh) ** 3

def F(u, h, ss):
    v = u / h
    p = lambda t: np.where(t > 0, t, 0.0) ** 2
    return np.exp(ss * u) * h * h * 0.5 * (p(v) - 3 * p(v - 1) + 3 * p(v - 2) - p(v - 3))

def consts(delta, h, ss, s0):
    P1 = Phi(np.array([1.0 + 0j]), h, ss)[0].real
    P0 = Phi(np.array([0.0 + 0j]), h, ss)[0].real
    t1 = s0 - 1 + delta; t2 = s0 - delta
    A1 = float(np.max(np.abs(Phi((1 - delta) + 1j * GAM, h, ss)) * (t1 ** 2 + GAM ** 2) / t1))
    A2 = float(np.max(np.abs(Phi(delta + 1j * GAM, h, ss)) * (t2 ** 2 + GAM ** 2) / t2))
    u = np.linspace(0, 3 * h, 6001)
    kn = float(np.max(F(u, h, ss)))
    return dict(P1=P1, P0=P0, A1=A1, A2=A2, kn=kn, U=3 * h, cB=1 / s0 + 1 / (s0 - 1), g=g_gamma(s0))

def triv_ratio(K, lx, h, ss, J1=30):
    tot = 0.0
    for j in range(1, J1 + 1):
        tot += math.exp(-j * lx) * Phi(np.array([-j + 0j]), h, ss)[0].real
    q = math.exp(3 * h - lx)
    return (tot + K['P0'] * q ** (J1 + 1) / (1 - q)) / K['P0']

def POW(lx, h, ss, kn):
    x = math.exp(lx)
    if x <= XEXACT:
        tot = 0.0; j = 2
        while 2.0 ** j <= x:
            ps = PR[PR <= x ** (1.0 / j) * (1 + 1e-12)]
            u = lx - j * np.log(ps.astype(float))
            tot += float(np.sum(np.log(ps.astype(float)) * F(u, h, ss)))
            j += 1
        return tot
    J = math.floor(lx / math.log(2)) + 1
    u = np.linspace(0, 3 * h, 6001)
    S2 = float(np.max(F(u, h, ss) * np.exp(-u / 2))); S3 = float(np.max(F(u, h, ss) * np.exp(-u / 3)))
    P2 = Phi(np.array([0.5 + 0j]), h, ss)[0].real; P3 = Phi(np.array([1 / 3 + 0j]), h, ss)[0].real
    P0 = Phi(np.array([0.0 + 0j]), h, ss)[0].real
    return 1.01624 * (math.exp(lx / 2) * (S2 + P2 / 2) + math.exp(lx / 3) * (S3 + P3 / 3) + max(0, J - 3) * math.exp(lx / 4) * (kn + P0 / 4))

def lam_b(delta, lx, h, ss, s0, margin):
    # largest Lambda_b with LHS(x, Lambda_b, n) >= margin*P1 for n in {2, n_max(Lambda_b)}; lx = log x fixed
    if h <= 0 or s0 <= 1.0001 or s0 > 6: return -1e9, None
    K = consts(delta, h, ss, s0)
    U = K['U']
    if lx <= U + 0.01: return -1e9, None
    xd = math.exp(-delta * lx); M = max(K['A1'], math.exp((2 * delta - 1) * lx) * K['A2'])
    ex = math.exp(-lx); eps = triv_ratio(K, lx, h, ss)
    Rf = min(K['kn'] * lx, K['kn'] * (U + math.log(2)), K['P0'] + K['kn'] * math.log(2)) / math.log(2)
    pw = POW(lx, h, ss, K['kn'])
    # LHS(n, Lam) = P1 - xd*M*(Lam/2 + cB - n g) - ex*[(n-1)P0 + n P0 eps + 2 kn Lam Lr/log2 + n pw]
    def coeffs(nfun_a, nfun_b):
        # n = nfun_a + nfun_b*Lam
        a0 = K['P1'] * (1 - margin) - xd * M * K['cB'] + xd * M * K['g'] * nfun_a - ex * ((nfun_a - 1) * K['P0'] + nfun_a * (K['P0'] * eps + pw))
        a1 = -xd * M / 2 + xd * M * K['g'] * nfun_b - ex * (nfun_b * K['P0'] * (1 + eps) + nfun_b * pw + 2 * Rf)
        return a0, a1
    res = []
    for (na, nb) in [(2.0, 0.0), (0.0, 1 / LP2)]:
        a0, a1 = coeffs(na, nb)
        if a1 >= 0:
            res.append(1e18 if a0 >= 0 else -1e9)
        else:
            res.append(-a0 / a1)
    return min(res), K

def tail_ok(delta, c, LamN, h, ss, s0, margin):
    if h <= 0 or s0 <= 1.0001 or s0 > 6: return -1e9
    K = consts(delta, h, ss, s0)
    ellN = LamN + LOG3; lx = math.log(c * ellN) / delta; U = K['U']
    if lx <= U + 0.5 or c * ellN < math.e: return -1e9
    M = max(K['A1'], math.exp((2 * delta - 1) * lx) * K['A2'])
    eps = triv_ratio(K, lx, h, ss)
    nmax = LamN / LP2
    ex = math.exp(-lx)
    if math.exp(lx) <= XEXACT:
        return -1e9  # tail only in the crude-power regime
    u = np.linspace(0, 3 * h, 6001)
    S2 = float(np.max(F(u, h, ss) * np.exp(-u / 2))); S3 = float(np.max(F(u, h, ss) * np.exp(-u / 3)))
    P2 = Phi(np.array([0.5 + 0j]), h, ss)[0].real; P3 = Phi(np.array([1 / 3 + 0j]), h, ss)[0].real
    powx = 1.01624 * (math.exp(-lx / 2) * (S2 + P2 / 2) + math.exp(-2 * lx / 3) * (S3 + P3 / 3) + 2 * lx * math.exp(-3 * lx / 4) * (K['kn'] + K['P0'] / 4))
    if 2 * delta * ellN / (ellN - LOG3) >= 1 or delta * (ellN / (ellN - LOG3) + 1 / math.log(c * ellN)) >= 0.75: return -1e9
    val = K['P1'] * (1 - margin) - M * (1 / (2 * c) + max(0.0, K['cB'] - LOG3 / 2 - 2 * K['g']) / (c * ellN)) \
        - ex * (nmax * K['P0'] * (1 + eps) + 2 * K['kn'] * LamN * lx / math.log(2)) - nmax * powx
    return val / K['P1']

def nm(f, x0, steps, iters):
    n = len(x0); pts = [np.array(x0, float)]
    for i in range(n):
        p = np.array(x0, float); p[i] += steps[i]; pts.append(p)
    vals = [f(p) for p in pts]
    for it in range(iters):
        o = np.argsort(vals); pts = [pts[i] for i in o]; vals = [vals[i] for i in o]
        cen = np.mean(pts[:-1], axis=0); xr = cen + (cen - pts[-1]); fr = f(xr)
        if fr < vals[0]:
            xe = cen + 2 * (cen - pts[-1]); fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]: pts[-1], vals[-1] = xr, fr
        else:
            xc = cen + 0.5 * (pts[-1] - cen); fc = f(xc)
            if fc < vals[-1]: pts[-1], vals[-1] = xc, fc
            else:
                for i in range(1, len(pts)): pts[i] = pts[0] + 0.5 * (pts[i] - pts[0]); vals[i] = f(pts[i])
    i = int(np.argmin(vals)); return pts[i], vals[i]

def chain(delta, c, Lam0, p0, margin=4e-3, maxsteps=400):
    out = []; La = Lam0; p = np.array(p0, float)
    for step in range(maxsteps):
        # try tail
        ft = lambda q: -tail_ok(delta, c, La, q[0], q[1], q[2], margin)
        qt, vt = nm(ft, p, [0.2 * p[0], 0.1, 0.2 * (p[2] - 1)], 80)
        if -vt > 0:
            out.append(dict(kind='tail', Lam_a=La, h=float(qt[0]), ss=float(qt[1]), s0=float(qt[2])))
            print(f'  tail from Lam={La:.5f} h={qt[0]:.4f} ss={qt[1]:.4f} s0={qt[2]:.4f} slack={-vt:.4f}', flush=True)
            return out
        lx = math.log(c * (La + LOG3)) / delta
        fi = lambda q: -lam_b(delta, lx, q[0], q[1], q[2], margin)[0]
        q, v = nm(fi, p, [0.2 * p[0], 0.1, 0.2 * (p[2] - 1)], 80)
        Lb = -v
        if Lb <= La + 1e-6:
            print(f'  STUCK at Lam={La:.5f} best Lb={Lb}', flush=True); return None
        Lb = min(Lb, 1e7)
        out.append(dict(kind='int', Lam_a=La, Lam_b=Lb * (1 - 1e-9), h=float(q[0]), ss=float(q[1]), s0=float(q[2])))
        print(f'  [{La:.5f}, {Lb:.5f}] h={q[0]:.4f} ss={q[1]:.4f} s0={q[2]:.4f} logx={lx:.3f}', flush=True)
        La = Lb * (1 - 1e-9); p = q
    return None

if __name__ == '__main__':
    name = sys.argv[1]; delta = eval(sys.argv[2]); c = float(sys.argv[3]); Lam0 = float(sys.argv[4])
    p0 = [float(v) for v in sys.argv[5].split(',')]
    res = chain(delta, c, Lam0, p0)
    json.dump(dict(name=name, delta_expr=sys.argv[2], c=c, Lam0=Lam0, intervals=res), open(f'chain_{name}.json', 'w'), indent=1)
    print('done', name, 'ok' if res else 'FAILED', len(res) if res else 0)
