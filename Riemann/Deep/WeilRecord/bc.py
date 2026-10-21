# pyright: reportMissingImports=false
"""Bounded-comparison certificate for Weil positivity on [-L, L] (half-width L, mu = e^{2L}).

Framework: Liu (2026, MCOM submission) / Mori (Zenodo 2026, papers VI-VIII); every step re-derived here.
Q(f) = 2|int f cosh(x/2)|^2 - 2|int f sinh(x/2)|^2 + (1/2pi) int a(t)|F_f(t)|^2 dt - <f, C_p f>,
a(t) = Re psi(1/4+it/2) - log pi,  C_p = sum_{n<mu} Lambda(n) n^{-1/2} (S_{log n} + S_{-log n}),  S_s f(x) = 1_I(x-s) f(x-s).
 (1) ||C_p|| <= s: weighted Schur test with a positive step weight w (C_p self-adjoint, nonnegative "kernel").
 (2) gamma = s+m, M = gamma - C_p, m <= M <= b = gamma+s.  a increasing on t>=0, so a(t)-gamma >= C := a(Omega)-gamma
     for |t| >= Omega; T_tail(f) = (1/2pi) int_{|t|>Omega} (a-gamma)|F|^2 >= C <f,(I-B)f>, B = band projection
     compressed to I (0 <= I-B <= I).  Cauchy-Schwarz for the form <.,(I-B).> on each parity:
     T_tail(f) >= sum_j u_j |<h_j,f>|^2, h_j = (I-B)v_j, u_j <= C/E_j, E_j >= <v_j,(I-B)v_j> (E_0=(c+1)/(pi c^2),
     E_1=(3c^2+6c+2)/(pi c^3), c = Omega L; integration by parts).  Hence Q(f) >= <f,(M+V)f>, V = K + U, K = poles
     + band part with symbol (a-gamma) 1_{|t|<=Omega}, U = sum u_j |h_j><h_j|.  [No prolate defect needed.]
 (3) If M - tau + V >= 0 then Q(f) >= tau ||f||^2.
 (4) Legendre basis v_j, P = proj on degree < N, A = E*ME, B2 = E*M^2E, J = E*VE, G = E*M^{-1}E, D = B2 - A^2.
     R = V - EJE*: ||(I-P)RP|| <= e, (I-P)R(I-P) >= -n, ||(I-P)R(I-P)|| <= h (Mori VII Lemma 3.2, re-derived).
     Block elimination (re-derived): G^{-1} + J >= alpha + beta D  =>  M + V >= 0, alpha = e + 2e^2/(m-n),
     beta = (e+n)/m^2 + 2h^2/(m^2(m-n)).  W(X) = X + X* - X*AX + m^{-1}(I - AX - X*A + X*B2X) >= G for every X
     (W - G = R_X*(m^{-1} - M^{-1})R_X, R_X = E - MEX), so W^{-1} <= G^{-1}.  With Z = J - eps_J - alpha - beta D:
     H := W + W Z W >= 0  <=>  W^{-1} + Z >= 0 (congruence by W > 0).  All of (4) applied to M - tau.
 (5) A0 = E*C_pE, B0 = E*C_p^2E: integrands are piecewise polynomials of degree <= 2N-2 between the breakpoints
     |1 - log n/L| (u = x/L in [0,1], parity symmetry), so (N+1)-point Gauss per piece is EXACT; Legendre values in
     ball arithmetic at raised precision.  Band integrals: composite Gauss with Bernstein-ellipse bounds (quad.py).
"""
import sys, os, math, time, json, struct
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
for _d in (os.path.join(HERE, '..', 'Weil', 'cert'), os.path.join(HERE, '..', 'weil')):
    if os.path.isdir(_d): sys.path.insert(1, _d)   # weil.py, cert.py, envelope.py (Riemann/Deep/Weil/cert)
from flint import arb, arb_mat, acb, ctx, fmpz
import weil, quad
from cert import chol_mid
from cert2 import phi_ext

def log(*a):
    print(*a, flush=True)

# ---------------- setup ----------------
def vm_int(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else None
    return None

class Setup:
    def __init__(self, tag):
        self.tag = tag
        if tag == 'e3':
            self.L = arb(3)/2; self.mu_f = math.exp(3.0)
        else:
            self.L = arb(int(tag)).log()/2; self.mu_f = float(int(tag))
        self.ns = []
        for n in range(2, int(math.ceil(self.mu_f)) + 1):
            p = vm_int(n)
            if p is None: continue
            if tag == 'e3': ok = n < self.mu_f
            else: ok = n < int(tag)
            if ok: self.ns.append((n, p))
    def L_(self): return self.L if self.tag == 'e3' else arb(int(self.tag)).log()/2
    def coef(self, n, p): return arb(p).log()/arb(n).sqrt()
    def tau(self, n):
        """log n / L in u-units; exact rational 2b/a when n^a = mu^b (e.g. mu=16, n=2,4,8)."""
        if self.tag != 'e3':
            mu = int(self.tag)
            for a_ in range(1, 65):
                for b_ in range(1, 2*a_ + 1):
                    if n**a_ == mu**b_:
                        from flint import fmpq
                        return arb(fmpq(2*b_, a_))
        return arb(n).log()/self.L_()

def a_fun(t):
    return acb(arb(1)/4, t/2).digamma().real - arb.pi().log()

# ---------------- serialization ----------------
def ser(x):
    m, e = x.mid().man_exp(); r = x.rad()
    rm, re_ = r.man_exp() if r != 0 else (fmpz(0), fmpz(0))
    return (int(m), int(e), int(rm), int(re_))
def deser(t):
    m, e, rm, re_ = t
    old = ctx.prec; ctx.prec = max(old, abs(m).bit_length() + 8)
    x = arb(m)*arb(2)**e
    ctx.prec = old
    if rm != 0:
        rr = float(arb(rm)*arb(2)**re_)*(1 + 2.0**-40) + 1e-310
        x = x + arb(0, rr)
    return x
# Safe binary format (no pickle): JSON header (kinds/shapes/plain values) + length-prefixed signed integers.
def _wint(f, x):
    x = int(x); b = abs(x).to_bytes(max(1, (abs(x).bit_length() + 7)//8), 'little')
    f.write(struct.pack('<bI', -1 if x < 0 else 1, len(b))); f.write(b)
def _rint(f):
    s, n = struct.unpack('<bI', f.read(5))
    return s*int.from_bytes(f.read(n), 'little')
def save_mats(path, d):
    head = {}
    for k, v in d.items():
        if isinstance(v, arb_mat): head[k] = ['M', v.nrows(), v.ncols()]
        elif isinstance(v, arb): head[k] = ['a']
        elif isinstance(v, list) and v and isinstance(v[0], arb): head[k] = ['l', len(v)]
        else: head[k] = ['o', [str(x) for x in v] if isinstance(v, tuple) else (str(v) if isinstance(v, arb) else v)]
    try:
        with open(path, 'wb') as f:
            hb = json.dumps(head).encode(); f.write(struct.pack('<Q', len(hb))); f.write(hb)
            for k, v in d.items():
                ents = ([v[i, j] for i in range(v.nrows()) for j in range(v.ncols())] if head[k][0] == 'M'
                        else [v] if head[k][0] == 'a' else v if head[k][0] == 'l' else [])
                for x in ents:
                    for t in ser(x): _wint(f, t)
    except OSError as ex:
        raise SystemExit(f'cannot write {path}: {ex}')
def load_mats(path):
    out = {}
    try:
        with open(path, 'rb') as f:
            (hl,) = struct.unpack('<Q', f.read(8)); head = json.loads(f.read(hl).decode())
            for k, h in head.items():
                if h[0] == 'M':
                    r, c = h[1], h[2]; M = arb_mat(r, c)
                    for i in range(r):
                        for j in range(c): M[i, j] = deser(tuple(_rint(f) for _ in range(4)))
                    out[k] = M
                elif h[0] == 'a': out[k] = deser(tuple(_rint(f) for _ in range(4)))
                elif h[0] == 'l': out[k] = [deser(tuple(_rint(f) for _ in range(4))) for _ in range(h[1])]
                else: out[k] = h[1]
    except (OSError, struct.error, ValueError) as ex:
        raise SystemExit(f'cannot read {path}: {ex}')
    return out

# ---------------- (1) Schur bound ----------------
def schur(S, wfile):
    try:
        with open(wfile) as fh: ws = [float(l) for l in fh]
    except (OSError, ValueError) as ex:
        raise SystemExit(f'bad weight file {wfile}: {ex}')
    M = len(ws); assert M % 2 == 0
    w = [arb(x) for x in ws]                 # exact binary values of the floats
    assert all(x > 0 for x in ws)
    shifts = []
    for n, p in S.ns:
        t = S.tau(n); c = S.coef(n, p)
        shifts.append((t, c)); shifts.append((-t, c))
    best = arb(0); argbest = None
    for i in range(M//2, M):
        ui = arb(-1) + arb(2*i)/M; ui1 = arb(-1) + arb(2*(i+1))/M
        tot = arb(0)
        for t, c in shifts:
            lo = ui + t; hi = ui1 + t
            if hi <= -1 or lo >= 1: continue
            jl = int(math.floor(float(((lo + 1)*M/2).lower()) - 1e-6)); jh = int(math.ceil(float(((hi + 1)*M/2).upper()) + 1e-6)) - 1
            jl = max(jl, 0); jh = min(jh, M - 1)
            tot += c*max(w[jl:jh+1])
        r = arb((tot/w[i]).upper())
        if r > best: best = r; argbest = i
    log(f"schur: mu={S.tag} cells={M} sup R_w <= {best.str(10, radius=False)} (cell {argbest})")
    return best

# ---------------- (5a) A0, B0 ----------------
def legendre_vec(y, N, a_c, b_c):
    out = [arb(1), y]
    p0, p1 = out[0], out[1]
    for j in range(1, N-1):
        p0, p1 = p1, a_c[j]*y*p1 - b_c[j]*p0
        out.append(p1)
    return out[:N]

def grams(S, N, prec, chunk=256, part=None):
    """part = (k, K): only pieces i with i % K == k (partial sums, unscaled; merge with 'merge')."""
    L = S.L_()
    bprec = prec + int(math.ceil(1.28*N)) + 64
    ctx.prec = bprec
    L = S.L_()
    taus = [(S.tau(n), S.coef(n, p)) for n, p in S.ns]
    bps = [arb(0), arb(1)]
    for t, c in taus:
        bps.append(abs(1 - t))
    bps = [b for b in bps if b.lower() >= -1e-30 and b.upper() <= 1 + 1e-30]
    bps.sort(key=lambda b: float(b.mid()))
    pieces = [(bps[i], bps[i+1]) for i in range(len(bps) - 1)]
    t0 = time.time()
    ng = N + 1
    gl = [arb.legendre_p_root(ng, k, weight=True) for k in range(ng)]
    log(f"grams: mu={S.tag} N={N} prec={prec} bprec={bprec} pieces={len(pieces)} gauss({ng}) {time.time()-t0:.0f}s")
    a_c = [arb(2*j+1)/(j+1) for j in range(N)]
    b_c = [arb(j)/(j+1) for j in range(N)]
    dims = [(N+1)//2, N//2]
    A0 = [arb_mat(dims[0], dims[0]), arb_mat(dims[1], dims[1])]
    B0 = [arb_mat(dims[0], dims[0]), arb_mat(dims[1], dims[1])]
    nsteps = 0
    for ip, (al, be) in enumerate(pieces):
        if part is not None and ip % part[1] != part[0]: continue
        ctx.prec = bprec
        md = (al + be)/2; hh = (be - al)/2
        if float(hh.mid()) < 1e-40:
            log(f"  piece [{float(al.mid()):.6f},{float(be.mid()):.6f}] degenerate: skipped (true length 0; contribution 0)")
            if not (hh.lower() > -1e-60 and abs(hh) < arb(2)**(-bprec + 40)): raise ValueError("bad degenerate piece")
            continue
        act = []
        for t, c in taus:
            if md + t < 1: act.append((t, c))
            elif not (md + t > 1): raise ValueError("ambiguous + shift")
            if md - t > -1: act.append((-t, c))
            elif not (md - t < -1): raise ValueError("ambiguous - shift")
        log(f"  piece [{float(al.mid()):.6f},{float(be.mid()):.6f}] active shifts {len(act)}")
        nodes = [(md + hh*x, hh*w) for (x, w) in gl]
        for c0 in range(0, len(nodes), chunk):
            ch = nodes[c0:c0+chunk]; m = len(ch)
            Pm = [arb_mat(dims[0], m), arb_mat(dims[1], m)]
            Cm = [arb_mat(dims[0], m), arb_mat(dims[1], m)]
            Wm = [arb_mat(dims[0], m), arb_mat(dims[1], m)]
            for q, (u, w) in enumerate(ch):
                ctx.prec = bprec
                pv = legendre_vec(u, N, a_c, b_c)
                cv = [arb(0)]*N
                for t, c in act:
                    sv = legendre_vec(u + t, N, a_c, b_c)          # y = u + t, t = +-log n / L
                    cv = [cv[j] + c*sv[j] for j in range(N)]
                nsteps += N*(1 + len(act))
                ctx.prec = prec
                for j in range(N):
                    pj = pv[j] + 0; cj = cv[j] + 0
                    if not (pj.is_finite() and cj.is_finite()): raise ValueError("non-finite Legendre value")
                    Pm[j % 2][j//2, q] = pj; Cm[j % 2][j//2, q] = cj; Wm[j % 2][j//2, q] = cj*w
            ctx.prec = prec
            for p in (0, 1):
                A0[p] += Pm[p]*Wm[p].transpose()
                B0[p] += Cm[p]*Wm[p].transpose()
            log(f"    nodes {c0+m}/{len(nodes)}  {time.time()-t0:.0f}s  (legendre steps {nsteps:.3e})")
    ctx.prec = prec
    if part is not None:
        return {'A0raw_0': A0[0], 'A0raw_1': A0[1], 'B0raw_0': B0[0], 'B0raw_1': B0[1], 'N': N, 'prec': prec, 'part': part}
    return finish_grams(A0, B0, N, prec)

def finish_grams(A0, B0, N, prec):
    ctx.prec = prec
    dims = [(N+1)//2, N//2]
    res = {}
    for p in (0, 1):
        d = dims[p]
        A = arb_mat(d, d); B = arb_mat(d, d)
        for i in range(d):
            ji = 2*i + p
            for k in range(d):
                jk = 2*k + p
                sc = (arb(2*ji+1)*(2*jk+1)).sqrt()
                A[i, k] = sc*(A0[p][i, k] + A0[p][k, i])/2
                B[i, k] = sc*(B0[p][i, k] + B0[p][k, i])/2
        rA = max(float(A[i, k].rad()) for i in range(d) for k in range(d))
        rB = max(float(B[i, k].rad()) for i in range(d) for k in range(d))
        log(f"  parity {p}: max radius A0 {rA:.2e}, B0 {rB:.2e}")
        res[f'A0_{p}'] = A; res[f'B0_{p}'] = B
    res['N'] = N; res['prec'] = prec
    return res

# ---------------- (5b) band matrices ----------------
def band(S, N, Omega, gamma, prec, lmax=64, target=None, chunk=400):
    ctx.prec = prec
    L = S.L_()
    Omega = arb(Omega); gamma = arb(gamma)
    target = target or arb(2)**(-prec + 40)
    x_max = L*Omega
    bprec = prec + int(float(x_max)) + 128
    ctx.prec = bprec; L_hi = S.L_(); ctx.prec = prec          # L at Bessel precision (x = L t must be sharp)
    # plan panels: weight g = a_ext - gamma (and weight 1 for the band Gram)
    pans = quad.panels(Omega, lmax)
    plan = []; errK = arb(0); errB = arb(0)
    per = target/len(pans)
    t0 = time.time()
    for (lo, hi) in pans:
        best = None
        for rho in ((1.5, 2) if float(hi) <= 2 else (2, 3, 4)):
            if not quad.digamma_poles_clear(lo, hi, rho): continue
            c, h, A_, B_ = quad.ellipse_geom(lo, hi, rho)
            mg = quad.max_abs_on_ellipse(lambda z: phi_ext(z) - gamma, lo, hi, rho, nbox=32)
            if mg is None: continue
            Mb = (2*L/arb.pi())*(2*L*B_).exp()
            Mk = Mb*mg; Mx = Mb*max(mg, arb(1))
            n = 4
            while not (quad.gauss_err_factor(h, rho, n)*Mx < per): n += 4
            if best is None or n < best[1]: best = (rho, n, Mk, Mb, h)
        if best is None: raise ValueError('no admissible ellipse')
        rho, n, Mk, Mb, h = best
        plan.append((lo, hi, n))
        errK += quad.gauss_err_factor(h, rho, n)*Mk; errB += quad.gauss_err_factor(h, rho, n)*Mb
    nn = sum(p[2] for p in plan)
    log(f"band: mu={S.tag} N={N} Omega={float(Omega)} gamma={float(gamma.mid()):.6f} prec={prec} bprec={bprec} panels={len(plan)} nodes={nn} entry err K<={errK.str(3,radius=False)} B<={errB.str(3,radius=False)} ({time.time()-t0:.0f}s)")
    # node rounding (Cauchy-free real-line bound): normalised |F'| <= (2L/pi)(|a'| + 2L|a-gamma|) (K) and (2L/pi)(2L) (B);
    # |a'(t)| <= psi'(1/4)/2 < 9; |a(t)-gamma| <= max(|a(0)-gamma|, |a(Omega)-gamma|) on [0,Omega] (a increasing).
    supg = max(abs(a_fun(arb(0)) - gamma), abs(a_fun(Omega) - gamma))
    dK = (2*L/arb.pi())*(9 + 2*L*supg); dB = (2*L/arb.pi())*2*L
    import cert2
    nodes = []; perr = arb(0)
    for (lo, hi, n) in plan:
        hh = (hi - lo)/2; md = (hi + lo)/2; mr = 0.0
        for (x, w) in weil.gauss(n):
            rr = md + hh*x
            nodes.append((arb(rr.mid()), w*hh)); mr = max(mr, float(rr.rad()))
        perr += 2*hh*arb(mr)
    flush = 3*arb(2)**(-cert2.TINY_EXP)*Omega*(2*L/arb.pi())*(2*N)*(supg + 1)
    errK += perr*dK + flush; errB += perr*dB + flush
    dims = [(N+1)//2, N//2]
    K = [arb_mat(dims[0], dims[0]), arb_mat(dims[1], dims[1])]
    Bg = [arb_mat(dims[0], dims[0]), arb_mat(dims[1], dims[1])]
    t0 = time.time()
    for c0 in range(0, len(nodes), chunk):
        ch = nodes[c0:c0+chunk]; m = len(ch)
        Bm = [arb_mat(dims[0], m), arb_mat(dims[1], m)]
        Wm = [arb_mat(dims[0], m), arb_mat(dims[1], m)]
        Vm = [arb_mat(dims[0], m), arb_mat(dims[1], m)]
        for q, (t, w) in enumerate(ch):
            g = a_fun(t) - gamma
            js = cert2.sbessel(N, L_hi, t, bprec, a_fn=S.L_)
            for j in range(N):
                s = js[j] if (j//2) % 2 == 0 else -js[j]
                Bm[j % 2][j//2, q] = s; Vm[j % 2][j//2, q] = s*w; Wm[j % 2][j//2, q] = s*(w*g)
        for p in (0, 1):
            K[p] += Wm[p]*Bm[p].transpose(); Bg[p] += Vm[p]*Bm[p].transpose()
            if not all(K[p][i, i].is_finite() for i in range(dims[p])): raise ValueError("non-finite band matrix")
        log(f"    nodes {c0+m}/{len(nodes)} {time.time()-t0:.0f}s")
    res = {}
    for p in (0, 1):
        d = dims[p]
        Kp = arb_mat(d, d); Bp = arb_mat(d, d)
        for i in range(d):
            ji = 2*i + p
            for k in range(d):
                jk = 2*k + p
                sc = (2*L/arb.pi())*(arb(2*ji+1)*(2*jk+1)).sqrt()
                Kp[i, k] = sc*(K[p][i, k] + K[p][k, i])/2
                Bp[i, k] = sc*(Bg[p][i, k] + Bg[p][k, i])/2
        res[f'K_{p}'] = Kp; res[f'Bb_{p}'] = Bp
    Av = [(2*(2*j+1)*L).sqrt()*(arb.pi()/L).sqrt()*(L/2).bessel_i(arb(j)+arb(1)/2) for j in range(N)]
    assert all(v.is_finite() for v in Av)
    res['Av'] = Av; res['errK'] = errK; res['errB'] = errB; res['N'] = N; res['Omega'] = Omega; res['gamma'] = gamma
    return res

# ---------------- (4) certificate ----------------
def mat_from(f, d):
    M = arb_mat(d, d)
    for i in range(d):
        for k in range(d): M[i, k] = f(i, k)
    return M
def eye(d, v=1):
    M = arb_mat(d, d)
    for i in range(d): M[i, i] = arb(v)
    return M
def frob_rad(M):
    s = arb(0)
    for i in range(M.nrows()):
        for k in range(M.ncols()): s += arb(M[i, k].rad())**2
    return arb(s.sqrt().upper())
def sym(M):
    d = M.nrows()
    return mat_from(lambda i, k: (M[i, k] + M[k, i])/2, d)
def lam_min_mid(M):
    return weil.lam_min(M.mid())

def sub_mat(M, r, c):
    X = arb_mat(r, c)
    for i in range(r):
        for k in range(c): X[i, k] = M[i, k]
    return X

DIAG = False
def certify(S, G, Bd, s, m, tau, u, prec, rho_list=(2, 2.5, 3, 3.5, 4, 5, 6), G2=None, N=None):
    """G2 (optional): Gram file at a larger order N2 for the two-level bound
    G <= W2(X) = X*F + F*X - X*A2X + m^{-1}(I - X*A2F - F*A2X + X*B2X), F = [I;0], X: C^N -> C^{N2}
    (proof: Y = E2 X, G = Y*E + E*Y - Y*MY + R*M^{-1}R with R = E - MY, and M^{-1} <= m^{-1})."""
    ctx.prec = prec
    L = S.L_()
    N = N or min(G['N'], Bd['N'])
    assert N <= G['N'] and N <= Bd['N']
    Omega = Bd['Omega']; gamma = Bd['gamma']
    s = arb(s); m = arb(m); tau = arb(tau)
    assert gamma >= s + m, "gamma must be >= s + m"
    m_eff = gamma - s                                   # M >= gamma - s
    b = gamma + s
    c = Omega*L
    C = a_fun(Omega) - gamma
    assert C > 0, "C <= 0"
    E0 = (c + 1)/(arb.pi()*c**2); E1 = (3*c**2 + 6*c + 2)/(arb.pi()*c**3)
    u0, u1 = arb(u[0]), arb(u[1])
    assert u0*E0 <= C and u1*E1 <= C, "u too large"
    # projection errors (Lemma 3.2), choose rho
    K0 = (Omega/arb.pi())*max(abs(a_fun(arb(0)) - gamma), abs(a_fun(Omega) - gamma))
    Hp = (2*L).sqrt()*(L/2).cosh()
    qN = (2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(N).fac()
    mt = m_eff - tau
    best = None
    for rho in rho_list:
        rho = arb(rho)
        rN = ((8*L*rho**2/(rho - 1)**2)*(c*(rho - 1/rho)).exp()*rho**(-2*N)).sqrt()
        e = (K0*(2*L).sqrt() + (Omega/arb.pi()).sqrt()*(u0*E0.sqrt() + u1*E1.sqrt()))*rN + 4*Hp*qN
        n = K0*rN**2 + 4*qN**2
        h = (K0 + (Omega/arb.pi())*(u0 + u1))*rN**2 + 4*qN**2
        if not (n < mt): continue
        alpha = e + 2*e**2/(mt - n)
        beta = (e + n)/mt**2 + 2*h**2/(mt**2*(mt - n))
        tot = alpha + beta*b**2
        if best is None or tot < best[0]: best = (tot, rho, e, n, h, alpha, beta)
    if best is None: raise ValueError('no admissible rho (n >= m)')
    tot, rho, e, n, h, alpha, beta = best
    log(f"cert: mu={S.tag} L={L.str(8,radius=False)} N={N} Omega={float(Omega.mid())} c={float(c.mid()):.2f} gamma={float(gamma.mid()):.6f} s={float(s.mid()):.6f} m_eff={float(m_eff.mid()):.6f} b={float(b.mid()):.4f} C={C.str(6,radius=False)}")
    log(f"  u=({u[0]},{u[1]}) E0={E0.str(4,radius=False)} E1={E1.str(4,radius=False)} rho={float(rho.mid())} e={e.str(3,radius=False)} n={n.str(3,radius=False)} h={h.str(3,radius=False)} alpha={alpha.str(3,radius=False)} beta={beta.str(3,radius=False)}")
    out = {}
    T0 = time.time()
    for p in (0, 1):
        dN = (N + 1 - p)//2
        A0 = sub_mat(G[f'A0_{p}'], dN, dN); B0 = sub_mat(G[f'B0_{p}'], dN, dN)
        Kb = sub_mat(Bd[f'K_{p}'], dN, dN); Bb = sub_mat(Bd[f'Bb_{p}'], dN, dN)
        d = dN
        idx = list(range(p, N, 2))
        Av = [Bd['Av'][j] for j in idx]
        sg = 2 if p == 0 else -2
        # J = band + poles + U ; errors
        ev = (Bb.nrows())*Bd['errB']                    # ||quadrature error of column of Bb|| <= sqrt(d)*entry <= d*entry
        Eh = [(1 if i == 0 else 0) - Bb[i, 0] for i in range(d)]
        uu = u0 if p == 0 else u1
        J = mat_from(lambda i, k: Kb[i, k] + sg*Av[i]*Av[k] + uu*Eh[i]*Eh[k], d)
        epsJ = arb(d)*Bd['errK'] + uu*(2*ev + ev**2)
        A = mat_from(lambda i, k: (gamma - tau if i == k else 0) - A0[i, k], d)            # E*(M - tau)E
        B2 = mat_from(lambda i, k: ((gamma - tau)**2 if i == k else 0) - 2*(gamma - tau)*A0[i, k] + B0[i, k], d)
        D = sym(B2 - A*A)
        I = eye(d)
        if G2 is None:
            # X = (B2 - mt A)^{-1} (A - mt I)  (midpoint; any X is valid)
            Am = A.mid(); B2m = B2.mid()
            X = (B2m - Am*mt.mid()).solve(Am - I*mt.mid()).mid()
            Xt = X.transpose()
            W = sym(X + Xt - Xt*A*X + (I - A*X - Xt*A + Xt*B2*X)*(1/mt))
        else:
            A0b = G2[f'A0_{p}']; B0b = G2[f'B0_{p}']; d2 = A0b.nrows()
            assert sub_mat(A0b, d, d).overlaps(A0) if hasattr(A0b, 'overlaps') else True
            Ab = mat_from(lambda i, k: (gamma - tau if i == k else 0) - A0b[i, k], d2)
            Bb2 = mat_from(lambda i, k: ((gamma - tau)**2 if i == k else 0) - 2*(gamma - tau)*A0b[i, k] + B0b[i, k], d2)
            F = arb_mat(d2, d)
            for i in range(d): F[i, i] = 1
            AF = Ab*F
            X = (Bb2.mid() - Ab.mid()*mt.mid()).solve(AF.mid() - F*mt.mid()).mid()
            Xt = X.transpose()
            W = sym(Xt*F + F.transpose()*X - Xt*Ab*X + (I - Xt*AF - AF.transpose()*X + Xt*Bb2*X)*(1/mt))
        Z = sym(J - D*beta - I*(epsJ + alpha))
        H = sym(W + W*Z*W)
        rH = frob_rad(H)
        Hm = H.mid()
        log(f"  parity {p}: d={d} epsJ={epsJ.str(3,radius=False)} ||rad H||_F<={rH.str(3,radius=False)}  ({time.time()-T0:.0f}s)")
        if DIAG:
            log(f"    diag: lam_min(A+J)~{lam_min_mid(sym(A.mid() + J.mid())).str(5,radius=False)} lam_min(T)~{lam_min_mid(sym(W.mid().inv() + Z.mid())).str(5,radius=False)}")
        # sigma search: largest sigma = 2^-k for which the midpoint Cholesky of H - sigma exists and certifies
        ok = False
        k0 = int(os.environ.get('BC_K0', '0'))      # first exponent to try (speed only; the check is rigorous)
        for k in range(k0, 2*prec, 4):
            sigma = arb(2)**(-k)
            Hs = mat_from(lambda i, kk: Hm[i, kk] - (sigma if i == kk else 0), d)
            Lm = chol_mid(Hs)
            if Lm is None: continue
            Res = Hs - Lm*Lm.transpose()
            res2 = arb(0)
            for i in range(d):
                for kk in range(d): res2 += arb(abs(Res[i, kk]).upper())**2
            bound = res2.sqrt() + rH
            ok = bool(bound < sigma)
            log(f"  parity {p}: H - 2^-{k} = LL^T + Res, ||Res||+||rad H|| <= {bound.str(3,radius=False)} -> {'PSD CERTIFIED' if ok else 'FAIL'}  ({time.time()-T0:.0f}s)")
            break
        out[p] = ok
    if out.get(0) and out.get(1):
        log(f"RESULT mu={S.tag}: CERTIFIED  Q(f) >= {tau.str(6,radius=False)} ||f||^2 for every f supported in [-L, L], L = {L.str(10,radius=False)}")
        return True
    log(f"RESULT mu={S.tag}: not certified"); return False

if __name__ == "__main__":
    cmd, tag = sys.argv[1], sys.argv[2]
    S = Setup(tag)
    if cmd == 'schur':
        ctx.prec = 128
        schur(S, sys.argv[3])
    elif cmd == 'grams':
        N, prec = int(sys.argv[3]), int(sys.argv[4])
        if len(sys.argv) > 5:
            k, K = [int(x) for x in sys.argv[5].split('/')]
            res = grams(S, N, prec, part=(k, K))
            save_mats(os.path.join(HERE, f"gramspart_{tag}_{N}_{prec}_{k}of{K}.bin"), res)
        else:
            res = grams(S, N, prec)
            save_mats(os.path.join(HERE, f"grams_{tag}_{N}_{prec}.bin"), res)
    elif cmd == 'merge':
        N, prec, K = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        ctx.prec = prec
        def load_part(k):
            pp = load_mats(os.path.join(HERE, f"gramspart_{tag}_{N}_{prec}_{k}of{K}.bin"))
            assert int(pp['part'][0]) == k and int(pp['part'][1]) == K and int(pp['N']) == N
            return pp
        p0 = load_part(0)
        A0 = [p0['A0raw_0'], p0['A0raw_1']]; B0 = [p0['B0raw_0'], p0['B0raw_1']]
        del p0
        for k in range(1, K):                   # accumulate one part at a time (memory)
            pp = load_part(k)
            A0 = [A0[0] + pp['A0raw_0'], A0[1] + pp['A0raw_1']]; B0 = [B0[0] + pp['B0raw_0'], B0[1] + pp['B0raw_1']]
            del pp
        res = finish_grams(A0, B0, N, prec)
        save_mats(os.path.join(HERE, f"grams_{tag}_{N}_{prec}.bin"), res)
    elif cmd == 'band':
        N, Om, gam, prec = int(sys.argv[3]), sys.argv[4], sys.argv[5], int(sys.argv[6])
        ctx.prec = prec
        from flint import fmpq
        gam_a = arb(fmpq(*[int(x) for x in gam.split('/')])) if '/' in gam else arb(gam)
        res = band(S, N, arb(Om), gam_a, prec)
        save_mats(os.path.join(HERE, f"band_{tag}_{N}_{Om}_{gam.replace('/','o')}_{prec}.bin"), res)
    elif cmd == 'cert':
        gfile, bfile, s, m, tau, u0, u1, prec = sys.argv[3:11]
        g2file = sys.argv[11] if len(sys.argv) > 11 and sys.argv[11] != '-' else None
        Nuse = int(sys.argv[12]) if len(sys.argv) > 12 else None
        DIAG = bool(int(os.environ.get('BC_DIAG', '0')))
        ctx.prec = int(prec)
        from flint import fmpq
        def q(x): return arb(fmpq(*[int(y) for y in x.split('/')])) if '/' in x else arb(x)
        G = load_mats(gfile); Bd = load_mats(bfile)
        G2 = load_mats(g2file) if g2file else None
        certify(S, G, Bd, q(s), q(m), q(tau), (int(u0), int(u1)), int(prec), G2=G2, N=Nuse)
