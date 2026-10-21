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
        self.rat = tag.startswith('L')          # tag 'L<p>o<q>': half-width L = p/q exactly (added 2026-10-20)
        if tag == 'e3':
            self.L = arb(3)/2; self.mu_f = math.exp(3.0)
        elif self.rat:
            from flint import fmpq
            p_, q_ = [int(x) for x in tag[1:].split('o')]
            self.Lq = fmpq(p_, q_); self.L = arb(self.Lq); self.mu_f = math.exp(2*p_/q_)
        else:
            self.L = arb(int(tag)).log()/2; self.mu_f = float(int(tag))
        self.ns = []
        for n in range(2, int(math.ceil(self.mu_f)) + 1):
            p = vm_int(n)
            if p is None: continue
            if tag == 'e3': ok = n < self.mu_f
            elif self.rat:
                lg = arb(n).log()
                if lg < 2*self.L: ok = True
                elif lg > 2*self.L: ok = False
                else: raise ValueError('cannot decide n < e^{2L}')
            else: ok = n < int(tag)
            if ok: self.ns.append((n, p))
    def L_(self):
        if self.tag == 'e3': return self.L
        if self.rat: return arb(self.Lq)
        return arb(int(self.tag)).log()/2
    def coef(self, n, p): return arb(p).log()/arb(n).sqrt()
    def tau(self, n):
        """log n / L in u-units; exact rational 2b/a when n^a = mu^b (e.g. mu=16, n=2,4,8)."""
        if self.tag != 'e3' and not self.rat:
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
    """Rigorous weighted Schur bound ||C_p|| <= max_i R_i, checked on EVERY cell i (no symmetry assumption).
    w = positive step function on M uniform cells of u = x/L in [-1,1] (exact binary values of the file's doubles).
    R_i = sum over shifts +-t_n of c_n * max{w_j : cell j meets (cell i + t) n [-1,1]} / w_i  >=  ess sup_{cell i} (C_p w)/w.
    Returns (bound, weight_exactly_even, M, argmax cell)."""
    try:
        with open(wfile) as fh: ws = [float(l) for l in fh]
    except (OSError, ValueError) as ex:
        raise SystemExit(f'bad weight file {wfile}: {ex}')
    M = len(ws); assert M % 2 == 0
    if not all(x > 0 and math.isfinite(x) for x in ws): raise ValueError("weight must be positive and finite")
    even = all(ws[i] == ws[M - 1 - i] for i in range(M))      # exact double comparison (informational only)
    w = [arb(x) for x in ws]                 # exact binary values of the floats
    shifts = []
    for n, p in S.ns:
        t = S.tau(n); c = S.coef(n, p)
        shifts.append((t, c)); shifts.append((-t, c))
    best = arb(0); argbest = None
    for i in range(M):                       # ALL cells (the old code checked i >= M/2 only)
        ui = arb(-1) + arb(2*i)/M; ui1 = arb(-1) + arb(2*(i+1))/M
        tot = arb(0)
        for t, c in shifts:
            lo = ui + t; hi = ui1 + t
            if hi <= -1 or lo >= 1: continue          # shifted cell provably outside I (up to a null set)
            jl = int(math.floor(float(((lo + 1)*M/2).lower()) - 1e-6)); jh = int(math.ceil(float(((hi + 1)*M/2).upper()) + 1e-6)) - 1
            jl = max(jl, 0); jh = min(jh, M - 1)
            if jl > jh: raise ValueError("empty cell range")
            tot += c*arb(max(ws[jl:jh+1]))            # exact: max of exact doubles
        r = arb((tot/w[i]).upper())
        if r > best: best = r; argbest = i
    log(f"schur: mu={S.tag} cells={M} (all {M} cells checked; weight exactly even: {even}) sup R_w <= {best.str(10, radius=False)} (cell {argbest})")
    return best, even, M, argbest

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
        if hh == 0:      # arb equality: certainly exactly 0 (coincident EXACT rational breakpoints, e.g. mu=16, n=2,8)
            log(f"  piece [{float(al.mid()):.6f},{float(be.mid()):.6f}] has exact length 0: skipped (contribution exactly 0)")
            continue
        if not (hh > 0): raise ValueError("piece of undetermined length (breakpoints not separated)")
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
    errK_parts = [errK, perr*dK, flush]; errB_parts = [errB, perr*dB, flush]     # Gauss, node rounding, flush
    errK += perr*dK + flush; errB += perr*dB + flush
    log(f"  per-entry error budget (normalised entries): K: gauss {errK_parts[0].str(3,radius=False)} + node rounding {errK_parts[1].str(3,radius=False)} + flush {flush.str(3,radius=False)} = {errK.str(3,radius=False)};  B: gauss {errB_parts[0].str(3,radius=False)} + {errB_parts[1].str(3,radius=False)} + {flush.str(3,radius=False)} = {errB.str(3,radius=False)}")
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
    res['errK_parts'] = errK_parts; res['errB_parts'] = errB_parts
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
def _s(x, d=4):
    return x.str(d, radius=False)

def certify(S, G, Bd, s, m, tau, u, prec, rho_list=(2, 2.5, 3, 3.5, 4, 5, 6), G2=None, N=None, wfile=None, summary=None):
    """G2 (optional): Gram file at a larger order N2 for the two-level bound
    G <= W2(X) = X*F + F*X - X*A2X + m^{-1}(I - X*A2F - F*A2X + X*B2X), F = [I;0], X: C^N -> C^{N2}
    (proof: Y = E2 X, G = Y*E + E*Y - Y*MY + R*M^{-1}R with R = E - MY, and M^{-1} <= m^{-1})."""
    ctx.prec = prec
    L = S.L_()
    N = N or min(G['N'], Bd['N'])
    assert N <= G['N'] and N <= Bd['N']
    Omega = Bd['Omega']; gamma = Bd['gamma']
    s = arb(s); m = arb(m); tau = arb(tau)
    summ = {'mu': S.tag, 'L': _s(L, 12), 'N': N, 'prec': prec, 'Omega': _s(Omega, 8), 'gamma': _s(gamma, 8),
            's': _s(s, 10), 'tau': _s(tau, 6), 'u': [int(u[0]), int(u[1])], 'parity': {}}
    # (H1) ||C_p|| <= s, re-verified inside this run from the weight file (every cell; no symmetry assumption)
    if wfile is not None:
        oldp = ctx.prec; ctx.prec = 128
        s_cert, w_even, Mc, _ = schur(S, wfile)
        ctx.prec = oldp
        if not (s_cert <= s): raise SystemExit(f"H1 FAILED: Schur bound {s_cert} > s = {s}")
        log(f"  H1 ||C_p|| <= {_s(s_cert, 10)} <= s = {_s(s, 10)}  [weighted Schur test, all {Mc} cells, weight exactly even: {w_even}]")
        summ['H1'] = {'schur_bound': _s(s_cert, 12), 'cells': Mc, 'all_cells_checked': True,
                      'weight_exactly_even': w_even, 'weight_file': os.path.basename(wfile)}
    else:
        log("  H1 WARNING: s is taken from the command line and NOT re-verified in this run (set BC_WFILE)")
        summ['H1'] = None
    assert gamma >= s + m, "gamma must be >= s + m"
    m_eff = gamma - s                                   # M >= gamma - s
    b = gamma + s
    c = Omega*L
    C = a_fun(Omega) - gamma
    if not (C > 0): raise SystemExit("H3 FAILED: C = a(Omega) - gamma is not > 0")
    E0 = (c + 1)/(arb.pi()*c**2); E1 = (3*c**2 + 6*c + 2)/(arb.pi()*c**3)
    u0, u1 = arb(u[0]), arb(u[1])
    if not (u0*E0 <= C and u1*E1 <= C): raise SystemExit("H4 FAILED: u_j E_j <= C violated")
    log(f"  H3 C = a(Omega) - gamma = {_s(C, 8)} > 0;  H4 u0 E0 = {_s(u0*E0, 8)} <= C, u1 E1 = {_s(u1*E1, 8)} <= C")
    summ['H3_C'] = _s(C, 10); summ['H4'] = [_s(u0*E0, 10), _s(u1*E1, 10)]
    # projection errors (Lemma 3.2), choose rho
    K0 = (Omega/arb.pi())*max(abs(a_fun(arb(0)) - gamma), abs(a_fun(Omega) - gamma))
    Hp = (2*L).sqrt()*(L/2).cosh()
    qN = (2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(N).fac()
    mt = m_eff - tau
    if not (mt > 0): raise SystemExit("H2 FAILED: m_t = gamma - s - tau is not > 0")
    best = None
    for rho in rho_list:
        rho = arb(rho)
        rN = ((8*L*rho**2/(rho - 1)**2)*(c*(rho - 1/rho)).exp()*rho**(-2*N)).sqrt()   # >= ||(I-P) e^{itx}||, |t| <= Omega
        e_band = K0*(2*L).sqrt()*rN                                    # ||(I-P) K_band P||
        e_U = (Omega/arb.pi()).sqrt()*(u0*E0.sqrt() + u1*E1.sqrt())*rN  # ||(I-P) U P||
        e_pole = 4*Hp*qN                                               # ||(I-P) K_pole P||
        e = e_band + e_U + e_pole
        n = K0*rN**2 + 4*qN**2                                         # (I-P)V(I-P) >= -n  (U part is >= 0)
        h = (K0 + (Omega/arb.pi())*(u0 + u1))*rN**2 + 4*qN**2           # ||(I-P)V(I-P)|| <= h
        if not (n < mt): continue
        alpha = e + 2*e**2/(mt - n)
        beta = (e + n)/mt**2 + 2*h**2/(mt**2*(mt - n))
        tot = alpha + beta*b**2
        if best is None or tot < best[0]: best = (tot, rho, e, n, h, alpha, beta, (rN, e_band, e_U, e_pole, K0*rN**2, 4*qN**2))
    if best is None: raise ValueError('no admissible rho (n >= m)')
    tot, rho, e, n, h, alpha, beta, (rN, e_band, e_U, e_pole, n_band, n_pole) = best
    corr = alpha + beta*(b - tau)**2       # >= ||alpha I + beta D||:  0 <= D <= ||M - tau||^2 <= (b - tau)^2
    log(f"cert: mu={S.tag} L={L.str(8,radius=False)} N={N} Omega={float(Omega.mid())} c={float(c.mid()):.2f} gamma={float(gamma.mid()):.6f} s={float(s.mid()):.6f} m_eff={float(m_eff.mid()):.6f} b={float(b.mid()):.4f} C={C.str(6,radius=False)}")
    log(f"  u=({u[0]},{u[1]}) E0={E0.str(4,radius=False)} E1={E1.str(4,radius=False)} rho={float(rho.mid())} e={e.str(3,radius=False)} n={n.str(3,radius=False)} h={h.str(3,radius=False)} alpha={alpha.str(3,radius=False)} beta={beta.str(3,radius=False)}")
    log(f"  H2 m_t = gamma - s - tau = {_s(mt, 10)} > 0;  H5 Legendre tail (rho={float(rho.mid())}): r_N = {_s(rN, 3)};  e = e_band {_s(e_band, 3)} + e_U {_s(e_U, 3)} + e_pole {_s(e_pole, 3)} = {_s(e, 3)};  n = {_s(n_band, 3)} + {_s(n_pole, 3)} = {_s(n, 3)} < m_t;  h = {_s(h, 3)}")
    log(f"     alpha = {_s(alpha, 3)}, beta = {_s(beta, 3)};  ||alpha I + beta D|| <= alpha + beta (b - tau)^2 = {_s(corr, 3)}  (vs tau = {_s(tau, 3)}; this correction is subtracted inside Z exactly)")
    summ['H2_mt'] = _s(mt, 12)
    summ['H5'] = {'rho': float(rho.mid()), 'rN': _s(rN, 4), 'e': _s(e, 4), 'e_band': _s(e_band, 4), 'e_U': _s(e_U, 4), 'e_pole': _s(e_pole, 4),
                  'n': _s(n, 4), 'h': _s(h, 4), 'alpha': _s(alpha, 4), 'beta': _s(beta, 4), 'correction_norm_bound': _s(corr, 4)}
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
        log(f"  parity {p}: d={d}  H6 quadrature: ||J_true - J_quad||_2 <= epsJ = d*errK {_s(arb(d)*Bd['errK'], 3)} + u(2 ev + ev^2) {_s(uu*(2*ev + ev**2), 3)} = {_s(epsJ, 3)} (ev = d*errB = {_s(ev, 3)});  ||rad H||_F <= {_s(rH, 3)}  ({time.time()-T0:.0f}s)")
        summ['parity'][p] = {'d': d, 'epsJ': _s(epsJ, 4), 'epsJ_band': _s(arb(d)*Bd['errK'], 4), 'epsJ_U': _s(uu*(2*ev + ev**2), 4), 'radH_F': _s(rH, 4)}
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
            resF = arb(res2.sqrt().upper())                 # >= ||Res_exact||_F >= ||Res_exact||_2
            bound = resF + rH                               # >= ||Res_exact + (H_true - mid H)||_2
            ok = bool(bound < sigma)
            log(f"  parity {p}: H7 H_mid - 2^-{k} I = L L^T + Res;  ||Res||_F + ||rad H||_F <= {_s(resF, 3)} + {_s(rH, 3)} = {_s(bound, 3)} {'<' if ok else 'NOT <'} 2^-{k} = {_s(sigma, 3)} -> {'PSD CERTIFIED' if ok else 'FAIL'}  ({time.time()-T0:.0f}s)")
            summ['parity'][p].update({'sigma': f'2^-{k}', 'resF': _s(resF, 4), 'chol_bound': _s(bound, 4), 'certified': ok})
            if ok:
                lamH = arb((sigma - bound).lower())         # lambda_min(H_true) >= sigma - ||Res + dH||_2
                Winf = arb(0); Wfro = arb(0)
                for i in range(d):
                    rs = arb(0)
                    for kk in range(d):
                        a_ = arb(abs(W[i, kk]).upper()); rs += a_; Wfro += a_**2
                    rs = arb(rs.upper())
                    if rs > Winf: Winf = rs
                Wfro = arb(Wfro.sqrt().upper())
                nW = Winf if Winf < Wfro else Wfro          # ||W||_2 <= ||W||_inf (W symmetric) and <= ||W||_F
                lamT = arb((lamH/nW**2).lower())            # T = W^-1 H W^-1  =>  lambda_min(T) >= lambda_min(H)/||W||_2^2
                log(f"  parity {p}: proved margins  lambda_min(H) >= {_s(lamH, 4)};  ||W||_2 <= {_s(nW, 4)};  lambda_min(W^-1 + Z) >= {_s(lamT, 4)} > 0")
                summ['parity'][p].update({'lambda_min_H_lower': _s(lamH, 4), 'normW_upper': _s(nW, 4), 'lambda_min_T_lower': _s(lamT, 4)})
            break
        out[p] = ok
    allok = bool(out.get(0) and out.get(1))
    summ['certified'] = allok and summ['H1'] is not None
    if summary:
        try:
            with open(summary, 'w') as fh: json.dump(summ, fh, indent=1)
        except OSError as ex:
            log(f"  cannot write summary {summary}: {ex}")
    if allok:
        if summ['H1'] is None: log("  (H1 not re-verified in this run: the RESULT below is conditional on ||C_p|| <= s)")
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
        wfile = os.environ.get('BC_WFILE'); summary = os.environ.get('BC_SUMMARY')
        import hashlib, flint as _fl
        def sha(pth):
            hs = hashlib.sha256()
            try:
                with open(pth, 'rb') as fh:
                    for blk in iter(lambda: fh.read(1 << 20), b''): hs.update(blk)
            except OSError as ex:
                raise SystemExit(f'cannot read {pth}: {ex}')
            return hs.hexdigest()[:16]
        log(f"inputs: {os.path.basename(gfile)} sha256:{sha(gfile)}  {os.path.basename(bfile)} sha256:{sha(bfile)}"
            + (f"  {os.path.basename(wfile)} sha256:{sha(wfile)}" if wfile else "") + f"  python-flint {_fl.__version__}  python {sys.version.split()[0]}")
        G = load_mats(gfile); Bd = load_mats(bfile)
        G2 = load_mats(g2file) if g2file else None
        ok = certify(S, G, Bd, q(s), q(m), q(tau), (int(u0), int(u1)), int(prec), G2=G2, N=Nuse, wfile=wfile, summary=summary)
        sys.exit(0 if ok else 1)
