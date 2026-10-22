# pyright: reportMissingImports=false
"""cert4.py -- driver for the bc.py certificate chain (PROOF.md sections 1-5), successor of cert3.py (2026-10-21).
Changes against cert3.py:
 (1) FIX (band flush): Bessel values below 2^-TINY_EXP are flushed to exact zeros and the dropped mass is charged
     to the entry error.  cert2.TINY_EXP = 400 made that charge 1.9e-112 per entry at prec 448, so
     epsJ = u0 (2 ev + ev^2) = 7.2e-106 exceeded the even margin at half-width 1.55 (the cert3 even FAIL).
     Now TINY_EXP = prec + 160 (env TINY_EXP overrides).
 (2) ev = sqrt(d) errB (2-norm of the quadrature error of a column; cert3 used d errB).
 (3) Z = J - (epsJ + alpha + beta (b - tau)^2) I.  Valid: 0 <= D = E*M_t(I-P)M_t E <= ||M_t||^2 <= (b - tau)^2,
     so Z <= J - epsJ - alpha - beta D; D is never formed (saves one product and two matrices).
 (4) W(X) = m^-1 I + X^T Q X + S X + (S X)^T, S = I - A/m, Q = B2/m - A, which is the same matrix as
     X + X^T - X^T A X + m^-1 (I - AX - X^T A + X^T B2 X) (A, B2 symmetric); three products instead of four.
 (5) Band: only K = Gram of (a - gamma) 1_{|t|<=Omega} and the FIRST column of the band Gram Bb (all that U needs);
     same panels, nodes and error budget as bc.band.  One parity of everything in memory at a time.
 (6) Grams: several grams_c part files (explicit piece lists, one or both parities) are summed per parity after
     checking same tag/N/prec and that each parity's pieces cover 0..npieces-1 exactly once.
 (7) Optional diagnostics (not part of the proof): DIAG_NEG=1 sizes a negative eigenvalue of H_mid when the
     Cholesky fails; SIGMA_K=k tries sigma = 2^-k first (skips the lambda_est factorisation); DIAG=1 estimates
     lambda_min(T) (T = W^-1 + Z, the certified lower bound) and lambda_min(A + J) (Galerkin, an upper bound) by
     inverse power iteration, to measure the slack of the chain at this N.
 (8) sigma = 2^-k is chosen below lambda_pow/2 (inverse power iteration with the explicit L^-1 of the first
     factorisation) instead of below the trace bound 1/||L^-1||_F^2 (still printed).
usage: cert4.py TAG GRAMS[,GRAMS2,...] OMEGA GAMMA WFILE S M TAU U0 U1 PREC [NUSE|-] [SUMMARY]
(S, M, TAU, GAMMA rationals p/q or decimals; OMEGA integer or rational)."""
import sys, os, time, json, struct, hashlib, resource
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, arb_mat, ctx, fmpq
import flint, bc, quad, weil, cert2
from cert3 import chol_inv

def load_keys(path, keys):
    """cert3.load_keys plus integrity checks (2026-10-21: a full disk once truncated a part file silently):
    every integer record must have sign byte +-1 and fit in the file, and the records must end exactly at EOF."""
    out = {}; size = os.path.getsize(path)
    def rint(f):
        b = f.read(5)
        if len(b) < 5: raise SystemExit(f"{path}: truncated (EOF inside a record)")
        s, n = struct.unpack('<bI', b)
        if s not in (1, -1) or n > 1 << 20: raise SystemExit(f"{path}: corrupt record at byte {f.tell() - 5}")
        d = f.read(n)
        if len(d) < n: raise SystemExit(f"{path}: truncated (EOF inside a record)")
        return s*int.from_bytes(d, 'little')
    def skip_int(f):
        b = f.read(5)
        if len(b) < 5: raise SystemExit(f"{path}: truncated (EOF inside a record)")
        s, n = struct.unpack('<bI', b)
        if s not in (1, -1) or n > 1 << 20: raise SystemExit(f"{path}: corrupt record at byte {f.tell() - 5}")
        f.seek(n, 1)
    with open(path, 'rb') as f:
        (hl,) = struct.unpack('<Q', f.read(8)); head = json.loads(f.read(hl).decode())
        for k, h in head.items():
            want = k in keys
            if h[0] == 'M':
                r, c = h[1], h[2]
                if want:
                    M = arb_mat(r, c)
                    for i in range(r):
                        for j in range(c): M[i, j] = bc.deser(tuple(rint(f) for _ in range(4)))
                    out[k] = M
                else:
                    for _ in range(4*r*c): skip_int(f)
            elif h[0] == 'a':
                if want: out[k] = bc.deser(tuple(rint(f) for _ in range(4)))
                else:
                    for _ in range(4): skip_int(f)
            elif h[0] == 'l':
                if want: out[k] = [bc.deser(tuple(rint(f) for _ in range(4))) for _ in range(h[1])]
                else:
                    for _ in range(4*h[1]): skip_int(f)
            else:
                out[k] = h[1]
        if f.tell() != size: raise SystemExit(f"{path}: records end at byte {f.tell()} but the file has {size} bytes")
    return out
log = bc.log

def rss():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return f"{r/2**30:.2f} GB" if sys.platform == 'darwin' else f"{r/2**20:.2f} GB"

def read_head(path):
    with open(path, 'rb') as f:
        (hl,) = struct.unpack('<Q', f.read(8)); return json.loads(f.read(hl).decode())

def q(x): return arb(fmpq(*[int(y) for y in x.split('/')])) if '/' in x else arb(x)

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''): h.update(b)
    return h.hexdigest()[:16]

def grams_parity(files, heads, S, p, prec):
    """Raw A0, B0 sums for parity p over the part files that carry it; coverage of the pieces checked."""
    key = f'A0raw_{p}'
    use = [(f, h) for f, h in zip(files, heads) if key in h]
    if not use: raise SystemExit(f"no grams part carries parity {p}")
    N = int(use[0][1]['N'][1]); pr = int(use[0][1]['prec'][1])
    for f, h in use:
        if int(h['N'][1]) != N or int(h['prec'][1]) != pr: raise SystemExit("grams parts disagree on N/prec")
        if 'tag' in h and h['tag'][1] != S.tag: raise SystemExit(f"grams part {f}: tag {h['tag'][1]} != {S.tag}")
    if all('pieces' in h for _, h in use):
        npc = int(use[0][1]['npieces'][1]); got = []
        for _, h in use:
            if int(h['npieces'][1]) != npc: raise SystemExit("grams parts disagree on npieces")
            got += [int(x) for x in h['pieces'][1]]
        if sorted(got) != list(range(npc)): raise SystemExit(f"parity {p}: pieces {sorted(got)} do not cover 0..{npc-1} exactly once")
        cov = f"pieces 0..{npc-1} covered once by {len(use)} part(s)"
    else:
        if len(use) != 1 or [int(x) for x in use[0][1]['part'][1]] != [0, 1]: raise SystemExit("old-format grams: need one 0/1 part")
        cov = "single 0/1 part"
    Ar = Br = None
    for f, _ in use:
        G = load_keys(f, {f'A0raw_{p}', f'B0raw_{p}'})
        Ar = G[f'A0raw_{p}'] if Ar is None else Ar + G[f'A0raw_{p}']
        Br = G[f'B0raw_{p}'] if Br is None else Br + G[f'B0raw_{p}']
        del G
    return Ar, Br, N, cov

def band_lean(S, N, Omega, gamma, prec, parity, lmax=64, chunk=400):
    """bc.band restricted to one parity: K (weight a - gamma on [0, Omega]) and the first column of Bb (weight 1).
    Panels, ellipses, node counts and the three-part entry error (Gauss + node rounding + flush) as in bc.band."""
    ctx.prec = prec
    L = S.L_(); Omega = arb(Omega); gamma = arb(gamma)
    target = arb(2)**(-prec + 40)
    x_max = L*Omega
    bprec = prec + int(float(x_max)) + 128
    ctx.prec = bprec; L_hi = S.L_(); ctx.prec = prec
    pans = quad.panels(Omega, lmax)
    plan = []; errK = arb(0); errB = arb(0); per = target/len(pans); t0 = time.time()
    for (lo, hi) in pans:
        best = None
        for rho in ((1.5, 2) if float(hi) <= 2 else (2, 3, 4)):
            if not quad.digamma_poles_clear(lo, hi, rho): continue
            c, h, A_, B_ = quad.ellipse_geom(lo, hi, rho)
            mg = quad.max_abs_on_ellipse(lambda z: cert2.phi_ext(z) - gamma, lo, hi, rho, nbox=32)
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
    supg = max(abs(bc.a_fun(arb(0)) - gamma), abs(bc.a_fun(Omega) - gamma))
    dK = (2*L/arb.pi())*(9 + 2*L*supg); dB = (2*L/arb.pi())*2*L
    nodes = []; perr = arb(0)
    for (lo, hi, n) in plan:
        hh = (hi - lo)/2; md = (hi + lo)/2; mr = 0.0
        for (x, w) in weil.gauss(n):
            rr = md + hh*x
            nodes.append((arb(rr.mid()), w*hh)); mr = max(mr, float(rr.rad()))
        perr += 2*hh*arb(mr)
    flush = 3*arb(2)**(-cert2.TINY_EXP)*Omega*(2*L/arb.pi())*(2*N)*(supg + 1)
    parts = (errK, perr*dK, flush, errB, perr*dB)
    errK = errK + perr*dK + flush; errB = errB + perr*dB + flush
    log(f"  band (parity {parity}): N={N} Omega={bc._s(Omega, 8)} gamma={bc._s(gamma, 8)} bprec={bprec} panels={len(plan)} nodes={len(nodes)} TINY_EXP={cert2.TINY_EXP}")
    log(f"    entry error (normalised): K: gauss {bc._s(parts[0], 3)} + node rounding {bc._s(parts[1], 3)} + flush {bc._s(parts[2], 3)} = {bc._s(errK, 3)};  B: gauss {bc._s(parts[3], 3)} + {bc._s(parts[4], 3)} + flush = {bc._s(errB, 3)}")
    idx = list(range(parity, N, 2)); d = len(idx)
    K = arb_mat(d, d); col = arb_mat(d, 1)
    for c0 in range(0, len(nodes), chunk):
        ch = nodes[c0:c0+chunk]; m = len(ch)
        Bm = arb_mat(d, m); Wm = arb_mat(d, m); V0 = arb_mat(m, 1)
        for qq, (t, w) in enumerate(ch):
            g = bc.a_fun(t) - gamma
            js = cert2.sbessel(N, L_hi, t, bprec, a_fn=S.L_)
            for i, j in enumerate(idx):
                s = js[j] if (j//2) % 2 == 0 else -js[j]
                Bm[i, qq] = s; Wm[i, qq] = s*(w*g)
                if i == 0: V0[qq, 0] = s*w
        K += Wm*Bm.transpose(); col += Bm*V0
        if not all(K[i, i].is_finite() for i in range(d)): raise ValueError("non-finite band matrix")
        if (c0//chunk) % 4 == 3 or c0 + m == len(nodes): log(f"    band nodes {c0+m}/{len(nodes)} {time.time()-t0:.0f}s")
    sc = [(2*L/arb.pi())*(arb(2*j + 1)).sqrt() for j in idx]
    Ks = bc.mat_from(lambda i, k: sc[i]*(arb(2*idx[k] + 1)).sqrt()*(K[i, k] + K[k, i])/2, d); del K
    b0 = [sc[i]*(arb(2*idx[0] + 1)).sqrt()*col[i, 0] for i in range(d)]      # = Bb[i, 0] (Bb symmetric in exact arithmetic)
    return Ks, b0, errK, errB

def inv_power(apply, d, iters=40):
    """Floating estimate of 1/||apply|| on its dominant eigenvector (apply = an SPD inverse); diagnostics only."""
    v = arb_mat(d, 1)
    for i in range(d): v[i, 0] = arb(1)/(i + 1)
    lam = None
    for it in range(iters):
        w = apply(v).mid()
        nrm = sum((w[i, 0]**2 for i in range(d)), arb(0)).sqrt().mid()
        v = (w*(1/nrm)).mid()
        new = (1/nrm).mid()
        if lam is not None and abs(new - lam) < 1e-6*abs(new): lam = new; break
        lam = new
    return lam

def main():
    tag, gfiles, Om_s, ga_s, wfile, s, m, tau, u0, u1, prec = sys.argv[1:12]
    Nuse = int(sys.argv[12]) if len(sys.argv) > 12 and sys.argv[12] != '-' else None
    summary = sys.argv[13] if len(sys.argv) > 13 else None
    prec = int(prec); prec_c = int(os.environ.get('PREC_C', prec + 128))
    cert2.TINY_EXP = int(os.environ.get('TINY_EXP', prec + 160))
    ctx.prec = prec
    S = bc.Setup(tag); L = S.L_()
    files = gfiles.split(','); heads = [read_head(f) for f in files]
    log(f"cert4 inputs: " + "  ".join(f"{os.path.basename(f)} sha256:{sha(f)} pieces={h.get('pieces', ['o', 'all'])[1]} parities={h.get('parities', ['o', [0, 1]])[1]}" for f, h in zip(files, heads)))
    log(f"  {os.path.basename(wfile)} sha256:{sha(wfile)}  python-flint {flint.__version__}  python {sys.version.split()[0]}  prec={prec} prec_c={prec_c} TINY_EXP={cert2.TINY_EXP}")
    s = q(s); m = q(m); tau = q(tau); u = (int(u0), int(u1)); Omega = q(Om_s); gamma = q(ga_s)
    ctx.prec = 128
    s_cert, w_even, Mc, _ = bc.schur(S, wfile)
    ctx.prec = prec
    if not (s_cert <= s): raise SystemExit(f"H1 FAILED: Schur bound {s_cert} > s = {s}")
    log(f"  H1 ||C_p|| <= {bc._s(s_cert, 10)} <= s = {bc._s(s, 10)}  [weighted Schur test, all {Mc} cells, weight exactly even: {w_even}]")
    Nmax = min(int(h['N'][1]) for h in heads); N = Nuse or Nmax; assert N <= Nmax
    summ = {'mu': tag, 'L': bc._s(L, 12), 'N': N, 'prec': prec, 'prec_c': prec_c, 'TINY_EXP': cert2.TINY_EXP, 'Omega': bc._s(Omega, 8),
            'gamma': bc._s(gamma, 8), 's': bc._s(s, 10), 'tau': bc._s(tau, 6), 'u': list(u), 'parity': {},
            'grams': [{'file': os.path.basename(f), 'pieces': h.get('pieces', ['o', None])[1], 'parities': h.get('parities', ['o', [0, 1]])[1]} for f, h in zip(files, heads)],
            'H1': {'schur_bound': bc._s(s_cert, 12), 'cells': Mc, 'all_cells_checked': True, 'weight_exactly_even': w_even, 'weight_file': os.path.basename(wfile)}}
    assert gamma >= s + m, "gamma must be >= s + m"
    m_eff = gamma - s; b = gamma + s; c = Omega*L
    C = bc.a_fun(Omega) - gamma
    if not (C > 0): raise SystemExit("H3 FAILED")
    E0 = (c + 1)/(arb.pi()*c**2); E1 = (3*c**2 + 6*c + 2)/(arb.pi()*c**3)
    u0, u1 = arb(u[0]), arb(u[1])
    if not (u0*E0 <= C and u1*E1 <= C): raise SystemExit("H4 FAILED: u_j E_j <= C violated")
    log(f"  H3 C = a(Omega) - gamma = {bc._s(C, 8)} > 0;  H4 u0 E0 = {bc._s(u0*E0, 8)} <= C, u1 E1 = {bc._s(u1*E1, 8)} <= C")
    summ['H3_C'] = bc._s(C, 10); summ['H4'] = [bc._s(u0*E0, 10), bc._s(u1*E1, 10)]
    K0 = (Omega/arb.pi())*max(abs(bc.a_fun(arb(0)) - gamma), abs(bc.a_fun(Omega) - gamma))
    Hp = (2*L).sqrt()*(L/2).cosh(); qN = (2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(N).fac()
    mt = m_eff - tau
    if not (mt > 0): raise SystemExit("H2 FAILED")
    best = None
    for rho in (1.4, 1.6, 1.8, 2, 2.2, 2.5, 3, 3.5, 4, 5, 6):
        rho = arb(rho)
        rN = ((8*L*rho**2/(rho - 1)**2)*(c*(rho - 1/rho)).exp()*rho**(-2*N)).sqrt()
        e_band = K0*(2*L).sqrt()*rN; e_U = (Omega/arb.pi()).sqrt()*(u0*E0.sqrt() + u1*E1.sqrt())*rN; e_pole = 4*Hp*qN
        e = e_band + e_U + e_pole
        n = K0*rN**2 + 4*qN**2
        h = (K0 + (Omega/arb.pi())*(u0 + u1))*rN**2 + 4*qN**2
        if not (n < mt): continue
        alpha = e + 2*e**2/(mt - n); beta = (e + n)/mt**2 + 2*h**2/(mt**2*(mt - n))
        tot = alpha + beta*b**2
        if best is None or tot < best[0]: best = (tot, rho, rN, e, n, h, alpha, beta)
    if best is None: raise SystemExit("H5 FAILED: no admissible rho")
    _, rho, rN, e, n, h, alpha, beta = best
    corr = alpha + beta*(b - tau)**2
    log(f"cert4: mu={tag} L={L.str(10, radius=False)} N={N} Omega={bc._s(Omega, 8)} c={float(c.mid()):.2f} N/c={N/float(c.mid()):.3f} gamma={bc._s(gamma, 8)} s={bc._s(s, 8)} b={float(b.mid()):.4f} tau={bc._s(tau, 3)}")
    log(f"  H2 m_t = gamma - s - tau = {bc._s(mt, 10)} > 0;  H5 (rho={float(rho.mid())}): r_N = {bc._s(rN, 3)}, e = {bc._s(e, 3)}, n = {bc._s(n, 3)} < m_t, h = {bc._s(h, 3)}, alpha = {bc._s(alpha, 3)}, beta = {bc._s(beta, 3)};  alpha + beta (b - tau)^2 = {bc._s(corr, 3)} (subtracted in Z)")
    summ['H2_mt'] = bc._s(mt, 12)
    summ['H5'] = {'rho': float(rho.mid()), 'rN': bc._s(rN, 4), 'e': bc._s(e, 4), 'n': bc._s(n, 4), 'h': bc._s(h, 4), 'alpha': bc._s(alpha, 4), 'beta': bc._s(beta, 4), 'correction_norm_bound': bc._s(corr, 4)}
    Av_all = [(2*(2*j+1)*L).sqrt()*(arb.pi()/L).sqrt()*(L/2).bessel_i(arb(j)+arb(1)/2) for j in range(N)]
    assert all(v.is_finite() for v in Av_all)
    T0 = time.time(); allok = True
    sig_first = os.environ.get('SIGMA_K'); diag_neg = os.environ.get('DIAG_NEG') == '1'
    parities = [int(x) for x in os.environ.get('PARITIES', '0,1').split(',')]
    for p in parities:
        d = (N + 1 - p)//2
        # --- grams: A = E*(M - tau)E, B2 = E*(M - tau)^2 E  (sym + Legendre scaling as bc.finish_grams)
        Ar, Br, Ng, cov = grams_parity(files, heads, S, p, prec)
        scs = [(arb(2*(2*i + p) + 1)).sqrt() for i in range(d)]
        gt = gamma - tau
        A = bc.mat_from(lambda i, k: (gt if i == k else 0) - scs[i]*scs[k]*(Ar[i, k] + Ar[k, i])/2, d)
        B2 = bc.mat_from(lambda i, k: (gt**2 if i == k else 0) - 2*gt*scs[i]*scs[k]*(Ar[i, k] + Ar[k, i])/2 + scs[i]*scs[k]*(Br[i, k] + Br[k, i])/2, d)
        del Ar, Br
        log(f"  parity {p}: d={d}  grams ({cov}) loaded  ({time.time()-T0:.0f}s, maxrss {rss()})")
        # --- X (any matrix is valid; midpoint of the optimal one)
        I = bc.eye(d)
        Am = A.mid(); B2m = B2.mid()
        X = (B2m - Am*mt.mid()).solve(Am - I*mt.mid(), algorithm="approx").mid()
        del Am, B2m
        diag = os.environ.get('DIAG') == '1'
        A_d = A.mid() if diag else None
        Sm = bc.mat_from(lambda i, k: (1 if i == k else 0) - A[i, k]/mt, d)
        Qm = bc.mat_from(lambda i, k: B2[i, k]/mt - A[i, k], d)
        del A, B2
        QX = Qm*X; del Qm
        T1 = X.transpose()*QX; del QX
        SX = Sm*X; del Sm, X
        imt = 1/mt
        W = bc.mat_from(lambda i, k: T1[i, k] + (SX[i, k] + SX[k, i]) + (imt if i == k else 0), d); del T1, SX
        W = bc.sym(W)
        Winf = arb(0); Wfro = arb(0)
        for i in range(d):
            rs = arb(0)
            for kk in range(d):
                a_ = arb(abs(W[i, kk]).upper()); rs += a_; Wfro += a_**2
            rs = arb(rs.upper())
            if rs > Winf: Winf = rs
        Wfro = arb(Wfro.sqrt().upper()); nW = Winf if Winf < Wfro else Wfro
        log(f"  parity {p}: W assembled, ||W||_2 <= {bc._s(nW, 4)}  ({time.time()-T0:.0f}s, maxrss {rss()})")
        # --- J = band + poles + U (one parity)
        Kb, b0, errK, errB = band_lean(S, N, Omega, gamma, prec, p)
        ctx.prec = prec
        Av = [Av_all[j] for j in range(p, N, 2)]
        sg = 2 if p == 0 else -2
        ev = arb(d).sqrt()*errB
        Eh = [(1 if i == 0 else 0) - b0[i] for i in range(d)]
        uu = u0 if p == 0 else u1
        epsJ = arb(d)*errK + uu*(2*ev + ev**2)
        shift = epsJ + alpha + beta*(b - tau)**2
        Z = bc.mat_from(lambda i, k: Kb[i, k] + sg*Av[i]*Av[k] + uu*Eh[i]*Eh[k] - (shift if i == k else 0), d); del Kb
        Z = bc.sym(Z)
        log(f"  parity {p}: H6 epsJ = d errK {bc._s(arb(d)*errK, 3)} + u(2 ev + ev^2) {bc._s(uu*(2*ev + ev**2), 3)} = {bc._s(epsJ, 3)} (ev = sqrt(d) errB = {bc._s(ev, 3)});  Z = J - (epsJ + alpha + beta(b-tau)^2) I  ({time.time()-T0:.0f}s, maxrss {rss()})")
        Gal = None; W_d = None
        if diag:
            Gal = bc.mat_from(lambda i, k: A_d[i, k] + Z[i, k].mid() + (shift.mid() if i == k else 0), d); del A_d
            W_d = W.mid()
        ZW = Z*W; del Z
        H = bc.sym(W + W*ZW); del ZW, W
        rH = bc.frob_rad(H); Hm = H.mid(); del H
        log(f"  parity {p}: H = W + W Z W;  ||rad H||_F <= {bc._s(rH, 3)}  ({time.time()-T0:.0f}s, maxrss {rss()})")
        summ['parity'][p] = {'d': d, 'grams': cov, 'epsJ': bc._s(epsJ, 4), 'errK': bc._s(errK, 4), 'errB': bc._s(errB, 4), 'radH_F': bc._s(rH, 4), 'normW_upper': bc._s(nW, 4)}
        ctx.prec = prec_c
        ks = []
        if sig_first: ks.append(int(sig_first))
        lam_done = False; ok = False
        while True:
            if not ks:
                if lam_done: break
                r0 = chol_inv(Hm)
                lam_done = True
                if r0 is None:
                    log(f"  parity {p}: midpoint Cholesky of H_mid fails (H_mid not numerically PD) -> FAIL  ({time.time()-T0:.0f}s)")
                    if diag_neg:
                        for kk in (300, 340, 380, 420, 460):
                            if chol_inv(bc.mat_from(lambda i, j: Hm[i, j] + (arb(2)**(-kk) if i == j else 0), d)) is not None:
                                log(f"  parity {p}: DIAG H_mid + 2^-{kk} I is numerically PD (lambda_min(H_mid) > -{bc._s(arb(2)**(-kk), 3)})"); break
                            log(f"  parity {p}: DIAG H_mid + 2^-{kk} I not PD")
                    break
                Linv = r0[1]; del r0
                fr = arb(0)
                for i in range(d):
                    for kk in range(i + 1): fr += Linv[i, kk].mid()**2
                lam_est = 1/fr
                LiT = Linv.transpose()
                lam_pow = inv_power(lambda v: LiT*(Linv*v), d)
                log(f"  parity {p}: lambda_est(H_mid) = 1/||L^-1||_F^2 = {bc._s(lam_est, 4)};  inverse power: lambda_min(H_mid) ~ {bc._s(lam_pow, 4)}  ({time.time()-T0:.0f}s)")
                summ['parity'][p]['lambda_est_H'] = bc._s(lam_est, 4); summ['parity'][p]['lambda_pow_H'] = bc._s(lam_pow, 4)
                if diag:
                    lam_T = inv_power(lambda v: W_d*(LiT*(Linv*(W_d*v))), d)
                    rg = chol_inv(Gal)
                    lam_G = None
                    if rg is not None:
                        Lg = rg[1]; LgT = Lg.transpose(); del rg
                        lam_G = inv_power(lambda v: LgT*(Lg*v), d); del Lg, LgT
                    log(f"  parity {p}: DIAG lambda_min(T) ~ {bc._s(lam_T, 4)} (T = W^-1 + Z, includes -tau and -epsJ);  lambda_min(A + J) ~ {bc._s(lam_G, 4) if lam_G is not None else 'not PD'} (Galerkin, includes -tau);  ratio {float((lam_T/lam_G).mid()) if lam_G is not None else float('nan'):.4f}  ({time.time()-T0:.0f}s)")
                    summ['parity'][p]['diag_lambda_T'] = bc._s(lam_T, 4); summ['parity'][p]['diag_lambda_AJ'] = bc._s(lam_G, 4) if lam_G is not None else None
                    del Gal, W_d; Gal = W_d = None
                del Linv, LiT
                kexp = int(-float((lam_pow/2).log().mid())/0.6931471805599453) + 1
                ks = [kexp, kexp + 2, kexp + 6]
            k = ks.pop(0)
            sigma = arb(2)**(-k)
            Hs = bc.mat_from(lambda i, kk: Hm[i, kk] - (sigma if i == kk else 0), d)
            r1 = chol_inv(Hs)
            if r1 is None: log(f"  parity {p}: Cholesky of H_mid - 2^-{k} I fails  ({time.time()-T0:.0f}s)"); del Hs; continue
            Lm = r1[0]; del r1
            Res = Hs - Lm*Lm.transpose(); del Lm, Hs
            res2 = arb(0)
            for i in range(d):
                for kk in range(d): res2 += arb(abs(Res[i, kk]).upper())**2
            del Res
            resF = arb(res2.sqrt().upper()); bound = resF + rH
            ok = bool(bound < sigma)
            log(f"  parity {p}: H7 H_mid - 2^-{k} I = L L^T + Res;  ||Res||_F + ||rad H||_F <= {bc._s(resF, 3)} + {bc._s(rH, 3)} = {bc._s(bound, 3)} {'<' if ok else 'NOT <'} 2^-{k} = {bc._s(sigma, 3)} -> {'PSD CERTIFIED' if ok else 'FAIL'}  ({time.time()-T0:.0f}s, maxrss {rss()})")
            summ['parity'][p].update({'sigma': f'2^-{k}', 'resF': bc._s(resF, 4), 'chol_bound': bc._s(bound, 4), 'certified': ok})
            if ok:
                lamH = arb((sigma - bound).lower()); lamT = arb((lamH/nW**2).lower())
                log(f"  parity {p}: proved margins  lambda_min(H) >= {bc._s(lamH, 4)};  ||W||_2 <= {bc._s(nW, 4)};  lambda_min(W^-1 + Z) >= {bc._s(lamT, 4)} > 0")
                summ['parity'][p].update({'lambda_min_H_lower': bc._s(lamH, 4), 'lambda_min_T_lower': bc._s(lamT, 4)})
                break
        ctx.prec = prec
        del Hm
        allok = allok and ok
        if summary:
            with open(summary, 'w') as fh: json.dump(summ, fh, indent=1)
    allok = allok and parities == [0, 1]
    summ['certified'] = allok
    if summary:
        with open(summary, 'w') as fh: json.dump(summ, fh, indent=1)
    if allok:
        log(f"RESULT mu={tag}: CERTIFIED  Q(f) >= {tau.str(6, radius=False)} ||f||^2 for every f supported in [-L, L], L = {L.str(10, radius=False)}")
        sys.exit(0)
    log(f"RESULT mu={tag}: not certified (parities run: {parities})"); sys.exit(1)

if __name__ == "__main__":
    main()
