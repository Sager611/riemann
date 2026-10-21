# pyright: reportMissingImports=false
"""cert3.py -- memory-lean driver for the bc.py certificate (same chain of inequalities, PROOF.md sections 1-5).
Changes against bc.certify (2026-10-20, for half-width >= 1.5, margins ~1e-97):
 (1) loads one parity at a time and frees intermediates (bc.certify kept 8 big matrices + ~12 per parity alive);
 (2) the floating-point Cholesky factor of H_mid - sigma I is computed at a raised working precision PREC_C
     (default prec + 128) by a recursive block algorithm with explicit inverse blocks (fast matrix products);
     bc's chol_mid at 384 bits leaves ||Res||_F ~ 1e-80, far above the ~1e-99 margins here.  Rigour is unchanged:
     Res = H_mid - sigma I - L L^T is enclosed in ball arithmetic, lambda_min(H_true) >= sigma - ||Res||_F - ||rad H||_F;
 (3) sigma is chosen from lambda_est = 1/||L0^-1||_F^2 (L0 L0^T ~ H_mid; a lower estimate of lambda_min(H_mid)),
     instead of a scan over 2^-k (each failed Cholesky costs a full factorisation);
 (4) H1 (weighted Schur test from the weight file, every cell) is mandatory.
usage: cert3.py TAG GRAMS BAND WFILE S M TAU U0 U1 PREC [NUSE] [SUMMARY]     (S, M, TAU rationals p/q or decimals)"""
import sys, os, time, json, struct, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, arb_mat, ctx, fmpq
import flint, bc
log = bc.log

def load_keys(path, keys):
    """bc.load_mats restricted to 'keys' (other entries are parsed and skipped, never materialised)."""
    out = {}
    def skip_int(f):
        s, n = struct.unpack('<bI', f.read(5)); f.seek(n, 1)
    with open(path, 'rb') as f:
        (hl,) = struct.unpack('<Q', f.read(8)); head = json.loads(f.read(hl).decode())
        for k, h in head.items():
            want = k in keys
            if h[0] == 'M':
                r, c = h[1], h[2]
                if want:
                    M = arb_mat(r, c)
                    for i in range(r):
                        for j in range(c): M[i, j] = bc.deser(tuple(bc._rint(f) for _ in range(4)))
                    out[k] = M
                else:
                    for _ in range(4*r*c): skip_int(f)
            elif h[0] == 'a':
                if want: out[k] = bc.deser(tuple(bc._rint(f) for _ in range(4)))
                else:
                    for _ in range(4): skip_int(f)
            elif h[0] == 'l':
                if want: out[k] = [bc.deser(tuple(bc._rint(f) for _ in range(4))) for _ in range(h[1])]
                else:
                    for _ in range(4*h[1]): skip_int(f)
            else:
                out[k] = h[1]
    return out

def blk(M, r0, r1, c0, c1):
    X = arb_mat(r1 - r0, c1 - c0)
    for i in range(r0, r1):
        for j in range(c0, c1): X[i - r0, j - c0] = M[i, j]
    return X

def chol_inv(G, leaf=48):
    """Floating-point (midpoint) Cholesky G ~ L L^T and Linv ~ L^-1 at the current ctx.prec; None if a pivot is
    not positive.  No rigour here: the caller encloses the residual."""
    n = G.nrows()
    if n <= leaf:
        A = [[G[i, j].mid() for j in range(n)] for i in range(n)]
        Lr = [[arb(0)]*n for _ in range(n)]
        for j in range(n):
            s = (A[j][j] - sum((Lr[j][k]**2 for k in range(j)), arb(0))).mid()
            if not (s > 0): return None
            Lr[j][j] = s.sqrt().mid()
            for i in range(j + 1, n):
                Lr[i][j] = ((A[i][j] - sum((Lr[i][k]*Lr[j][k] for k in range(j)), arb(0)))/Lr[j][j]).mid()
        Ir = [[arb(0)]*n for _ in range(n)]                  # forward substitution for L^-1
        for c in range(n):
            for i in range(c, n):
                t = (arb(1) if i == c else arb(0)) - sum((Lr[i][k]*Ir[k][c] for k in range(c, i)), arb(0))
                Ir[i][c] = (t/Lr[i][i]).mid()
        Lm = arb_mat(n, n); Im = arb_mat(n, n)
        for i in range(n):
            for j in range(i + 1): Lm[i, j] = Lr[i][j]; Im[i, j] = Ir[i][j]
        return Lm, Im
    k = n//2
    r1 = chol_inv(blk(G, 0, k, 0, k), leaf)
    if r1 is None: return None
    L11, I11 = r1
    L21 = (blk(G, k, n, 0, k)*I11.transpose()).mid()        # L21 = G21 L11^-T
    S = (blk(G, k, n, k, n) - L21*L21.transpose()).mid()
    r2 = chol_inv(S, leaf)
    if r2 is None: return None
    L22, I22 = r2
    I21 = (-(I22*(L21*I11))).mid()
    Lm = arb_mat(n, n); Im = arb_mat(n, n)
    for i in range(k):
        for j in range(i + 1): Lm[i, j] = L11[i, j]; Im[i, j] = I11[i, j]
    for i in range(n - k):
        for j in range(k): Lm[k + i, j] = L21[i, j]; Im[k + i, j] = I21[i, j]
        for j in range(i + 1): Lm[k + i, k + j] = L22[i, j]; Im[k + i, k + j] = I22[i, j]
    return Lm, Im

def q(x): return arb(fmpq(*[int(y) for y in x.split('/')])) if '/' in x else arb(x)

def main():
    tag, gfile, bfile, wfile, s, m, tau, u0, u1, prec = sys.argv[1:11]
    Nuse = int(sys.argv[11]) if len(sys.argv) > 11 and sys.argv[11] != '-' else None
    summary = sys.argv[12] if len(sys.argv) > 12 else None
    prec = int(prec); prec_c = int(os.environ.get('PREC_C', prec + 128))
    ctx.prec = prec
    S = bc.Setup(tag); L = S.L_()
    def sha(p):
        h = hashlib.sha256()
        with open(p, 'rb') as fh:
            for b in iter(lambda: fh.read(1 << 20), b''): h.update(b)
        return h.hexdigest()[:16]
    log(f"cert3 inputs: {os.path.basename(gfile)} sha256:{sha(gfile)}  {bfile if bfile.startswith('compute:') else os.path.basename(bfile) + ' sha256:' + sha(bfile)}  {os.path.basename(wfile)} sha256:{sha(wfile)}  python-flint {flint.__version__}  python {sys.version.split()[0]}  prec={prec} prec_c={prec_c}")
    s = q(s); m = q(m); tau = q(tau); u = (int(u0), int(u1))
    ctx.prec = 128
    s_cert, w_even, Mc, _ = bc.schur(S, wfile)
    ctx.prec = prec
    if not (s_cert <= s): raise SystemExit(f"H1 FAILED: Schur bound {s_cert} > s = {s}")
    log(f"  H1 ||C_p|| <= {bc._s(s_cert, 10)} <= s = {bc._s(s, 10)}  [weighted Schur test, all {Mc} cells, weight exactly even: {w_even}]")
    Gh = load_keys(gfile, {'N', 'part'})
    raw = 'part' in Gh                       # grams_c part file (K = 1): symmetrise and scale per parity here
    if raw: assert int(Gh['part'][0]) == 0 and int(Gh['part'][1]) == 1
    if bfile.startswith('compute:'):        # band matrices computed in this process (no band file on disk)
        _, Om_s, ga_s = bfile.split(':')
        Bfull = bc.band(S, Nuse or int(Gh['N']), q(Om_s), q(ga_s), prec); Bh = Bfull
    else:
        Bfull = None
        Bh = load_keys(bfile, {'N', 'Omega', 'gamma', 'errK', 'errB', 'Av'})
    Nmax = min(int(Gh['N']), int(Bh['N'])); N = Nuse or Nmax; assert N <= Nmax
    Omega = Bh['Omega'] if isinstance(Bh['Omega'], arb) else q(str(Bh['Omega']))
    gamma = Bh['gamma']; errK = Bh['errK']; errB = Bh['errB']; Avall = Bh['Av']
    assert isinstance(gamma, arb) and isinstance(errK, arb), "band file: gamma/errK must be balls"
    summ = {'mu': tag, 'L': bc._s(L, 12), 'N': N, 'prec': prec, 'prec_c': prec_c, 'Omega': bc._s(Omega, 8), 'gamma': bc._s(gamma, 8),
            's': bc._s(s, 10), 'tau': bc._s(tau, 6), 'u': list(u), 'parity': {},
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
    for rho in (1.6, 1.8, 2, 2.2, 2.5, 3, 3.5, 4, 5, 6):
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
    log(f"cert3: mu={tag} L={L.str(10, radius=False)} N={N} Omega={bc._s(Omega, 8)} c={float(c.mid()):.2f} gamma={bc._s(gamma, 8)} s={bc._s(s, 8)} b={float(b.mid()):.4f} tau={bc._s(tau, 3)}")
    log(f"  H2 m_t = gamma - s - tau = {bc._s(mt, 10)} > 0;  H5 (rho={float(rho.mid())}): r_N = {bc._s(rN, 3)}, e = {bc._s(e, 3)}, n = {bc._s(n, 3)} < m_t, h = {bc._s(h, 3)}, alpha = {bc._s(alpha, 3)}, beta = {bc._s(beta, 3)}, ||alpha I + beta D|| <= {bc._s(corr, 3)}")
    summ['H2_mt'] = bc._s(mt, 12)
    summ['H5'] = {'rho': float(rho.mid()), 'rN': bc._s(rN, 4), 'e': bc._s(e, 4), 'n': bc._s(n, 4), 'h': bc._s(h, 4), 'alpha': bc._s(alpha, 4), 'beta': bc._s(beta, 4), 'correction_norm_bound': bc._s(corr, 4)}
    T0 = time.time(); allok = True
    for p in (0, 1):
        dN = (N + 1 - p)//2; d = dN
        if raw:
            Gp = load_keys(gfile, {f'A0raw_{p}', f'B0raw_{p}'}); Ar = Gp[f'A0raw_{p}']; Br = Gp[f'B0raw_{p}']; del Gp
            scs = [(arb(2*(2*i + p) + 1)).sqrt() for i in range(d)]      # = bc.finish_grams for this parity
            A0 = bc.mat_from(lambda i, k: scs[i]*scs[k]*(Ar[i, k] + Ar[k, i])/2, d)
            B0 = bc.mat_from(lambda i, k: scs[i]*scs[k]*(Br[i, k] + Br[k, i])/2, d); del Ar, Br
        else:
            Gp = load_keys(gfile, {f'A0_{p}', f'B0_{p}'})
            A0 = bc.sub_mat(Gp[f'A0_{p}'], d, d); B0 = bc.sub_mat(Gp[f'B0_{p}'], d, d); del Gp
        Bp = Bfull if Bfull is not None else load_keys(bfile, {f'K_{p}', f'Bb_{p}'})
        Kb = bc.sub_mat(Bp[f'K_{p}'], d, d); Bb = bc.sub_mat(Bp[f'Bb_{p}'], d, d); del Bp
        if Bfull is not None: del Bfull[f'K_{p}'], Bfull[f'Bb_{p}']
        Av = [Avall[j] for j in range(p, N, 2)]
        sg = 2 if p == 0 else -2
        ev = Bb.nrows()*errB
        Eh = [(1 if i == 0 else 0) - Bb[i, 0] for i in range(d)]; del Bb
        uu = u0 if p == 0 else u1
        J = bc.mat_from(lambda i, k: Kb[i, k] + sg*Av[i]*Av[k] + uu*Eh[i]*Eh[k], d); del Kb
        epsJ = arb(d)*errK + uu*(2*ev + ev**2)
        A = bc.mat_from(lambda i, k: (gamma - tau if i == k else 0) - A0[i, k], d)
        B2 = bc.mat_from(lambda i, k: ((gamma - tau)**2 if i == k else 0) - 2*(gamma - tau)*A0[i, k] + B0[i, k], d)
        del A0, B0
        I = bc.eye(d)
        Am = A.mid(); B2m = B2.mid()
        X = (B2m - Am*mt.mid()).solve(Am - I*mt.mid(), algorithm="approx").mid()   # any X is valid (PROOF.md (E))
        del Am, B2m
        Xt = X.transpose()
        AX = A*X
        W = bc.sym(X + Xt - Xt*AX + (I - AX - AX.transpose() + Xt*(B2*X))*(1/mt))
        del X, Xt, AX
        D = bc.sym(B2 - A*A); del A, B2
        Winf = arb(0); Wfro = arb(0)
        for i in range(d):
            rs = arb(0)
            for kk in range(d):
                a_ = arb(abs(W[i, kk]).upper()); rs += a_; Wfro += a_**2
            rs = arb(rs.upper())
            if rs > Winf: Winf = rs
        Wfro = arb(Wfro.sqrt().upper()); nW = Winf if Winf < Wfro else Wfro
        Z = bc.sym(J - D*beta - I*(epsJ + alpha)); del J, D
        H = bc.sym(W + W*(Z*W)); del Z, W
        rH = bc.frob_rad(H); Hm = H.mid(); del H
        log(f"  parity {p}: d={d}  H6 epsJ = {bc._s(epsJ, 3)} (ev = {bc._s(ev, 3)});  ||W||_2 <= {bc._s(nW, 4)};  ||rad H||_F <= {bc._s(rH, 3)}  ({time.time()-T0:.0f}s)")
        summ['parity'][p] = {'d': d, 'epsJ': bc._s(epsJ, 4), 'radH_F': bc._s(rH, 4), 'normW_upper': bc._s(nW, 4)}
        ctx.prec = prec_c
        r0 = chol_inv(Hm)
        if r0 is None:
            log(f"  parity {p}: midpoint Cholesky of H_mid fails (H_mid not numerically PD) -> FAIL"); allok = False; ctx.prec = prec; continue
        Linv = r0[1]; del r0
        fr = arb(0)
        for i in range(d):
            for kk in range(i + 1): fr += Linv[i, kk].mid()**2
        del Linv
        lam_est = 1/fr                                       # ~ 1/trace(H_mid^-1) <= lambda_min(H_mid)
        log(f"  parity {p}: lambda_est(H_mid) = 1/||L^-1||_F^2 = {bc._s(lam_est, 4)}  ({time.time()-T0:.0f}s)")
        ok = False
        kexp = int(-float(lam_est.log().mid())/0.6931471805599453) + 1
        for k in (kexp, kexp + 2, kexp + 6):
            sigma = arb(2)**(-k)
            Hs = bc.mat_from(lambda i, kk: Hm[i, kk] - (sigma if i == kk else 0), d)
            r1 = chol_inv(Hs)
            if r1 is None: log(f"  parity {p}: Cholesky of H_mid - 2^-{k} I fails"); continue
            Lm = r1[0]; del r1
            Res = Hs - Lm*Lm.transpose(); del Lm, Hs
            res2 = arb(0)
            for i in range(d):
                for kk in range(d): res2 += arb(abs(Res[i, kk]).upper())**2
            del Res
            resF = arb(res2.sqrt().upper()); bound = resF + rH
            ok = bool(bound < sigma)
            log(f"  parity {p}: H7 H_mid - 2^-{k} I = L L^T + Res;  ||Res||_F + ||rad H||_F <= {bc._s(resF, 3)} + {bc._s(rH, 3)} = {bc._s(bound, 3)} {'<' if ok else 'NOT <'} 2^-{k} = {bc._s(sigma, 3)} -> {'PSD CERTIFIED' if ok else 'FAIL'}  ({time.time()-T0:.0f}s)")
            summ['parity'][p].update({'sigma': f'2^-{k}', 'resF': bc._s(resF, 4), 'chol_bound': bc._s(bound, 4), 'certified': ok})
            if ok:
                lamH = arb((sigma - bound).lower()); lamT = arb((lamH/nW**2).lower())
                log(f"  parity {p}: proved margins  lambda_min(H) >= {bc._s(lamH, 4)};  ||W||_2 <= {bc._s(nW, 4)};  lambda_min(W^-1 + Z) >= {bc._s(lamT, 4)} > 0")
                summ['parity'][p].update({'lambda_min_H_lower': bc._s(lamH, 4), 'lambda_min_T_lower': bc._s(lamT, 4)})
                break
        ctx.prec = prec
        del Hm
        allok = allok and ok
    summ['certified'] = allok
    if summary:
        with open(summary, 'w') as fh: json.dump(summ, fh, indent=1)
    if allok:
        log(f"RESULT mu={tag}: CERTIFIED  Q(f) >= {tau.str(6, radius=False)} ||f||^2 for every f supported in [-L, L], L = {L.str(10, radius=False)}")
        sys.exit(0)
    log(f"RESULT mu={tag}: not certified"); sys.exit(1)

if __name__ == "__main__":
    main()
