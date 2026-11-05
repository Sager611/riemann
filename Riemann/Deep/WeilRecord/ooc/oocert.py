"""oocert.py -- memory-lean (out-of-core) driver of the WeilRecord certificate chain (weilooc, 2026-11-05).

Mathematics: PROOF.md sections 1-5 and cert4.py, unchanged, with three reformulations (DERIVATION.md):
  * dyadic m = gamma - s - tau (gamma, s, tau dyadic), and W' = m W(X) = I + G1 X + (G1 X)^T + X^T G2 X with
    G1 = A0s - s I = m I - A,  G2 = (gamma-tau) s I - (gamma-tau+s) A0s + B0 = B2 - m A  (A0s = sym E*C_pE, B0 = E*C_p^2E),
    so that no division occurs; H' = m^2 H = m W' + W' Z W'.
  * every matrix is a fixed-point disk matrix (fxmm.py) with a rigorous 2-norm error bound carried here in arb;
  * PSD of H': block Cholesky H~' - sigma' I = L L^T + Res with Res computed EXACTLY; lambda_min(H'_true) >= sigma' -
    ||Res||_F - e_H', and lambda_min(W^-1 + Z) >= (sigma' - ||Res||_F - e_H') / (||W~'|| + e_W')^2  (m cancels).
Every stage is resumable; a process stops cleanly (exit 3) when its time budget is used, and the chain runner restarts it.
usage: oocert.py RUNDIR        (RUNDIR/config.json; see make_config.py)
"""
import sys, os, json, time, math, subprocess, hashlib, resource
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); WR = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(1, WR)
from flint import arb, ctx, fmpq, fmpz
from fxmm import Fx
import oocla, fxkern, bandplan
import bc
from oocla import View, Stack, FileMat, OutOfTime

BIN = os.environ.get('OOC_BIN', os.getcwd())

def rat(s):
    s = str(s).strip()
    if s.startswith('2^'): k = int(s[2:]); return fmpq(1, 2**(-k)) if k < 0 else fmpq(2**k)
    if '/' in s: a, b = s.split('/'); return fmpq(int(a), int(b))
    return fmpq(int(s))
def dyadic_bits(q):
    d = int(q.q)
    if d & (d - 1): return None
    return d.bit_length() - 1
def A(q): return arb(q)
def ub_store(x):
    """exact dyadic upper bound of an arb, as [man, exp] (json)."""
    u = x.upper(); m, e = u.mid().man_exp(); return [int(m), int(e)]
def ub_load(t): return arb(t[0])*arb(2)**t[1]
def ceil_scaled(x, k):
    """smallest integer >= upper(x) 2^k."""
    m, e = ub_store(x); e += k
    return m << e if e >= 0 else -((-m) >> (-e))
def fl(x): return float(x.upper()) if isinstance(x, arb) else float(x)
def S_(x, n=4): return x.str(n, radius=False)
def rss():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return f"{r/2**20:.0f} MB" if sys.platform == 'darwin' else f"{r/2**10:.0f} MB"

class Run:
    def __init__(self, rundir):
        self.dir = rundir
        with open(os.path.join(rundir, 'config.json')) as fh: self.cfg = json.load(fh)
        self.stp = os.path.join(rundir, 'state.json')
        self.st = json.load(open(self.stp)) if os.path.exists(self.stp) else {'done': [], 'v': {}}
        c = self.cfg
        self.F, self.T, self.N = int(c['F']), int(c['T']), int(c['N'])
        self.fx = Fx(self.F, self.T)
        self.b = int(c.get('block', 400))
        self.t0 = time.time(); self.deadline = self.t0 + float(c.get('budget', 2000))
        self.logf = open(os.path.join(rundir, 'log.txt'), 'a')
        ctx.prec = self.F + 64
        self.Sx = bc.Setup(c['tag'])
        self.gamma, self.s, self.tau = rat(c['gamma']), rat(c['s']), rat(c['tau'])
        self.m = self.gamma - self.s - self.tau
        for nm, q in (('gamma', self.gamma), ('s', self.s), ('tau', self.tau)):
            k = dyadic_bits(q)
            if k is None or k > self.F - 2: raise SystemExit(f"{nm} must be dyadic with <= F-2 fractional bits")
        self.rho_out = arb(2)**(-self.F - 1)       # rounding of one stored entry
    def log(self, *a):
        msg = ' '.join(str(x) for x in a)
        line = f"[{time.strftime('%H:%M:%S')} +{time.time()-self.t0:.0f}s {rss()}] {msg}"
        print(line, flush=True); self.logf.write(line + '\n'); self.logf.flush()
    def save(self):
        tmp = self.stp + '.tmp'
        with open(tmp, 'w') as fh: json.dump(self.st, fh, indent=1)
        os.replace(tmp, self.stp)
    def v(self, k): return ub_load(self.st['v'][k])
    def setv(self, k, x): self.st['v'][k] = ub_store(x)
    def path(self, name): return os.path.join(self.dir, name)
    def stage(self, name, fn):
        if name in self.st['done']: return
        self.log(f"stage {name}")
        ts = time.time()
        fn()
        self.st['done'].append(name); self.save()
        self.log(f"stage {name} done ({time.time() - ts:.0f}s)")
        if time.time() - ts > 30 or oocla.rss_now_mb() > 300:
            raise OutOfTime(f"stage boundary after {name}")      # fresh process for the next stage (memory)
    def d(self, p): return (self.N + 1 - p)//2
    def int2F(self, q, k=1):
        """q * 2^(kF) as an exact integer (q dyadic)."""
        x = q*fmpq(2**(k*self.F))
        if int(x.q) != 1: raise ValueError("not representable at this scale")
        return int(x.p)

    # ---------------- stage: constants (H1-H5), pole vector ----------------
    def consts(self):
        c = self.cfg; S = self.Sx; L = S.L_()
        ctx.prec = 128
        s_cert, w_even, Mc, _ = bc.schur(S, c['wfile'])
        ctx.prec = self.F + 64
        s = A(self.s); gamma = A(self.gamma); tau = A(self.tau); m = A(self.m)
        if not (s_cert <= s): raise SystemExit(f"H1 FAILED: Schur bound {s_cert} > s = {s}")
        if not (m > 0): raise SystemExit("H2 FAILED: m_t <= 0")
        Omega = A(rat(c['Omega'])); b = gamma + s; cc = Omega*L
        C = bc.a_fun(Omega) - gamma
        if not (C > 0): raise SystemExit("H3 FAILED")
        E0 = (cc + 1)/(arb.pi()*cc**2); E1 = (3*cc**2 + 6*cc + 2)/(arb.pi()*cc**3)
        u0, u1 = arb(int(c['u0'])), arb(int(c['u1']))
        if not (u0*E0 <= C and u1*E1 <= C): raise SystemExit("H4 FAILED")
        N = self.N
        K0 = (Omega/arb.pi())*max(abs(bc.a_fun(arb(0)) - gamma), abs(bc.a_fun(Omega) - gamma))
        Hp = (2*L).sqrt()*(L/2).cosh(); qN = (2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(N).fac()
        best = None
        for rho in (1.4, 1.6, 1.8, 2, 2.2, 2.5, 3, 3.5, 4, 5, 6):
            rho = arb(rho)
            rN = ((8*L*rho**2/(rho - 1)**2)*(cc*(rho - 1/rho)).exp()*rho**(-2*N)).sqrt()
            e = K0*(2*L).sqrt()*rN + (Omega/arb.pi()).sqrt()*(u0*E0.sqrt() + u1*E1.sqrt())*rN + 4*Hp*qN
            n = K0*rN**2 + 4*qN**2
            h = (K0 + (Omega/arb.pi())*(u0 + u1))*rN**2 + 4*qN**2
            if not (n < m): continue
            alpha = e + 2*e**2/(m - n); beta = (e + n)/m**2 + 2*h**2/(m**2*(m - n))
            tot = alpha + beta*b**2
            if best is None or tot < best[0]: best = (tot, rho, rN, e, n, h, alpha, beta)
        if best is None: raise SystemExit("H5 FAILED: no admissible rho")
        _, rho, rN, e, n, h, alpha, beta = best
        corr = alpha + beta*(b - tau)**2
        self.log(f"H1 ||C_p|| <= {S_(s_cert, 10)} <= s = {S_(s, 10)} (all {Mc} cells, weight even: {w_even});  H2 m_t = {S_(m, 12)};  "
                 f"H3 C = a(Omega) - gamma = {S_(C, 8)};  H4 u0 E0 = {S_(u0*E0, 8)}, u1 E1 = {S_(u1*E1, 8)}")
        self.log(f"H5 rho={float(rho.mid())} r_N={S_(rN, 3)} e={S_(e, 3)} n={S_(n, 3)} h={S_(h, 3)} alpha={S_(alpha, 3)} beta={S_(beta, 3)} corr={S_(corr, 3)};  N/c={N/float(cc.mid()):.3f}")
        for k, x in (('s_cert', s_cert), ('C', C), ('E0', E0), ('E1', E1), ('rN', rN), ('e', e), ('n', n), ('h', h), ('alpha', alpha), ('beta', beta), ('corr', corr)):
            self.setv(k, x)
        self.st['H'] = {'H1': S_(s_cert, 12), 'cells': Mc, 'weight_even': w_even, 'H3_C': S_(C, 10), 'H4': [S_(u0*E0, 10), S_(u1*E1, 10)],
                        'H5': {'rho': float(rho.mid()), 'rN': S_(rN, 4), 'e': S_(e, 4), 'n': S_(n, 4), 'h': S_(h, 4), 'alpha': S_(alpha, 4), 'beta': S_(beta, 4), 'corr': S_(corr, 4)}}
        # pole vector Av_j = (2(2j+1)L)^(1/2) (pi/L)^(1/2) I_{j+1/2}(L/2), stored per parity (fixed point) + error
        ctx.prec = self.F + 64 + 64
        fx = self.fx
        for p in (0, 1):
            d = self.d(p); D = np.zeros((d, 1, fx.ND), dtype=np.uint16); err = arb(0)
            for i in range(d):
                j = 2*i + p
                v = (2*(2*j + 1)*L).sqrt()*(arb.pi()/L).sqrt()*(L/2).bessel_i(arb(j) + arb(1)/2)
                assert v.is_finite()
                mm, ee = v.mid().man_exp(); mm, ee = int(mm), int(ee) + self.F
                X = mm << ee if ee >= 0 else (mm + (1 << (-ee - 1))) >> (-ee)
                D[i, 0] = fx.int_digits(X)
                err = max(err, arb(v.rad()) + arb(2)**(-self.F - 1))
            np.save(self.path(f'Av_{p}.npy'), D)
            self.setv(f'Av_err_{p}', err)
        ctx.prec = self.F + 64
        # dyadic shift checks
        self.setv('m', A(self.m))

    # ---------------- stage: Legendre values (gvals) ----------------
    def gvals(self):
        c = self.cfg; gdir = self.path('gv'); os.makedirs(gdir, exist_ok=True)
        cmd = [os.path.join(BIN, 'gvals'), c['tag'], str(self.N), str(self.F), str(self.T), str(int(c.get('FV', self.F + 64))), gdir]
        done = self.st.setdefault('gv_done', [])
        npieces = len(self.Sx.ns) + 1          # breakpoints 0, 1 and |1 - log n / L| (all in [0, 1] since n < e^{2L})
        self.st['gv_npieces'] = npieces
        for ip in range(npieces):
            if ip in done: continue
            if time.time() > self.deadline - 300: raise OutOfTime('gvals')     # a piece must finish before the alarm
            out = subprocess.run(cmd + [str(ip)], capture_output=True, text=True)
            if out.returncode != 0: raise SystemExit(f"gvals piece {ip} failed: {out.stdout[-2000:]} {out.stderr[-2000:]}")
            hdr = [l for l in out.stdout.split('\n') if l.startswith('gvals:')]
            if not hdr or int(hdr[0].split('pieces=')[1].split()[0]) != npieces: raise SystemExit("gvals: piece count mismatch")
            done.append(ip); self.save()
            self.log(out.stdout.strip().split('\n')[-2])
        self.gv_budget()

    def gv_metas(self):
        ms = []
        for ip in range(self.st['gv_npieces']):
            with open(self.path(f'gv/gv_meta_pc{ip}.json')) as fh: ms.append(json.load(fh))
        return ms

    def gv_budget(self):
        """Gram error budget (DERIVATION.md section 2): 2-norm bounds E_A, E_B of A0s, B0 for each parity."""
        ms = [m for m in self.gv_metas() if m.get('nodes', 0) > 0]
        N, F = self.N, self.F; FV = ms[0]['FV']
        sm = lambda key: sum((arb(m[key]) for m in ms), arb(0))
        mx = lambda key: max((arb(m[key]) for m in ms), key=lambda x: float(x.mid()))
        eA, eB = sm('eA'), sm('eB'); W1, SW = sm('sumw'), sm('sumsqrtw')
        nq = sum(m['nodes'] for m in ms)
        du, dy, dsig, cferr = mx('du'), mx('dy'), mx('dsig'), mx('cferr')
        csum = mx('csum'); nact = max(m['active_shifts'] for m in ms)
        two = arb(2); Nn = arb(N)
        rec = arb('2.1')*(Nn + 1)**2*two**(-FV)                 # rotation recurrence error, j < N
        etaP = rec + Nn*(Nn + 1)/2*du
        etaC = csum*(rec + Nn*(Nn + 1)/2*dy) + nact*cferr + two**(-FV - 1)
        rhoP = two**(-F - 1) + dsig*(1 + etaP)
        rhoC = two**(-F - 1) + dsig*(csum + etaC)
        for p in (0, 1):
            d = self.d(p)
            S1 = arb(sum(2*(2*i + p) + 1 for i in range(d)))
            def gram_err(etaX, MX, rhoX, etaY, MY, rhoY):
                return (S1*W1*(etaX*(MY + etaY) + MX*etaY) + (arb(d)*S1).sqrt()*SW*(etaX*rhoY + MX*rhoY + rhoX*(MY + etaY))
                        + arb(d)*arb(nq)*rhoX*rhoY)
            EA = S1*eA + gram_err(etaP, arb(1), rhoP, etaC, csum, rhoC) + arb(d)*self.rho_out
            EB = S1*eB + gram_err(etaC, csum, rhoC, etaC, csum, rhoC) + arb(d)*self.rho_out
            self.setv(f'EA_{p}', EA); self.setv(f'EB_{p}', EB)
            self.log(f"gram budget parity {p}: nodes {nq}, quad {S_(S1*eA, 3)} / {S_(S1*eB, 3)}, value eta_P {S_(etaP, 3)} eta_C {S_(etaC, 3)}; E_A {S_(EA, 3)}, E_B {S_(EB, 3)}")
        self.st['gv'] = {'nodes': nq, 'pieces': len(ms), 'eA': S_(eA, 4), 'eB': S_(eB, 4), 'etaP': S_(etaP, 3), 'etaC': S_(etaC, 3)}

    # ---------------- stage: band values ----------------
    def band(self):
        c = self.cfg; bdir = self.path('bv'); os.makedirs(bdir, exist_ok=True)
        if 'bplan' not in self.st:
            ctx.prec = self.F
            nodes, pl, eK, eB, parts, supg = bandplan.plan(self.Sx, self.N, A(rat(c['Omega'])), A(self.gamma), self.F,
                                                           target_bits=int(c.get('band_target_bits', self.F - 40)))
            ctx.prec = self.F
            split, wmax = bandplan.write_nodes(os.path.join(bdir, 'band_nodes.txt'), nodes, A(self.gamma))
            ctx.prec = self.F + 64
            self.setv('errK', eK); self.setv('errB', eB); self.setv('supg', supg)
            self.st['bplan'] = {'nodes': len(nodes), 'split': split, 'panels': len(pl), 'errK': S_(eK, 4), 'errB': S_(eB, 4), 'wmax': wmax,
                                'parts': {k: S_(v, 3) for k, v in parts.items()}}
            self.save(); self.log(f"band plan: {len(nodes)} nodes ({len(pl)} panels), split {split}, errK {S_(eK, 3)}, errB {S_(eB, 3)}")
        if time.time() > self.deadline - float(self.cfg.get('bvals_time', 900)): raise OutOfTime('band')
        L = self.Sx.L_(); xmax = float((L*A(rat(c['Omega']))).upper())
        bprec = int(c.get('bprec', self.F + int(xmax) + 128))
        cmd = [os.path.join(BIN, 'bvals'), c['tag'], str(self.N), str(self.F), str(self.T), str(bprec), os.path.join(bdir, 'band_nodes.txt'), bdir]
        out = subprocess.run(cmd, capture_output=True, text=True)
        if out.returncode != 0: raise SystemExit(f"bvals failed: {out.stdout[-2000:]} {out.stderr[-2000:]}")
        self.log(out.stdout.strip().split('\n')[-1])
        with open(os.path.join(bdir, 'bv_meta.json')) as fh: bm = json.load(fh)
        self.st['bv'] = bm

    # ---------------- per-parity stages ----------------
    def gv_srcs(self, p):
        d = self.d(p); fx = self.fx; Ps, Cs = [], []
        for m in self.gv_metas():
            if m.get('nodes', 0) == 0: continue
            ip, nq = m['piece'], m['nodes']
            Ps.append(View(FileMat(self.path(f"gv/gv_p{p}_pc{ip}_P.bin"), shape=(nq, d, fx.ND))))
            Cs.append(View(FileMat(self.path(f"gv/gv_p{p}_pc{ip}_C.bin"), shape=(nq, d, fx.ND))))
        return Stack([s.T for s in Ps], 1), Stack(Cs, 0), Stack([s.T for s in Cs], 1)

    def gram(self, p):
        d = self.d(p); fx = self.fx; b = self.b
        PT, Cst, CT = self.gv_srcs(p)
        oocla.run_blocks(fx, self.path(f'A0_{p}.npy'), d, d, b, lambda r0, r1, c0, c1: [('prod', PT, Cst, 1)], deadline=self.deadline, log=self.log, tag=f'A0 p{p}')
        oocla.run_blocks(fx, self.path(f'B0_{p}.npy'), d, d, b, lambda r0, r1, c0, c1: [('prod', CT, Cst, 1)], sym=True, deadline=self.deadline, log=self.log, tag=f'B0 p{p}')

    def Gmats(self, p):
        d = self.d(p); fx = self.fx; b = self.b; F = self.F
        A0 = oocla.view(self.path(f'A0_{p}.npy')); B0 = oocla.view(self.path(f'B0_{p}.npy'))
        half = 1 << (F - 1)
        sI = -self.int2F(self.s, 2)
        oocla.run_blocks(fx, self.path(f'G1_{p}.npy'), d, d, b, lambda *a: [('mat', A0, half), ('mat', A0.T, half), ('diag', sI)],
                         sym=True, deadline=self.deadline, log=self.log, tag=f'G1 p{p}')
        gt = self.gamma - self.tau
        cA = -self.int2F((gt + self.s)/2); cD = self.int2F(gt*self.s, 2)
        oocla.run_blocks(fx, self.path(f'G2_{p}.npy'), d, d, b, lambda *a: [('mat', A0, cA), ('mat', A0.T, cA), ('mat', B0, 1 << F), ('diag', cD)],
                         sym=True, deadline=self.deadline, log=self.log, tag=f'G2 p{p}')
        dd = arb(d)
        EG1 = self.v(f'EA_{p}') + dd*self.rho_out
        EG2 = (A(gt) + A(self.s))*self.v(f'EA_{p}') + self.v(f'EB_{p}') + dd*self.rho_out
        self.setv(f'EG1_{p}', EG1); self.setv(f'EG2_{p}', EG2)
        self.log(f"parity {p}: ||G1 - G1~|| <= {S_(EG1, 3)}, ||G2 - G2~|| <= {S_(EG2, 3)}")

    def solveX(self, p):
        """X ~ G2^-1 (-G1) by mixed-precision iterative refinement (double inverse, exact residuals).  X is DATA."""
        d = self.d(p); fx = self.fx; b = self.b; F = self.F
        key = f'X_{p}'
        stx = self.st.setdefault('Xref', {}).setdefault(str(p), {'it': 0, 'hist': []})
        G1a = oocla.fm(self.path(f'G1_{p}.npy')); G2a = oocla.fm(self.path(f'G2_{p}.npy'))
        G2d = oocla.to_float_mat(fx, G2a); inv = np.linalg.inv(G2d)
        if stx['it'] == 0:
            ev = np.linalg.eigvalsh(G2d); self.log(f"parity {p}: G2 eigenvalues (double) [{ev[0]:.3e}, {ev[-1]:.3e}], cond {ev[-1]/ev[0]:.2e}")
            stx['G2_eig'] = [float(ev[0]), float(ev[-1])]
        del G2d
        maxit = int(self.cfg.get('X_maxit', 14)); tol = 2.0**(-int(self.cfg.get('X_tolbits', 300)))
        while True:
            it = stx['it']
            Xp = self.path(f'X_{p}_{it}.npy')
            if it == 0 and not os.path.exists(Xp):
                G1d = oocla.to_float_mat(fx, G1a)
                X0 = -(inv @ G1d); del G1d
                oocla.write_float(fx, Xp, X0); del X0
            if it >= maxit: break
            oocla.check_budget(self.deadline, 'X')
            Xa = oocla.view(Xp); G1 = View(G1a); G2 = View(G2a)
            Rp = self.path(f'R_{p}_{it}.npy')
            oocla.run_blocks(fx, Rp, d, d, b, lambda *a: [('mat', G1, -(1 << F)), ('prod', G2, Xa, -1)], deadline=self.deadline, log=self.log, tag=f'R{it} p{p}')
            Rd = oocla.to_float_mat(fx, oocla.fm(Rp)); rmax = float(np.abs(Rd).max())
            delta = inv @ Rd; del Rd
            dmax = float(np.abs(delta).max())
            stx['hist'].append([rmax, dmax]); self.log(f"parity {p}: refinement {it}: max|R| {rmax:.3e}, max|dX| {dmax:.3e}")
            if dmax < tol or (it > 2 and stx['hist'][-1][1] > 0.5*stx['hist'][-2][1]):
                os.remove(Rp); break
            Dp = self.path(f'dX_{p}.npy'); oocla.write_float(fx, Dp, delta); del delta
            Da = oocla.view(Dp)
            Xn = self.path(f'X_{p}_{it + 1}.npy')
            oocla.run_blocks(fx, Xn, d, d, b, lambda *a: [('mat', Xa, 1 << F), ('mat', Da, 1 << F)], deadline=None, log=self.log, tag=f'X{it+1} p{p}')
            stx['it'] = it + 1; self.save()
            for f in (Dp, Rp, Xp): os.remove(f)
        os.replace(self.path(f'X_{p}_{stx["it"]}.npy'), self.path(f'X_{p}.npy'))
        nX = arb(oocla.frob_bound(fx, oocla.fm(self.path(f'X_{p}.npy')))); self.setv(f'normX_{p}', nX)
        self.log(f"parity {p}: X done ({stx['it']} refinements), ||X||_F <= {S_(nX, 4)}")

    def Wmat(self, p):
        d = self.d(p); fx = self.fx; b = self.b; F = self.F
        Xa = oocla.view(self.path(f'X_{p}.npy'))
        G1 = oocla.view(self.path(f'G1_{p}.npy')); G2 = oocla.view(self.path(f'G2_{p}.npy'))
        oocla.run_blocks(fx, self.path(f'Y1_{p}.npy'), d, d, b, lambda *a: [('prod', G1, Xa, 1)], deadline=self.deadline, log=self.log, tag=f'Y1 p{p}')
        oocla.run_blocks(fx, self.path(f'Y2_{p}.npy'), d, d, b, lambda *a: [('prod', G2, Xa, 1)], deadline=self.deadline, log=self.log, tag=f'Y2 p{p}')
        Y1 = oocla.view(self.path(f'Y1_{p}.npy')); Y2 = oocla.view(self.path(f'Y2_{p}.npy'))
        oocla.run_blocks(fx, self.path(f'W_{p}.npy'), d, d, b,
                         lambda *a: [('diag', 1 << (2*F)), ('mat', Y1, 1 << F), ('mat', Y1.T, 1 << F), ('prod', Xa.T, Y2, 1)],
                         sym=True, deadline=self.deadline, log=self.log, tag=f'W p{p}')
        Wa = oocla.fm(self.path(f'W_{p}.npy'))
        nX = self.v(f'normX_{p}'); dd = arb(d); r = self.rho_out
        eW = 2*(nX*self.v(f'EG1_{p}') + dd*r) + dd*r + arb(2).sqrt()*(nX*(nX*self.v(f'EG2_{p}') + dd*r) + dd*r)
        nWinf = arb(oocla.inf_bound(fx, Wa)); nWF = arb(oocla.frob_bound(fx, Wa))
        nW = nWinf if nWinf < nWF else nWF
        self.setv(f'eW_{p}', eW); self.setv(f'normW_{p}', nW)
        self.log(f"parity {p}: W' = m W(X) assembled: ||W'~||_2 <= {S_(nW, 4)} (inf {S_(nWinf, 4)}, F {S_(nWF, 4)}), ||W'_true - W'~||_2 <= {S_(eW, 3)}")
        for f in ('Y1', 'Y2'): os.remove(self.path(f'{f}_{p}.npy'))

    def Zmat(self, p):
        d = self.d(p); fx = self.fx; b = self.b; F = self.F; c = self.cfg
        bm = self.st['bv']; nq = bm['nodes']; split = self.st['bplan']['split']
        Bv = FileMat(self.path(f'bv/bv_p{p}.bin'), shape=(nq, d, fx.ND))
        Bneg, Bpos = View(Bv, 0, 0, split, d), View(Bv, split, 0, nq - split, d)
        # rank-2 term: V = [Av, Eh], Vs = [sg Av, u Eh];  Eh = e_0 - b0
        Av = np.load(self.path(f'Av_{p}.npy'))
        lines = open(self.path(f'bv/bv_b0_p{p}.txt')).read().split('\n')
        Eh = np.zeros((d, 1, fx.ND), dtype=np.uint16); b0rad = arb(0); b0r = arb(0)
        for i in range(d):
            mm, ee, rr = lines[i].split(); mm, ee = int(mm), int(ee) + F
            X = mm << ee if ee >= 0 else (mm + (1 << (-ee - 1))) >> (-ee)
            X = (1 << F if i == 0 else 0) - X
            Eh[i, 0] = fx.int_digits(X); b0rad = max(b0rad, arb(rr))
        sg = 2 if p == 0 else -2; u = int(c['u0'] if p == 0 else c['u1'])
        V = np.concatenate([Av, Eh], axis=1)
        Vs = np.zeros_like(V)
        for i in range(d):
            Vs[i, 0] = fx.int_digits(sg*fx.digits_int(V[i, 0])); Vs[i, 1] = fx.int_digits(u*fx.digits_int(V[i, 1]))
        Vsrc, Vss = View(V), View(Vs)
        # epsJ: ||J_true - (K~ + rank2~)||_2 (before the final rounding of Z)
        dd = arb(d); L = self.Sx.L_()
        errK, errB, supg = self.v('errK'), self.v('errB'), self.v('supg')
        maxerr = arb(bm[f'maxerr{p}']); wrad, grad = arb(bm['wrad']), arb(bm['grad'])
        Omega = A(rat(c['Omega'])); c2Lpi = 2*L/arb.pi()
        S1 = arb(sum(2*(2*i + p) + 1 for i in range(d)))
        sigsum = (c2Lpi*arb(nq)*Omega*(supg + grad)).sqrt()*(1 + arb(2)**-30)    # >= sum_q sqrt(c2Lpi w_q |g_q|)
        E_valK = maxerr*(2*sigsum*(dd*S1).sqrt() + dd*arb(nq)*maxerr)
        wmax = arb(self.st['bplan']['wmax'])
        E_wg = c2Lpi*S1*arb(nq)*(wrad*(supg + grad) + wmax*grad)
        Kerr = dd*errK + E_valK + E_wg
        dAv = self.v(f'Av_err_{p}')*dd.sqrt()
        Avnorm = arb(math.sqrt(sum(fx.to_float(V[i, 0])**2 for i in range(d))))*(1 + arb(2)**-30) + dAv
        Eh_err = dd.sqrt()*(errB + b0rad + self.rho_out)
        epsJ = Kerr + abs(sg)*(2*Avnorm*dAv + dAv**2) + u*(2*Eh_err + Eh_err**2)
        shift = epsJ + self.v('corr') + 2*dd*self.rho_out
        sh_int = ceil_scaled(shift, 2*F) + 1
        self.setv(f'epsJ_{p}', epsJ)
        self.log(f"parity {p}: epsJ = d errK {S_(dd*errK, 3)} + values {S_(E_valK, 3)} + w,g {S_(E_wg, 3)} + poles {S_(abs(sg)*(2*Avnorm*dAv + dAv**2), 3)} + U {S_(u*(2*Eh_err + Eh_err**2), 3)} = {S_(epsJ, 3)};  shift = epsJ + corr {S_(self.v('corr'), 3)} + rounding")
        oocla.run_blocks(fx, self.path(f'Z_{p}.npy'), d, d, b,
                         lambda *a: [('prod', Bpos.T, Bpos, 1), ('prod', Bneg.T, Bneg, -1), ('prod', Vss, Vsrc.T, 1), ('diag', -sh_int)],
                         sym=True, deadline=self.deadline, log=self.log, tag=f'Z p{p}')
        Za = oocla.fm(self.path(f'Z_{p}.npy'))
        nZ = arb(min(oocla.inf_bound(fx, Za), oocla.frob_bound(fx, Za))); self.setv(f'normZ_{p}', nZ)
        self.log(f"parity {p}: Z = J~ - shift I stored exactly (it is the Z of the proof); ||Z||_2 <= {S_(nZ, 4)}")

    def Hmat(self, p):
        d = self.d(p); fx = self.fx; b = self.b; F = self.F
        Wa = oocla.view(self.path(f'W_{p}.npy')); Za = oocla.view(self.path(f'Z_{p}.npy'))
        oocla.run_blocks(fx, self.path(f'ZW_{p}.npy'), d, d, b, lambda *a: [('prod', Za, Wa, 1)], deadline=self.deadline, log=self.log, tag=f'ZW p{p}')
        ZW = oocla.view(self.path(f'ZW_{p}.npy'))
        mI = self.int2F(self.m)
        oocla.run_blocks(fx, self.path(f'H_{p}.npy'), d, d, b, lambda *a: [('mat', Wa, mI), ('prod', Wa, ZW, 1)],
                         sym=True, deadline=self.deadline, log=self.log, tag=f'H p{p}')
        eW, nW, nZ = self.v(f'eW_{p}'), self.v(f'normW_{p}'), self.v(f'normZ_{p}')
        dd = arb(d); r = self.rho_out; m = A(self.m)
        eH = m*eW + arb(2).sqrt()*(eW*nZ*(2*nW + eW) + nW*dd*r + dd*r)
        self.setv(f'eH_{p}', eH)
        self.log(f"parity {p}: H' = m W' + W' Z W' assembled; ||H'_true - H'~||_2 <= {S_(eH, 3)}")
        os.remove(self.path(f'ZW_{p}.npy'))

    # ---------------- Cholesky, inverse iteration, residual ----------------
    def chol(self, p, sig_int, name):
        """block Cholesky of H~' - (sig_int/2^2F) I into file name; progress in state.  Returns True / False (not PD)."""
        d = self.d(p); fx = self.fx; b = self.b; F = self.F; prec = F + 64
        Ha = oocla.view(self.path(f'H_{p}.npy'))
        Lp = self.path(name)
        stc = self.st.setdefault('chol', {}).setdefault(name, {'done': [], 'sig': sig_int})
        if not os.path.exists(Lp):
            FileMat.create(Lp + '.tmp.npy', d, d, fx.ND, fill=fx.int_digits(0)).close(); os.replace(Lp + '.tmp.npy', Lp)
        Lf = FileMat(Lp, 'w')
        B = oocla.blocks(d, b); nb = len(B)
        for k in range(nb):
            (k0, k1) = B[k]
            Lk = View(Lf, 0, 0, d, k0) if k0 > 0 else None
            if ['d', k] not in stc['done']:
                try: oocla.check_budget(self.deadline, 'chol')
                except OutOfTime: Lf.close(); raise
                tb = time.time()
                terms = [('mat', Ha, 1 << F), ('diag', -sig_int)] + ([('prod', Lk, Lk.T, -1)] if Lk else [])
                S = oocla.combine_block(fx, k0, k1, k0, k1, terms)
                Lkk, ret = fxkern.chol(fx, S, prec)
                if ret != 0:
                    self.log(f"  {name}: pivot failure in diagonal block {k} (ret {ret}): H~' - sigma' I is not numerically PD")
                    stc['fail'] = [k, ret]; self.save(); Lf.close(); return False
                Lf.write(k0, k0, Lkk); Lf.flush()
                stc['done'].append(['d', k]); self.save()
                self.log(f"    {name} diag block {k}/{nb}: {time.time()-tb:.1f}s")
            Lkk = Lf.read(k0, k1, k0, k1)
            for i in range(k + 1, nb):
                if ['o', i, k] in stc['done']: continue
                try: oocla.check_budget(self.deadline, 'chol')
                except OutOfTime: Lf.close(); raise
                tb = time.time()
                (i0, i1) = B[i]
                terms = [('mat', Ha, 1 << F)] + ([('prod', Lk, Lk.T, -1)] if Lk else [])
                S = oocla.combine_block(fx, i0, i1, k0, k1, terms)
                Lf.write(i0, k0, fxkern.trsm(fx, S, Lkk, prec)); Lf.flush()
                stc['done'].append(['o', i, k]); self.save()
                self.log(f"    {name} block ({i},{k}): {time.time()-tb:.1f}s")
        Lf.close()
        return True

    def tri_solve(self, Lf, v, trans):
        """high-precision solve L y = v (trans: L^T y = v) for one vector (digits d x 1), block substitution."""
        fx = self.fx; F = self.F; d = Lf.rows; B = oocla.blocks(d, self.b); prec = F + 64
        y = np.zeros_like(v)
        order = range(len(B)) if not trans else range(len(B) - 1, -1, -1)
        for I in order:
            (r0, r1) = B[I]
            if not trans and r0 > 0:
                rhs = oocla.combine_block(fx, r0, r1, 0, 1, [('mat', View(v), 1 << F), ('prod', View(Lf, 0, 0, d, r0), View(y, 0, 0, r0, 1), -1)])
            elif trans and r1 < d:
                rhs = oocla.combine_block(fx, r0, r1, 0, 1, [('mat', View(v), 1 << F), ('prod', View(Lf, r1, 0, d - r1, d).T, View(y, r1, 0, d - r1, 1), -1)])
            else:
                rhs = np.ascontiguousarray(v[r0:r1])
            y[r0:r1] = fxkern.trsm_left(fx, Lf.read(r0, r1, r0, r1), rhs, prec, trans=trans)
        return y

    def lam_est(self, p, name, iters=10):
        """Diagnostic only: inverse power iteration with the factor L (L L^T ~ H~'), estimate of lambda_min(H~').
        Vectors are rescaled by powers of 2 before each solve so that the fixed-point range is never exceeded."""
        fx = self.fx; d = self.d(p)
        Lf = FileMat(self.path(name))
        piv = min(float(np.min(np.abs(fx.to_float(np.ascontiguousarray(np.diagonal(Lf.read(r0, r1, r0, r1), axis1=0, axis2=1).T)))))
                  for r0, r1 in oocla.blocks(d, self.b))
        self.log(f"    smallest Cholesky pivot^2 {piv**2:.3e} (>= lambda_min)")
        rng = np.random.default_rng(7); x = rng.standard_normal(d); x /= np.linalg.norm(x)
        sc = max(piv*2.0**-8, 2.0**-400); lam = None
        def solve(vec, trans, s):
            while True:
                try:
                    return fx.to_float(self.tri_solve(Lf, fx.from_float((vec*s).reshape(d, 1)), trans)).reshape(d), s
                except OverflowError:
                    s *= 2.0**-24
                    if s < 2.0**-420: raise
        for it in range(iters):
            oocla.check_budget(self.deadline, 'lam')
            y, s1 = solve(x, False, sc)
            ny = float(np.abs(y).max()); yn = y/ny
            z, s2 = solve(yn, True, sc)
            nz = float(np.linalg.norm(z))       # z = s1 s2 (L L^T)^-1 x / ny
            new = s1*s2/(ny*nz)
            x = z/nz
            self.log(f"    inverse iteration {it}: lambda_min(H~') ~ {new:.6e}")
            if lam is not None and abs(new - lam) < 1e-5*abs(new): lam = new; break
            lam = new
            sc = math.sqrt(lam)*2.0**-4
        Lf.close()
        return lam

    def residual(self, p, name, sig_int):
        """||H~' - sigma' I - L L^T||_F (rigorous upper bound; the residual is computed exactly by CRT)."""
        d = self.d(p); fx = self.fx; b = self.b; F = self.F
        Ha = oocla.view(self.path(f'H_{p}.npy')); Lf = FileMat(self.path(name))
        B = oocla.blocks(d, b)
        str_ = self.st.setdefault('res', {}).setdefault(name, {'done': [], 'sumsq': 0.0})
        for I in range(len(B)):
            for J in range(I + 1):
                if [I, J] in str_['done']: continue
                oocla.check_budget(self.deadline, 'res')
                (r0, r1), (c0, c1) = B[I], B[J]
                Lr = View(Lf, 0, 0, d, c1)
                acc = oocla.combine_block(fx, r0, r1, c0, c1, [('mat', Ha, 1 << F), ('diag', -sig_int), ('prod', Lr, Lr.T, -1)], raw=True)
                vals, rel = fx.crt_float(acc.reshape(fx.K, -1), 2*F)
                ss = float(np.sum(vals*vals))*(1 + rel)**2*(1 + 1e-9)
                str_['sumsq'] += ss*(2 if I != J else 1); str_['done'].append([I, J]); self.save()
        return math.sqrt(str_['sumsq'])*(1 + 1e-9)

    def certify(self, p):
        F = self.F; fx = self.fx
        stp = self.st.setdefault('cert', {}).setdefault(str(p), {})
        if 'lam' not in stp:
            ok = self.chol(p, 0, f'L0_{p}.npy')
            if not ok:
                stp['result'] = 'FAIL: H~\' not PD'; self.save(); return
            if 'lam_pending' not in stp:
                stp['lam_pending'] = True; self.save(); raise OutOfTime('after chol L0')
            lam = self.lam_est(p, f'L0_{p}.npy')
            stp['lam'] = lam; self.save()
            os.remove(self.path(f'L0_{p}.npy'))
            raise OutOfTime('after lambda estimate')
        lam = stp['lam']
        k = int(math.floor(-math.log2(lam/2))) + 1 if 'k' not in stp else stp['k']
        while True:
            stp['k'] = k; self.save()
            sig_int = 1 << (2*F - k)
            name = f'L_{p}_{k}.npy'
            ok = self.chol(p, sig_int, name)
            if ok: break
            os.remove(self.path(name)); k += 2
            if k > 2*F - 8: stp['result'] = 'FAIL'; self.save(); return
        if not stp.get('chol_ok'):
            stp['chol_ok'] = True; self.save(); raise OutOfTime('after chol L')
        resF = self.residual(p, name, sig_int)
        if not self.cfg.get('keep', False):
            for f in (name, f'H_{p}.npy'):
                if os.path.exists(self.path(f)): os.remove(self.path(f))
        sigma = arb(2)**(-k); eH = self.v(f'eH_{p}')
        lamH = sigma - arb(resF) - eH
        nW = self.v(f'normW_{p}'); eW = self.v(f'eW_{p}')
        lamT = lamH/(nW + eW)**2
        ok = bool(lamH > 0)
        stp.update({'sigma_prime': f'2^-{k}', 'resF': resF, 'eH': S_(eH, 4), 'lamH_prime': S_(lamH, 4), 'normW': S_(nW, 4), 'eW': S_(eW, 3),
                    'lambda_min_T_lower': S_(lamT, 4) if ok else None, 'certified': ok})
        self.save()
        self.log(f"parity {p}: H7 H~' - 2^-{k} I = L L^T + Res;  ||Res||_F <= {resF:.3e}, e_H' <= {S_(eH, 3)}:  lambda_min(H'_true) >= {S_(lamH, 4)} -> "
                 + (f"CERTIFIED;  lambda_min(W^-1 + Z) >= {S_(lamT, 4)} > 0" if ok else "FAIL"))

    def run(self):
        c = self.cfg
        self.stage('consts', self.consts)
        self.stage('gvals', self.gvals)
        self.stage('band', self.band)
        for p in c.get('parities', [0, 1]):
            self.stage(f'gram{p}', lambda: self.gram(p))
            self.stage(f'G{p}', lambda: self.Gmats(p))
            if f'G{p}' in self.st['done'] and not c.get('keep', False):
                for f in ('A0', 'B0'):
                    if os.path.exists(self.path(f'{f}_{p}.npy')): os.remove(self.path(f'{f}_{p}.npy'))
            self.stage(f'X{p}', lambda: self.solveX(p))
            self.stage(f'W{p}', lambda: self.Wmat(p))
            if f'W{p}' in self.st['done'] and not c.get('keep', False):
                for f in ('G1', 'G2', 'X'):
                    if os.path.exists(self.path(f'{f}_{p}.npy')): os.remove(self.path(f'{f}_{p}.npy'))
            self.stage(f'Z{p}', lambda: self.Zmat(p))
            self.stage(f'H{p}', lambda: self.Hmat(p))
            if f'H{p}' in self.st['done'] and not c.get('keep', False):
                for f in ('Z', 'W'):
                    if os.path.exists(self.path(f'{f}_{p}.npy')): os.remove(self.path(f'{f}_{p}.npy'))
            self.stage(f'cert{p}', lambda: self.certify(p))
        if not c.get('keep', False):
            import shutil
            for sub in ('gv', 'bv'):
                if os.path.isdir(self.path(sub)): shutil.rmtree(self.path(sub))
        res = {str(p): self.st['cert'][str(p)] for p in c.get('parities', [0, 1])}
        allok = all(r.get('certified') for r in res.values()) and c.get('parities', [0, 1]) == [0, 1]
        self.st['RESULT'] = 'CERTIFIED' if allok else 'not certified'; self.save()
        tau = A(self.tau)
        self.log(f"RESULT {c['tag']}: {'CERTIFIED  Q(f) >= ' + S_(tau, 6) + ' ||f||^2 for every f supported in [-L, L], L = ' + S_(self.Sx.L_(), 10) if allok else 'not certified'}")
        return 0 if allok else 1

if __name__ == "__main__":
    R = Run(sys.argv[1])
    try:
        sys.exit(R.run())
    except OutOfTime as ex:
        R.save(); R.log(f"time budget used in {ex}; exit 3 (resume)"); sys.exit(3)
