#!/usr/bin/env python3
"""explore36 / lecam_multi: Assouad-type test.  In a block of B consecutive zeta ordinates, merge K disjoint consecutive
pairs (every m-th pair) into off-line pairs 1/2 +- delta + i gbar_j and move all other ordinates of the block (real moves,
|move| <= CAP_SP mean spacings) to minimise the Hilbert-Schmidt Le Cam distance on [-L, L]; report the relative
operator-norm distance of the Weil-window forms.  usage: lecam_multi.py <set> <B> <K-list> <fracs> [delta|rel]"""
import sys, os, math
import numpy as np
os.environ.setdefault('CAP_SP', '1.0')
import lecam1 as lc

def build(gam, c0, B, K):
    blk = np.arange(c0 - B // 2, c0 + B // 2)
    # K pairs spread evenly over the inner 80% of the block
    inner = blk[B // 10: B - B // 10 - 1]
    starts = np.unique(np.linspace(inner[0], inner[-1], K).astype(int))
    starts = [s for i, s in enumerate(starts) if i == 0 or s > starts[i - 1] + 1]
    removed = set()
    for s in starts:
        removed.add(s); removed.add(s + 1)
    movable = np.array([j for j in blk if j not in removed])
    return starts, movable

def sym(y, centre, gam, starts, movable, eps, delta):
    h = np.zeros(len(y), dtype=complex)
    for s in starts:
        gb = 0.5 * (gam[s] + gam[s + 1]) - centre
        h += 2.0 * np.cosh(delta * y) * np.exp(1j * gb * y) - np.exp(1j * (gam[s] - centre) * y) - np.exp(1j * (gam[s + 1] - centre) * y)
    th = gam[movable] - centre
    nz = np.nonzero(eps)[0]
    if len(nz):
        E = np.exp(1j * np.outer(y, th[nz]))
        h += (E * (np.exp(1j * np.outer(y, eps[nz])) - 1.0)).sum(axis=1)
    return h

def opn(a, centre, gam, starts, movable, eps, delta, nx):
    x, wx = lc.gl(nx, -a, a)
    d = x[:, None] - x[None, :]
    h = np.zeros((nx, nx), dtype=complex)
    for s in starts:
        gb = 0.5 * (gam[s] + gam[s + 1]) - centre
        h += 2.0 * np.cosh(delta * d) * np.exp(1j * gb * d) - np.exp(1j * (gam[s] - centre) * d) - np.exp(1j * (gam[s + 1] - centre) * d)
    th = gam[movable] - centre
    for k in np.nonzero(eps)[0]:
        h += np.exp(1j * (th[k] + eps[k]) * d) - np.exp(1j * th[k] * d)
    s = np.sqrt(wx)
    K = s[:, None] * h * s[None, :]
    K = 0.5 * (K + K.conj().T)
    lc.guard()
    return float(np.max(np.abs(np.linalg.eigvalsh(K))))

def main():
    which, B = sys.argv[1], int(sys.argv[2])
    Ks = [int(k) for k in sys.argv[3].split(',')]
    fracs = [float(f) for f in sys.argv[4].split(',')]
    dmode = sys.argv[5] if len(sys.argv) > 5 else 'rel'
    if which == 'z1e8':
        gam = lc.load_gram(f'{lc.ROOT}/ExploreDPP/logs/zeros_1e8.txt'); c0 = len(gam) // 2
    else:
        gam = lc.load_zeta_low(); c0 = int(np.searchsorted(gam, float(which[4:])))
    centre = gam[c0]
    Ls = 0.5 * math.log(centre / (2 * math.pi)); D = 2 * Ls / (2 * math.pi)
    cap = float(os.environ['CAP_SP']) / D
    delta = (0.5 / Ls) if dmode == 'rel' else float(dmode)
    print(f"# {which} t={centre:.4e} L*={Ls:.4f} block B={B} cap={os.environ['CAP_SP']} sp delta={delta:.4f}", flush=True)
    for K in Ks:
        starts, movable = build(gam, c0, B, K)
        for fr in fracs:
            L = fr * Ls
            thm = float(np.max(np.abs(gam[c0 - B // 2: c0 + B // 2] - centre))) + 1.0
            ny = int(max(400, 4.0 * thm * L / math.pi + 200))
            y, wq = lc.gl(ny, -L, L); w = wq * (L - np.abs(y)); sw = np.sqrt(w)
            eps = np.zeros(len(movable))
            cur = math.sqrt(float(np.sum(w * np.abs(sym(y, centre, gam, starts, movable, eps, delta)) ** 2)))
            n0 = cur; mu = None
            th = gam[movable] - centre
            for it in range(40):
                h = sym(y, centre, gam, starts, movable, eps, delta)
                J = (1j * y)[:, None] * np.exp(1j * np.outer(y, th + eps))
                A = np.vstack([(sw[:, None] * J).real, (sw[:, None] * J).imag]); r = np.concatenate([(sw * h).real, (sw * h).imag])
                AtA = A.T @ A; g = A.T @ r; del J, A
                if mu is None: mu = 1e-6 * np.trace(AtA) / len(th)
                ok = False
                for _ in range(25):
                    trial = eps + np.linalg.solve(AtA + mu * np.eye(len(th)), -g)
                    if np.max(np.abs(trial)) <= cap:
                        nt = math.sqrt(float(np.sum(w * np.abs(sym(y, centre, gam, starts, movable, trial, delta)) ** 2)))
                        if nt < cur:
                            eps, cur, ok = trial, nt, True; mu = max(mu / 4, 1e-14 * np.trace(AtA) / len(th)); break
                    mu *= 8.0
                lc.guard()
                if not ok: break
            a = 0.5 * L
            nx = int(min(500, max(160, 2.0 * thm * a / math.pi + 100)))
            o0 = opn(a, centre, gam, starts, movable, np.zeros_like(eps), delta, nx)
            o1 = opn(a, centre, gam, starts, movable, eps, delta, nx)
            print(f"K={len(starts):3d} offline zeros {2*len(starts)}/{B} ({2*len(starts)/B:.2f}) L/L*={fr:.2f} c={math.exp(L):9.2f} "
                  f"rho_HS={cur/n0:.2e} rho_op={o1/o0:.2e} maxmove={np.max(np.abs(eps))*D:.2f} sp rss={lc.rss_mb():.0f}MB", flush=True)

if __name__ == '__main__':
    main()
