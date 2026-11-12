#!/usr/bin/env python3
"""explore36 / lecam1: Le Cam two-point indistinguishability of an off-line pair under prime data up to X = e^L.

Data D_L (primes <= e^L, explicit formula) <=> zero-side Toeplitz symbol G(y) = sum_rho e^{i gamma_rho y} on [-L, L]
(Weil window a = L/2: Q(f) = int int G(x - x') f(x) conj f(x') dx dx', f in L^2(-a, a)).
Alternative Z1: two consecutive ordinates g1 < g2 near height t merged into the off-line pair 1/2 +- delta + i gbar
(gbar = (g1 + g2)/2), plus real moves eps_k of the M nearest other ordinates. Demodulated difference symbol
  h(y) = 2 cosh(delta y) - e^{i th1 y} - e^{i th2 y} + sum_k (e^{i (th_k + eps_k) y} - e^{i th_k y}),  th = gamma - gbar.
Le Cam distance (window a): eta = ||T_h|| (operator norm on L^2(-a, a)); proxy for fitting: HS norm,
  ||T_h||_HS^2 = int_{-L}^{L} |h(y)|^2 (L - |y|) dy.
Prediction: eta_min(M) -> 0 iff L < L*(t) = (1/2) log(t / 2 pi)  (X < sqrt(t / 2 pi)).
"""
import json, sys, resource, math
import numpy as np

CAP_MB = 280.0

def rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 1e6 if sys.platform == 'darwin' else r / 1e3

def guard():
    m = rss_mb()
    if m > CAP_MB:
        print(f"ABORT: RSS {m:.0f} MB > {CAP_MB}", flush=True)
        sys.exit(3)

ROOT = '/Users/adriansager/playground/riemann/Riemann/Deep'

def load_gram(path):
    """gram_win -zeros output: 'ZG n=.. tc=..' then 'ZZ n=.. offsets' -> absolute ordinates."""
    out = []
    tc = None
    with open(path) as fh:
        for line in fh:
            if line.startswith('ZG'):
                for tok in line.split():
                    if tok.startswith('tc='):
                        tc = float(tok[3:])
            elif line.startswith('ZZ') and tc is not None:
                vals = [float(x) for x in line.split()[2:]]
                out.extend(tc + v for v in vals)
    z = np.unique(np.round(np.array(sorted(out)), 7))
    return z

def load_zeta_low():
    return np.sort(np.load(f'{ROOT}/ExploreHelio/results/zeta_zeros.npy'))

def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def symbol(y, th, eps, delta, th1, th2):
    """demodulated difference symbol h(y) at nodes y (vector)."""
    h = 2.0 * np.cosh(delta * y) - np.exp(1j * th1 * y) - np.exp(1j * th2 * y)
    for k in range(len(th)):
        if eps[k] != 0.0:
            h += np.exp(1j * (th[k] + eps[k]) * y) - np.exp(1j * th[k] * y)
    return h

def hs_norm(y, w, h):
    return math.sqrt(float(np.sum(w * np.abs(h) ** 2)))

def opnorm(a, th, eps, delta, th1, th2, nx):
    """operator norm of T_h on L^2(-a, a) by Gauss-Legendre Nystrom (Hermitian kernel)."""
    x, wx = gl(nx, -a, a)
    d = x[:, None] - x[None, :]
    K = 2.0 * np.cosh(delta * d) - np.exp(1j * th1 * d) - np.exp(1j * th2 * d)
    for k in range(len(th)):
        if eps[k] != 0.0:
            K += np.exp(1j * (th[k] + eps[k]) * d) - np.exp(1j * th[k] * d)
    s = np.sqrt(wx)
    K = s[:, None] * K * s[None, :]
    K = 0.5 * (K + K.conj().T)
    ev = np.linalg.eigvalsh(K)
    guard()
    return float(np.max(np.abs(ev))), float(ev.min()), float(ev.max())

def compensate(gam, i0, delta, L, M, lam=1e-10, gn_iter=40):
    g1, g2 = gam[i0], gam[i0 + 1]
    gbar = 0.5 * (g1 + g2)
    th1, th2 = g1 - gbar, g2 - gbar
    # M nearest other ordinates (by index distance), excluding the pair
    lo = max(0, i0 - M // 2)
    hi = min(len(gam), i0 + 2 + M // 2)
    idx = [j for j in range(lo, hi) if j not in (i0, i0 + 1)]
    th = gam[idx] - gbar
    thmax = float(np.max(np.abs(th))) + 1.0
    ny = int(max(400, 4.0 * thmax * L / math.pi + 200))
    y, wq = gl(ny, -L, L)
    w = wq * (L - np.abs(y))
    eps = np.zeros(len(th))
    h0 = symbol(y, th, eps, delta, th1, th2)
    n0 = hs_norm(y, w, h0)
    sw = np.sqrt(w)
    D = math.log(max(gbar, 7.0) * float(__import__('os').environ.get('COND', '1')) / (2 * math.pi)) / (2 * math.pi)
    cap = float(__import__('os').environ.get('CAP_SP', '0.45')) / D   # moves below CAP_SP mean spacings (ordering kept)
    cur = n0
    mu = None
    for it in range(gn_iter):
        h = symbol(y, th, eps, delta, th1, th2)
        J = (1j * y)[:, None] * np.exp(1j * np.outer(y, th + eps))
        A = np.vstack([(sw[:, None] * J).real, (sw[:, None] * J).imag])
        r = np.concatenate([(sw * h).real, (sw * h).imag])
        AtA = A.T @ A
        g = A.T @ r
        del J, A
        if mu is None:
            mu = 1e-6 * np.trace(AtA) / len(th)
        accepted = False
        for _ in range(25):
            step = np.linalg.solve(AtA + mu * np.eye(len(th)), -g)
            trial = eps + step
            if np.max(np.abs(trial)) <= cap:
                nt = hs_norm(y, w, symbol(y, th, trial, delta, th1, th2))
                if nt < cur:
                    eps, cur = trial, nt
                    mu = max(mu / 4.0, 1e-14 * np.trace(AtA) / len(th))
                    accepted = True
                    break
            mu *= 8.0
        guard()
        if not accepted:
            break
    h1 = symbol(y, th, eps, delta, th1, th2)
    n1 = hs_norm(y, w, h1)
    return dict(gbar=float(gbar), gap=float(g2 - g1), M=len(th), ny=ny, hs0=n0, hs1=n1,
                rho_hs=n1 / n0, maxeps=float(np.max(np.abs(eps))), rmseps=float(np.sqrt(np.mean(eps ** 2)))), (th, eps, th1, th2)

def main():
    which = sys.argv[1]
    delta_mode = sys.argv[2] if len(sys.argv) > 2 else 'rel'
    Ms = [int(m) for m in (sys.argv[3] if len(sys.argv) > 3 else '40,100,200').split(',')]
    fracs = [float(f) for f in (sys.argv[4] if len(sys.argv) > 4 else '0.5,0.7,0.85,1.0,1.15,1.3,1.6,2.0').split(',')]
    do_op = (len(sys.argv) > 5 and sys.argv[5] == 'op')
    if which == 'z1e8':
        gam = load_gram(f'{ROOT}/ExploreDPP/logs/zeros_1e8.txt')
    elif which == 'z1e11':
        gam = load_gram(f'{ROOT}/ExploreDPP/logs/zeros_1e11.txt')
    elif which.startswith('zlow'):
        z = load_zeta_low()
        tt = float(which[4:])
        gam = z
    else:
        raise SystemExit('which?')
    if which.startswith('zlow'):
        i0 = int(np.searchsorted(gam, tt))
    else:
        i0 = len(gam) // 2
    i0 += int(sys.argv[6]) if len(sys.argv) > 6 else 0
    t = gam[i0]
    D = math.log(t / (2 * math.pi)) / (2 * math.pi)
    Lstar = 0.5 * math.log(t / (2 * math.pi))
    print(f"# set={which} N={len(gam)} range=[{gam[0]:.3f},{gam[-1]:.3f}] i0={i0} t={t:.6f} density D={D:.4f} "
          f"L*={Lstar:.4f} X*=sqrt(t/2pi)={math.exp(Lstar):.3f} gap={gam[i0+1]-gam[i0]:.4f} (mean {1/D:.4f})", flush=True)
    for M in Ms:
        if M > len(gam) - 4:
            continue
        for fr in fracs:
            L = fr * Lstar
            delta = (0.5 / Lstar) if delta_mode == 'rel' else float(delta_mode)
            res, (th, eps, th1, th2) = compensate(gam, i0, delta, L, M)
            line = (f"M={res['M']:4d} L/L*={fr:5.2f} L={L:7.3f} c=e^L={math.exp(L):10.3f} delta={delta:.4f} "
                    f"rho_HS={res['rho_hs']:.3e} maxmove/spacing={res['maxeps']*D:.3e} rms={res['rmseps']*D:.2e}")
            if do_op:
                a = 0.5 * L
                thmax = float(np.max(np.abs(th))) + 1.0
                nx = int(min(900, max(200, 2.0 * thmax * a / math.pi + 120)))
                op0 = opnorm(a, th, np.zeros_like(eps), delta, th1, th2, nx)
                op1 = opnorm(a, th, eps, delta, th1, th2, nx)
                line += f" op0={op0[0]:.3e} op1={op1[0]:.3e} rho_op={op1[0]/op0[0]:.3e} nx={nx}"
            print(line + f" rss={rss_mb():.0f}MB", flush=True)

if __name__ == '__main__':
    main()
