# explore18 test T3: bulk symbol B, companion E = G*B, symmetrized model F = G(s)B(s) + G(1-s)B(1-s)
#   (functional equation F(s) = F(1-s), F = 2 Re(G B) on the critical line).
#   Index (edge) bound: N_0(F) >= W - 1, W = (1/pi) * change of arg(G B) along the line.
#   Bulk: zeros of B in Re s > 1/2 (kappa_B) are the defects. Euler products (finite X) have kappa_B = 0 and
#   bounded arg on the line; Dirichlet polynomials and the Davenport-Heilbronn combination of two Euler products do not.
import sys, time
from topo_common import *
from topo_common import _seg
mp.mp.dps = 15

def primes(X):
    return [p for p in range(2, X + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]

def chi5(n):
    return {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}[n % 5]

TAU = KAP

def make(kind, par):
    if kind == 'euler':
        P = primes(par)
        B = lambda s: mp.fprod([1 / (1 - mp.power(p, -s)) for p in P])
        G = lambda s: s * (s - 1) / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2)
    elif kind == 'zpart':
        B = lambda s: mp.fsum([mp.power(n, -s) for n in range(1, par + 1)])
        G = lambda s: s * (s - 1) / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2)
    elif kind == 'dhpart':
        a = {1: 1, 2: TAU, 3: -TAU, 4: -1, 0: 0}
        B = lambda s: mp.fsum([a[n % 5] * mp.power(n, -s) for n in range(1, par + 1)])
        G = lambda s: mp.power(5 / mp.pi, (s + 1) / 2) * mp.gamma((s + 1) / 2)
    elif kind == 'dheuler':
        P = [p for p in primes(par) if p != 5]
        def B(s):
            e1 = mp.fprod([1 / (1 - chi5(p) * mp.power(p, -s)) for p in P])
            e2 = mp.fprod([1 / (1 - mp.conj(chi5(p)) * mp.power(p, -s)) for p in P])
            return (1 - 1j * TAU) / 2 * e1 + (1 + 1j * TAU) / 2 * e2
        G = lambda s: mp.power(5 / mp.pi, (s + 1) / 2) * mp.gamma((s + 1) / 2)
    F = lambda s: G(s) * B(s) + G(1 - s) * B(1 - s)
    return G, B, F

def unwrap_arg(vals):
    out = [mp.arg(vals[0])]
    for k in range(1, len(vals)):
        out.append(out[-1] + mp.im(mp.log(vals[k] / vals[k - 1])))
    return out

def run(kind, par, T=250, s_lo=0.05, s_hi=0.95):
    G, B, F = make(kind, par)
    T0 = time.time()
    ts = [1 + mp.mpf(k) / 50 for k in range(int((T - 1) * 50) + 1)]
    Fl = [mp.re(F(mp.mpc(0.5, t))) for t in ts]
    n0 = sum(1 for k in range(len(Fl) - 1) if Fl[k] * Fl[k + 1] < 0)
    Ev = [G(mp.mpc(0.5, t)) * B(mp.mpc(0.5, t)) for t in ts]
    W = (unwrap_arg(Ev)[-1] - unwrap_arg(Ev)[0]) / mp.pi
    ab = unwrap_arg([B(mp.mpc(0.5, t)) for t in ts[::5]])
    osc = (max(ab) - min(ab)) / mp.pi
    mm = (ab[-1] - ab[0]) / (T - 1)
    ntot = zeros_in_rect(F, s_lo, s_hi, 1, T, 0.1)
    kB = zeros_in_rect(B, 0.5, 3, 1, T, 0.1)
    off = int(mp.nint(ntot)) - n0
    print(f"{kind:8s} {par:4d}: N_0 = {n0:4d}  N_total = {mp.nstr(ntot, 6):>7s}  off-line = {off:3d}  W = {mp.nstr(W, 6):>8s}  "
          f"kappa_B (zeros of B, Re>1/2) = {mp.nstr(kB, 5):>6s}  osc arg B/pi = {mp.nstr(osc, 4):>6s}  mean motion of arg B = {mp.nstr(mm, 3):>8s}  "
          f"[N_total - W = {mp.nstr(ntot - W, 4)}]  ({time.time() - T0:.0f} s)", flush=True)

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'euler':
        for X in (3, 7, 13, 29, 59, 101):
            run('euler', X)
    elif which == 'zpart':
        for N in (2, 5, 10, 30):
            run('zpart', N, s_lo=-1, s_hi=2)
    elif which == 'dhpart':
        for N in (5, 10, 30):
            run('dhpart', N, s_lo=-1, s_hi=2)
    elif which == 'dheuler':
        for X in (7, 13, 29, 59, 101):
            run('dheuler', X)
