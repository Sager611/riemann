# t3_turing.py (explore21): the "Turing critical mode" of a neural field whose lateral kernel is Weil's distribution on the
# window [-L/2, L/2], L = ln c, i.e. the ground state of the window Weil form (CCM's even ground state), for zeta and for the
# Davenport-Heilbronn function (DH; explicit formula, Gamma((s+1)/2), conductor 5; Weil form from explore19's wmat.py).
# Krein / Caratheodory-Toeplitz: Q - lambda_min ||.||^2 >= 0 is singular, so the exact ground state's Fourier transform has only
# real zeros, whatever the sign of lambda_min. Test: zeros of the Galerkin ground state's FT (real sign changes vs argument
# principle in a box), their match to the true zeros, and DH beyond its first failure (c = 40 > 30.7), whose off-line pair
# 0.808517 + 85.699348i (z = 85.699348 -+ 0.308517i) cannot be matched by real zeros.
# Usage: t3_turing.py kind c Zmax [Nadd]
import sys, os, math, time
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreAntenna/code')
from flint import arb, acb, ctx, fmpq
import mpmath as mp
import wmat

def run(kind, c, Zmax, Nadd=20):
    t0 = time.time()
    bits = (int(4*math.pi*c*1.4427) if kind == 'zeta' else int(4*math.pi*c/5*1.4427)) + 160
    ctx.prec = bits; mp.mp.prec = bits
    N = int(math.ceil(0.6*c*math.log(c))) + Nadd
    F = wmat.Form(fmpq(int(round(c*64)), 64), kind)
    T, idx = F.block(N, 'even')
    M = wmat.to_mp(T, mp.mp)
    E, Q = mp.mp.eigsy(M)
    order = sorted(range(len(E)), key=lambda i: E[i])
    i0 = order[0]
    v = [Q[r, i0] for r in range(len(E))]
    L = mp.log(mp.mpf(int(round(c*64)))/64)
    w = [2*mp.pi*j/L for j in idx]
    Nj = [mp.sqrt(1/L) if j == 0 else mp.sqrt(2/L) for j in idx]
    cj = [v[k]*Nj[k]*(1 if idx[k] % 2 == 0 else -1) for k in range(len(idx))]
    w2 = [x*x for x in w]
    def FT(z):
        z = mp.mpmathify(z)
        s = mp.mpf(0)
        for k in range(len(cj)):
            s += cj[k]/(z*z - w2[k])
        return mp.sin(z*L/2)*2*z*s
    print('%s c=%g L=%.4f N=%d bits=%d  lambda_min=%s  lambda_2=%s  build+eig %.0fs  w_N=%.1f'
          % (kind, c, float(L), N, bits, mp.nstr(E[i0], 6), mp.nstr(E[order[1]], 6), time.time()-t0, float(w[-1])), flush=True)
    # real zeros by sign changes (avoid exact w_j), refine by bisection
    zs = []; step = mp.mpf('0.01'); z = mp.mpf('0.0123'); fz = FT(z)
    while z < Zmax:
        z2 = z + step; f2 = FT(z2)
        if mp.sign(f2) != mp.sign(fz):
            a, b, fa = z, z2, fz
            for _ in range(40):
                m = (a + b)/2; fm = FT(m)
                if mp.sign(fm) == mp.sign(fa): a, fa = m, fm
                else: b = m
            zs.append((a + b)/2)
        z, fz = z2, f2
    # argument principle on the box [0.0123, Zmax] x [-H, H]
    H = mp.mpf(2)
    def seg(p, q, n):
        return [p + (q - p)*k/n for k in range(n + 1)]
    path = seg(mp.mpc(0.0123, -H), mp.mpc(Zmax, -H), int(Zmax/0.05)) + seg(mp.mpc(Zmax, -H), mp.mpc(Zmax, H), 200)[1:] \
         + seg(mp.mpc(Zmax, H), mp.mpc(0.0123, H), int(Zmax/0.05))[1:] + seg(mp.mpc(0.0123, H), mp.mpc(0.0123, -H), 200)[1:]
    tot = mp.mpf(0); prev = FT(path[0])
    for p in path[1:]:
        cur = FT(p); d = mp.arg(cur/prev)
        tot += d; prev = cur  # step 0.05 on the box edges, at distance >= 2 from the real zeros
    nbox = float(tot/(2*mp.pi))
    print('  real zeros in (0, %g): %d   zeros in box |Im z|<2 (argument principle): %.3f   %.0fs' % (Zmax, len(zs), nbox, time.time()-t0), flush=True)
    return zs

def true_zeros(kind, Zmax):
    if kind == 'zeta':
        out = []; n = 1
        while True:
            g = mp.im(mp.zetazero(n))
            if g > Zmax: break
            out.append(g); n += 1
        return out
    # DH on-line zeros: sign changes of the completed real function (5/pi)^{s/2} Gamma((s+1)/2) f(s) on s = 1/2 + it
    mp.mp.dps = 20
    kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
    A = [1, kap, -kap, -1]
    def lam(t):
        s = mp.mpc(0.5, t)
        f = sum(A[r-1]*mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5))*mp.power(5, -s)
        return (mp.power(5/mp.pi, s/2)*mp.gamma((s + 1)/2)*f)
    out = []; t = mp.mpf('0.01'); prev = mp.re(lam(t)); imax = 0
    while t < Zmax:
        t2 = t + mp.mpf('0.02'); cur = lam(t2); imax = max(imax, abs(mp.im(cur))/(abs(cur) + 1e-300))
        if mp.sign(mp.re(cur)) != mp.sign(prev):
            out.append(mp.findroot(lambda x: mp.re(lam(x)), (t, t2), solver='anderson'))
        t, prev = t2, mp.re(cur)
    print('  DH on-line zeros below %g: %d (max |Im|/|.| on the line %.1e)' % (Zmax, len(out), float(imax)))
    return out

if __name__ == '__main__':
    kind, c, Zmax = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    Nadd = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    zs = run(kind, c, Zmax, Nadd)
    bits = mp.mp.prec
    tz = true_zeros(kind, Zmax)
    mp.mp.prec = bits
    print('  k   ground-state FT zero   nearest true on-line zero   difference')
    for k, z in enumerate(zs):
        if tz:
            nz = min(tz, key=lambda g: abs(g - z))
            print('  %2d  %20s  %20s  %10s' % (k+1, mp.nstr(z, 12), mp.nstr(nz, 12), mp.nstr(z - nz, 3)))
        else:
            print('  %2d  %20s' % (k+1, mp.nstr(z, 12)))
