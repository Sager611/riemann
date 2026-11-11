# dh2.py -- jerbyconf: fast (double precision) scan of the DH section family S_N(t;r) = cos th(t) + r sum_{k=2}^N c_k k^-1/2 cos(th(t) - t log k)
# for N = 2..NMAX and two theta variants: exact th(t) = Im log Gamma(3/4+it/2) - (t/2) log(pi/5), and Jerby's (Zenodo g_n^DH satisfy
# th_J(g_n) = n pi with th_J = th - 1/(48 t), i.e. Stirling without the 1/t term).  Reports mu_sec = |S'/S|(g_44;1), the same with the
# Gamma-modulus drift d/dt log|G| added (unnormalised), and r_* (first zero of Delta(r) = S(g(r);r), g(r) the continued extremum).
import math, sys, mpmath as mp
mp.mp.dps = 20
s5 = math.sqrt(5); kap = (math.sqrt(10 - 2*s5) - 2)/(s5 - 1); c = [0.0, 1.0, kap, -kap, -1.0]
def th_ex(t): return float(mp.im(mp.loggamma(mp.mpf(3)/4 + 0.5j*t)) - t/2*mp.log(mp.pi/5))
def thp_ex(t): return float(mp.re(mp.digamma(mp.mpf(3)/4 + 0.5j*t))/2 - mp.log(mp.pi/5)/2)
def thpp_ex(t): return float(-mp.im(mp.psi(1, mp.mpf(3)/4 + 0.5j*t))/4)
VAR = {'exact': (th_ex, thp_ex, thpp_ex), 'jerby': (lambda t: th_ex(t) - 1/(48*t), lambda t: thp_ex(t) + 1/(48*t*t), lambda t: thpp_ex(t) - 1/(24*t**3))}
def S(t, r, N, th, thp, thpp):
    a, b, d = th(t), thp(t), thpp(t)
    v = math.cos(a); v1 = -math.sin(a)*b; v2 = -math.cos(a)*b*b - math.sin(a)*d
    for k in range(2, N + 1):
        w = r*c[k % 5]/math.sqrt(k)
        if w == 0: continue
        L = math.log(k); ph = a - t*L; q = b - L
        v += w*math.cos(ph); v1 += -w*math.sin(ph)*q; v2 += -w*(math.cos(ph)*q*q + math.sin(ph)*d)
    return v, v1, v2
G_drift = float(mp.diff(lambda t: mp.log(abs(mp.gamma((mp.mpf(3)/2 + 1j*t)/2))), 85.5848))
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
for name, (th, thp, thpp) in VAR.items():
    g = 85.58
    for _ in range(30): g -= (th(g) - 44*math.pi)/thp(g)
    print(f'## theta variant {name}: g_44 = {g:.12f}')
    print('| N | S(g44;1) | mu_sec | mu_sec + dlog|G| | r_* |')
    for N in range(2, NMAX + 1):
        if N % 5 == 0: continue          # c_N = 0: same as N-1
        v, v1, _ = S(g, 1.0, N, th, thp, thpp)
        # track the extremum: Newton on S' = 0
        x = g; prev = None; rstar = None
        R = [i/400 for i in range(401)]
        for r in R:
            for _ in range(20):
                _, d1, d2 = S(x, r, N, th, thp, thpp); x -= d1/d2
            D = S(x, r, N, th, thp, thpp)[0]
            if prev is not None and prev[1]*D < 0:
                lo, hi, xl = prev[0], r, prev[2]
                for _ in range(40):
                    m = 0.5*(lo + hi); xm = xl
                    for _ in range(20):
                        _, d1, d2 = S(xm, m, N, th, thp, thpp); xm -= d1/d2
                    Dm = S(xm, m, N, th, thp, thpp)[0]
                    if Dm*prev[1] > 0: lo, xl = m, xm
                    else: hi = m
                rstar = 0.5*(lo + hi); break
            prev = (r, D, x)
        print(f'| {N} | {v:.6f} | {abs(v1/v):.5f} | {abs(v1/v + G_drift):.5f} | ' + (f'{rstar:.5f}' if rstar else '–') + ' |'); sys.stdout.flush()
