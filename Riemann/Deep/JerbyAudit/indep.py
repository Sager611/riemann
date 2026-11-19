# JerbyAudit independent recomputation (no Arb, no JerbyFals code): own Stirling theta, own Riemann-Siegel Z and Z'
# (C0..C2 terms, termwise derivative), mpmath siegelz / siegelz(derivative=1), Richardson central differences of siegelz.
# usage: indep.py cases | fidelity
import sys, time
import mpmath as mp
mp.mp.dps = 40
PI = mp.pi
def theta_S(t):   # Stirling series of arg Gamma(1/4+it/2) - (t/2) log pi, branch with theta(0)=0 (valid for large t)
    t = mp.mpf(t)
    return t/2*mp.log(t/(2*PI)) - t/2 - PI/8 + 1/(48*t) + mp.mpf(7)/(5760*t**3) + mp.mpf(31)/(80640*t**5) + mp.mpf(127)/(430080*t**7)
def dtheta_S(t):
    t = mp.mpf(t)
    return mp.log(t/(2*PI))/2 - 1/(48*t**2) - mp.mpf(21)/(5760*t**4) - mp.mpf(155)/(80640*t**6)
def gram(n, th=theta_S):
    x0 = 2*PI*mp.exp(1 + mp.lambertw((8*n + 1)/(8*mp.e)).real)
    return mp.findroot(lambda t: th(t) - n*PI, x0)
def Psi(p): return mp.cos(2*PI*(p*p - p - mp.mpf(1)/16))/mp.cos(2*PI*p)
def ZRS(t):
    """own Riemann-Siegel: returns (Z, Z') with C0, C1, C2 and their termwise t-derivatives"""
    t = mp.mpf(t); a = mp.sqrt(t/(2*PI)); N = int(mp.floor(a)); p = a - N
    th = theta_S(t); thp = dtheta_S(t)
    S = mp.mpf(0); dS = mp.mpf(0)
    for k in range(1, N + 1):
        lk = mp.log(k); ph = th - t*lk; w = 1/mp.sqrt(k)
        S += w*mp.cos(ph); dS -= w*mp.sin(ph)*(thp - lk)
    c = mp.taylor(Psi, p, 7); d = [c[j]*mp.factorial(j) for j in range(8)]
    C0, C0p = d[0], d[1]
    C1, C1p = -d[3]/(96*PI**2), -d[4]/(96*PI**2)
    C2, C2p = d[2]/(64*PI**2) + d[6]/(18432*PI**4), d[3]/(64*PI**2) + d[7]/(18432*PI**4)
    sg = (-1)**(N - 1)
    R = sg*(a**mp.mpf(-0.5)*C0 + a**mp.mpf(-1.5)*C1 + a**mp.mpf(-2.5)*C2)
    dRda = sg*(-a**mp.mpf(-1.5)*C0/2 + a**mp.mpf(-0.5)*C0p - 3*a**mp.mpf(-2.5)*C1/2 + a**mp.mpf(-1.5)*C1p - 5*a**mp.mpf(-3.5)*C2/2 + a**mp.mpf(-2.5)*C2p)
    return 2*S + R, 2*dS + dRda/(4*PI*a)
def ZFD(t, h=mp.mpf('1e-6')):   # Richardson-extrapolated central difference of mpmath siegelz
    D = lambda h: (mp.siegelz(t + h) - mp.siegelz(t - h))/(2*h)
    return (4*D(h/2) - D(h))/3
def gb(k, Z): return 'good' if (-1)**k*Z > 0 else 'BAD'
mode = sys.argv[1]
if mode == 'cases':
    # page / certify_cx*.log values: n, printed g_n midpoint, Arb mu midpoint and radius
    cases = [(1565734083, '568009639.698978187581685032969', '3.968215684906503', 2.30e-16),
             (2565007238, '906089969.752693821200280999797', '3.944872238693607', 4.89e-16),
             (3315629581, '1155454823.10758508544277960703', '3.922880048686967', 1.05e-16)]
    for n, gstr, mustr, murad in cases:
        t0 = time.time()
        print(f'=== n = {n} ({"odd" if n % 2 else "even"})')
        for k in range(n - 2, n + 3):
            g = gram(k); gm = mp.findroot(lambda t: mp.siegeltheta(t) - k*PI, g)
            Z, dZ = ZRS(g); Zm = mp.siegelz(g)
            print(f'  k={k} g(own theta)={mp.nstr(g, 30)} |g - g(siegeltheta)|={mp.nstr(abs(g - gm), 3)} Z_RS={mp.nstr(Z, 18)} Z_mp={mp.nstr(Zm, 18)} |diff|={mp.nstr(abs(Z - Zm), 3)} {gb(k, Zm)}')
            if k == n:
                gc = mp.mpf(gstr); th_err = theta_S(gc) - n*PI; th_err2 = mp.siegeltheta(gc) - n*PI
                dZm = mp.siegelz(g, derivative=1); dZf = ZFD(g)
                mu_rs, mu_mp, mu_fd = abs(dZ/Z), abs(dZm/Zm), abs(dZf/Zm)
                print(f'  log g_n midpoint: |g - g_log|={mp.nstr(abs(g - gc), 3)}; theta_S(g_log)-n*pi={mp.nstr(th_err, 3)}; siegeltheta(g_log)-n*pi={mp.nstr(th_err2, 3)}; theta_S(g_log)-siegeltheta(g_log)={mp.nstr(theta_S(gc) - mp.siegeltheta(gc), 3)}')
                print(f"  Z'_RS={mp.nstr(dZ, 20)}  Z'_mp(derivative=1)={mp.nstr(dZm, 20)}  Z'_FD={mp.nstr(dZf, 20)}")
                print(f'  mu_RS={mp.nstr(mu_rs, 20)} mu_mp={mp.nstr(mu_mp, 20)} mu_FD={mp.nstr(mu_fd, 20)}')
                print(f'  vs page mu={mustr} +/- {murad:.2e}: |mu_mp - page|={mp.nstr(abs(mu_mp - mp.mpf(mustr)), 3)}  4 - mu = {mp.nstr(4 - mu_mp, 10)}  4|Z| - |Z\'| = {mp.nstr(4*abs(Zm) - abs(dZm), 12)}')
                print(f"  theta'(g_n) = {mp.nstr(dtheta_S(g), 15)}")
        print(f'  [{time.time() - t0:.1f} s]', flush=True)
elif mode == 'fidelity':
    t0 = time.time()
    bad = []
    for k in range(0, 141):
        g = mp.findroot(lambda t: mp.siegeltheta(t) - k*PI, gram(k) if k > 0 else mp.mpf(17.8))
        if (-1)**k*mp.siegelz(g) <= 0: bad.append(k)
    print(f'bad Gram points with 0 <= n <= 140 (standard theta, theta(g_n) = n pi): {bad}   g_0 = {mp.nstr(mp.findroot(lambda t: mp.siegeltheta(t), 17.8), 12)}')
    dth7 = mp.diff(mp.siegeltheta, 7); tstar = mp.findroot(lambda t: mp.diff(mp.siegeltheta, t), 6.3)
    print(f"theta'(7) = {mp.nstr(dth7, 8)} > 0; theta' = 0 at t* = {mp.nstr(tstar, 12)}, theta(t*) = {mp.nstr(mp.siegeltheta(tstar), 8)}")
    for k, note in [(126, 'first bad (2406 l.214, 743)'), (134, 'second bad (2406 l.743)'), (6708, 'isolated bad (2310 l.1265)'),
                    (730119, 'isolated bad, mu=4.46023 (2310 l.1270)'), (9807962, 'corrupt, mu=0.0750883, block {9807960..9807963} (2406 l.757, 771)'),
                    (69784845, 'conjfals min isolated mu=4.21553 (RECORDS)')]:
        gs = [mp.findroot(lambda t: mp.siegeltheta(t) - j*PI, gram(j)) for j in range(k - 2, k + 3)]
        Zs = [mp.siegelz(x) for x in gs]
        pat = ''.join('G' if (-1)**j*z > 0 else 'B' for j, z in zip(range(k - 2, k + 3), Zs))
        mu = abs(mp.siegelz(gs[2], derivative=1)/Zs[2])
        print(f'n={k}: pattern g_(n-2..n+2) = {pat}; mu(g_n) = {mp.nstr(mu, 9)}   [{note}]', flush=True)
    print(f'[{time.time() - t0:.1f} s]')
