
# End-to-end numerical check of Casana-Giner Lemma 3.1, Lemma 3.2 and Proposition 4.1 at sample heights,
# using Odlyzko's first 1e5 zeros (RH holds there, so b_rho = 0) and mpmath for zeta, loggamma, polygamma.
import json, math, sys
import numpy as np
import mpmath as mp
mp.mp.dps = 40
K = json.load(open('Riemann/Deep/STCAudit/kernels.json'))
g = np.loadtxt('.scratch/stcaudit/zeros1')
GMAX = g[-1]
def M(T): return T/(2*math.pi)*math.log(T/(2*math.pi*math.e))
def A(s):
    s = mp.mpc(s); T = s.imag; sig = s.real
    return mp.atan(T/sig) + mp.atan(T/(sig-1)) - T/2*mp.log(mp.pi) + mp.im(mp.loggamma(s/2))
def imlogzeta(s):
    z = mp.zeta(s); return mp.arg(z)   # |Im log zeta| <= log zeta(sigma) < pi for sigma >= 1.35: principal arg is the branch
def zsum(fun, T, d):
    # sum over all zeros rho (both signs of gamma) of fun(t_rho) with b_rho = 0, plus density tail beyond GMAX
    t1 = T - g; t2 = T + g
    s = float(np.sum(fun(t1)) + np.sum(fun(t2)))
    # tail: zeros with |gamma| > GMAX, density (1/2pi) log(gamma/2pi)
    tail = mp.quad(lambda x: mp.log(x/(2*mp.pi))/(2*mp.pi)*(fun_mp(fun, T - x) + fun_mp(fun, T + x)), [GMAX, 10*GMAX, mp.inf])
    return s, float(tail)
def fun_mp(fun, t):
    return mp.mpf(float(fun(np.array([float(t)]))[0]))
out = []
for T in [30.5, 100.3, 1000.7, 5000.2, 20000.4]:
    N = int(np.sum(g <= T))
    for name, k in K.items():
        d = k['d']; s1 = mp.mpc(0.5+d, T); s2 = mp.mpc(0.5+2*d, T)
        argxi1 = A(s1) + imlogzeta(s1); argxi2 = A(s2) + imlogzeta(s2)
        main = float((2*argxi1 - argxi2)/mp.pi)
        f0 = lambda t: 4*np.arctan(d/t) - 2*np.arctan(2*d/t)
        fs, ft = zsum(f0, T, d)
        ident = main + (fs+ft)/(2*math.pi)
        # Lemma 3.2 check and Prop 4.1
        zz = lambda s: mp.zeta(s, 1, 1)/mp.zeta(s)
        dz = {}
        Zform = []; Zzero = []
        for (X, j), c in zip(k['B'], k['c']):
            if j not in dz:
                dz[j] = mp.diff(zz, s1, j-1) if j > 1 else zz(s1)
            Dj = (1/s1 + 1/(s1-1) - mp.log(mp.pi)/2 + mp.digamma(s1/2)/2) if j == 1 else (s1**(-j) + (s1-1)**(-j) + (-1)**(j-1)/mp.factorial(j-1)/2**j*mp.psi(j-1, s1/2))
            dser = (-1)**j*dz[j]          # = sum Lambda(n) (log n)^(j-1) n^(-s1)
            PhiX = (lambda w: mp.re(w)) if X == 'R' else (lambda w: -mp.im(w))
            zf = 2*PhiX(Dj) - 2/mp.factorial(j-1)*PhiX(dser)
            Gx = (lambda t, j=j, X=X: 2*(((d+1j*t)**(-j)).real if X == 'R' else -((d+1j*t)**(-j)).imag))
            zs, zt = zsum(Gx, T, d)
            Zform.append(float(zf)); Zzero.append(zs+zt)
        c = np.array(k['c']); eps = np.array([1 if X == 'R' else -1 for X, j in k['B']])
        up = main - M(T) + float(np.dot(c, Zform))/(2*math.pi)
        lo = main - M(T) - float(np.dot(eps*c, Zform))/(2*math.pi)
        maxdev = max(abs(a-b) for a, b in zip(Zform, Zzero))
        A1 = math.ceil(k['C1']*1e5)/1e5
        rec = dict(T=T, kernel=name, N=N, ident_resid=ident-N, lemma32_maxdev=maxdev, NminusM=N-M(T), upper=up, lower=lo,
                   ok=(lo <= N-M(T) <= up), slack_up=up-(N-M(T)), slack_lo=(N-M(T))-lo)
        out.append(rec)
        print('T=%-8g %-5s N=%-5d ident_resid=%+.2e L32dev=%.1e  N-M=%.4f in [%.4f, %.4f] %s' % (T, name, N, ident-N, maxdev, N-M(T), lo, up, 'OK' if rec['ok'] else 'FAIL'))
        sys.stdout.flush()
json.dump(out, open('Riemann/Deep/STCAudit/identity_check.json', 'w'), indent=1)
