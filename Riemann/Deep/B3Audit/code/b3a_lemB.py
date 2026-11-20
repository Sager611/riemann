# b3audit F3: independent direct-quadrature falsification of B3Cert Lemma B (exact-profile dual term T(n)) and Lemma R.
# Independent of B3Cert's routes (FFT profile / Laguerre / Lambda-coordinates): A_{-1}(xi) by direct trapezoid in y with an
# aliasing-controlled step (Poisson: error ~ exp(-margin^2/4)); T_sigma(n) by direct trapezoid in delta = 2 pi xi/T - 1;
# Hankel functions of complex order d by their asymptotic series (checked against mpmath). Cutoff chi~ = 1 on the computed
# range (|y0| <= 6.8 => u in [e^{-13.6/H}, e^{13.6/H}] subset [2/3, 3/2] for H >= 34; truncation ~ e^{-46}).
# T(n) = -pi int A(xi) xi^{-z/2} B(X) dxi, B = sin(pi d/2) J_d + cos(pi d/2) Y_d = c_+ H1_d + c_- H2_d, X = 4 pi sqrt(w xi), w = n/(pq).
# usage: b3a_lemB.py T:theta[,T:theta...]  (prints one block per (T, theta, shift set))
import sys, math, time, resource, numpy as np, mpmath as mp

def hankel_asym(nu, X, K=12):
    mu = 4*nu*nu; a = [1.0+0j]
    for k in range(1, K): a.append(a[-1]*(mu - (2*k-1)**2)/(k*8))
    s1 = np.zeros(X.shape, complex); s2 = np.zeros(X.shape, complex); Xk = np.ones(X.shape)
    for k in range(K):
        s1 += (1j)**k*a[k]/Xk; s2 += (-1j)**k*a[k]/Xk; Xk = Xk*X
    ph = X - nu*np.pi/2 - np.pi/4; pre = np.sqrt(2/(np.pi*X))
    return pre*np.exp(1j*ph)*s1, pre*np.exp(-1j*ph)*s2

def profile(T, H, beta, dl, Y=6.5, margin=30.0):
    K = T/H; dmax = float(np.max(np.abs(dl)))
    wmax = 2*K*max(abs(1-(1+s)*math.exp(2*t/H)) for s in (-dmax, dmax) for t in (-Y, Y))
    h = 2*math.pi/(wmax + margin); ny = int(2*Y/h) + 2
    y = np.linspace(-Y, Y, ny); hy = y[1]-y[0]; e = 2*y/H; em1 = np.expm1(e); c0 = -(em1 - e)
    amp = np.exp(-y*y + (1-2*beta)*y/H)*hy/math.sqrt(math.pi)
    out = np.empty(len(dl), complex); CH = max(1, int(1.5e6//ny))
    for i in range(0, len(dl), CH):
        d = dl[i:i+CH, None]
        out[i:i+CH] = np.exp(1j*(T*(c0[None, :] - d*em1[None, :]))) @ amp
    xi = T*(1+dl)/(2*math.pi)
    return out*xi**(-beta), ny

def profile_mp(T, H, beta, dl):   # mpmath reference for one delta (y-integral on 400 subintervals)
    mp.mp.dps = 25; T = mp.mpf(T); H = mp.mpf(H); b = mp.mpc(beta); d = mp.mpf(dl)
    f = lambda y: mp.exp(-y*y + (1-2*b)*y/H + 1j*(T*(2*y/H - mp.expm1(2*y/H)) - T*d*mp.expm1(2*y/H)))
    s = mp.quad(f, mp.linspace(-6.5, 6.5, 401))
    return complex(s*(T*(1+d)/(2*mp.pi))**(-b)/mp.sqrt(mp.pi))

def run(T, th, al, be, rs):
    H = T**th; K = T/H; tau = T/H**2; d = al + be; z = al - be
    yb = 6.8; dlo, dhi = math.exp(-2*yb/H) - 1, math.exp(2*yb/H) - 1
    X0max = 2*math.sqrt(2*math.pi)*math.sqrt(max(rs)*tau*T)
    wd = T*max(abs(dlo), dhi)/(1 - abs(dlo)) + X0max/(2*math.sqrt(1 + dlo))
    hd = 2*math.pi/(wd + 12*K); nd = int((dhi - dlo)/hd) + 2
    dl = np.linspace(dlo, dhi, nd); hdl = dl[1] - dl[0]
    t0 = time.time(); A, ny = profile(T, H, be, dl)
    xi = T*(1+dl)/(2*math.pi); dxi = T*hdl/(2*math.pi)
    lr = abs(np.sum(A*xi**(be))*dxi/(H/(2*math.sqrt(math.pi))) - 1)   # Lemma R at s = 1 + beta (w = 1): M/M_lead - 1
    edge = float(max(abs(A[0]), abs(A[-1]))/np.max(np.abs(A)))
    print(f"T={T:.0e} th={th} H={H:.2f} K={K:.1f} tau={tau:.3f} shifts=({al},{be}) n_delta={nd} n_y={ny} prof {time.time()-t0:.1f}s "
          f"edge|A|/max={edge:.1e} max|A|sqrt(tau)={np.max(np.abs(A))*math.sqrt(tau):.3f} |M/M_lead-1|={lr:.2e} (H/T={H/T:.1e})", flush=True)
    cp, cm = -0.5j*np.exp(1j*np.pi*d/2), 0.5j*np.exp(-1j*np.pi*d/2)
    rows = []
    for r in rs:
        w = r*tau; X0 = 2*math.sqrt(2*math.pi)*math.sqrt(w*T); cr = X0/K
        X = 4*math.pi*np.sqrt(w*xi); H1, H2 = hankel_asym(d, X)
        base = -math.pi*A*xi**(-z/2)*dxi
        Tp = np.sum(base*cp*H1); Tm = np.sum(base*cm*H2); scale = H/math.sqrt(X0)
        pred = math.exp(-math.pi*r/2)/(2*math.sqrt(2))*abs((T/(2*math.pi))**(-d/2))
        rows.append((r, X0, cr, cr <= H/46, abs(Tp)/scale, abs(Tm)/scale, abs(Tp+Tm)/scale, pred))
        print(f"  r={r:<6} X0={X0:.4g} c_r={cr:.3f} bulk={'Y' if cr <= H/46 else 'n'} |T+|/HX0^-1/2={abs(Tp)/scale:.5f} "
              f"|T-|={abs(Tm)/scale:.5f} |T|={abs(Tp+Tm)/scale:.5f} pred/branch={pred:.5f}", flush=True)
    return A, dl, rows

if __name__ == "__main__":
    jobs = [tuple(map(float, j.split(':'))) for j in sys.argv[1].split(',')]
    shifts = [(0.0, 0.0), (0.05+0.1j, -0.03+0.07j), (0.2-0.1j, 0.03+0.2j)]
    rs = [0.003, 0.03, 0.1, 0.3, 1.0, 2.0, 4.0, 8.0]
    if len(sys.argv) > 2 and sys.argv[2] == 'check':   # mpmath checks of the y-trapezoid and of the Hankel series
        T, th = jobs[0]; H = T**th
        for be in [0.0, -0.03+0.07j]:
            dl = np.array([-3/H, -0.5/H, 0.0, 1/H, 4/H]); A, _ = profile(T, H, be, dl)
            for i, x in enumerate(dl):
                Am = profile_mp(T, H, be, x); print(f"check A: T={T:.0e} th={th} beta={be} delta*H={x*H:.2f} trap={A[i]:.10e} mp={Am:.10e} rel={abs(A[i]-Am)/abs(Am):.1e}", flush=True)
        for nu in [0.0, 0.02+0.17j]:
            for Xv in [60.0, 800.0, 2e4]:
                h1, h2 = hankel_asym(nu, np.array([Xv])); m1 = complex(mp.hankel1(nu, Xv)); m2 = complex(mp.hankel2(nu, Xv))
                print(f"check H: nu={nu} X={Xv} relH1={abs(h1[0]-m1)/abs(m1):.1e} relH2={abs(h2[0]-m2)/abs(m2):.1e}", flush=True)
        sys.exit(0)
    for T, th in jobs:
        for al, be in shifts:
            run(T, th, al, be, rs)
    print(f"maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20:.1f} MB", flush=True)
