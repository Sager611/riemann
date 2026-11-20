"""nbdir: checks of Burnol's liminf for real chi (DERIVATION §1-5) and the finite-lambda lower bound.
A identities; B Lemma B-chi (closed form vs quadrature); C orthogonality identity (D_1 a_chi, Y_s) = L(s,chi)/s;
D zeros and the constant C_chi; E finite-lambda Burnol bound L*||P_R 1||^2 (R = zeros up to Tz).
Usage: python3 nbdir_burnol.py <m4|5> <Tz> <parts e.g. ABCDE>
"""
import sys, time, json
import mpmath as mp
import numpy as np
mp.mp.dps = 20
CHI = {"m4": (4, [0, 1, 0, -1], 1), "5": (5, [0, 1, -1, -1, 1], 0)}
name = sys.argv[1]; Tz = float(sys.argv[2]); parts = sys.argv[3]
q, chil, a = CHI[name]
pi = mp.pi
def Lfun(s): return mp.dirichlet(s, chil)
def gam_a(s): return (q / pi)**(s - 0.5) * mp.gamma((s + a) / 2) / mp.gamma((1 - s + a) / 2)
def gpm(s):  # Mellin of 2cos(2pi x) (a=0) / 2sin(2pi x) (a=1)
    return 2 * (2 * pi)**(-s) * mp.gamma(s) * (mp.cos(pi * s / 2) if a == 0 else mp.sin(pi * s / 2))
def U(w): return w / (1 - w) * gpm(w)
def g(x): return mp.sin(2 * pi * x) / (pi * x) if a == 0 else (1 - mp.cos(2 * pi * x)) / (pi * x)
def J01(w, x):  # int_0^1 u^{w-2} c_a(2 pi x u) du by series
    y = 2 * pi * x
    if a == 0: return mp.nsum(lambda j: (-1)**j * y**(2*j+1) / (mp.factorial(2*j+1) * (w + 2*j)), [0, mp.inf])
    return mp.nsum(lambda j: (-1)**(j+1) * y**(2*j) / (mp.factorial(2*j) * (w + 2*j - 1)), [1, mp.inf])
def Etail(w, Y, sgn):  # int_Y^inf y^{w-2} e^{sgn i y} dy
    p = -sgn * 1j
    return p**(1 - w) * mp.gammainc(w - 1, p * Y)
def Ttail(w, Y):
    if a == 0: return (Etail(w, Y, 1) - Etail(w, Y, -1)) / (2j)
    return Y**(w - 1) / (1 - w) - (Etail(w, Y, 1) + Etail(w, Y, -1)) / 2
def phi(w, x):  # Lemma B-chi, k = 0
    x = mp.mpf(x)
    if x <= 1.5: return U(w) * x**(-w) + g(x) - w / (pi * x) * J01(w, x)
    return g(x) + w / (pi * x) * (2 * pi * x)**(1 - w) * Ttail(w, 2 * pi * x)
out = {}
if "A" in parts:
    s = mp.mpc(0.3, 2.7)
    closed = pi**(0.5 - s) * mp.gamma(s / 2) / mp.gamma((1 - s) / 2) if a == 0 else pi**(0.5 - s) * mp.gamma((1 + s) / 2) / mp.gamma(1 - s / 2)
    e1 = abs(gpm(s) - closed)
    e2 = abs(Lfun(1 - s) / Lfun(s) - gam_a(s))
    e3 = abs(gam_a(s) - q**(s - 0.5) * gpm(s))
    print(f"A: |gamma_pm closed forms| {mp.nstr(e1,3)}  |L(1-s)/L(s) - gamma_a(s)| {mp.nstr(e2,3)}  |gamma_a - q^(s-1/2) gamma_pm| {mp.nstr(e3,3)}", flush=True)
if "B" in parts:
    w = mp.mpc(0.4, 5)
    for x in [0.3, 1.2, 3.0, 7.5]:
        if a == 0: Iq = mp.quadosc(lambda u: u**(w - 2) * mp.sin(2*pi*x*u), [1, mp.inf], omega=2*pi*x)
        else: Iq = 1 / (1 - w) - mp.quadosc(lambda u: u**(w - 2) * mp.cos(2*pi*x*u), [1, mp.inf], omega=2*pi*x)
        direct = g(x) + w / (pi * x) * Iq
        print(f"B: x={x} phi closed {mp.nstr(phi(w, x), 12)} direct {mp.nstr(direct, 12)} |U x^-w| {mp.nstr(abs(U(w)*mp.mpf(x)**(-w)),6)}", flush=True)
if "C" in parts:  # (D_1 a_chi, Y_s) = int_1^inf A(t)/t * conj(q^{-1/2} phi_s(t/q)) dt = L(s,chi)/s
    Aper = np.cumsum(chil)
    for s in [mp.mpc(0.5, 3.0)]:
        T = 600; tot = mp.mpf(0); xs, ws = np.polynomial.legendre.leggauss(12)
        for n in range(1, T):
            An = Aper[n % q]
            if An == 0: continue
            for xi, wi in zip(xs, ws):
                t = n + 0.5 + 0.5 * xi
                tot += 0.5 * wi * An / t * mp.conj(phi(s, t / q))
        tot /= mp.sqrt(q)
        # tail beyond T: A ~ c0, conj phi ~ g(t/q): c0 q^{-1/2} int_T^inf g(t/q)/t dt ~ c0 q^{-1/2} * q/(pi T) (a=1; a=0: oscillatory ~0)
        # tail: A(t) * g(t/q) ~ q-periodic numerator / (pi t/q): int_T^inf A(t) t^-1 q^-1/2 g(t/q) dt ~ q^-1/2 (q/pi) mA / T
        mA = mp.quad(lambda t: Aper[int(mp.floor(t)) % q] * (mp.sin(2*pi*t/q) if a == 0 else 1 - mp.cos(2*pi*t/q)), mp.linspace(0, q, q + 1)) / q
        tail = mp.sqrt(q) / pi * mA / T
        if a == 1: tail += float(np.mean(Aper)) * mp.conj(s / (pi * (1 - s))) * q / T / mp.sqrt(q)  # phi ~ g + kap/x (a = 1)
        print(f"C: s={s} int={mp.nstr(tot,10)} +tail~{mp.nstr(tail,4)} -> {mp.nstr(tot+tail,10)} ; L(s)/s={mp.nstr(Lfun(s)/s,10)}", flush=True)
zeros = []
if "D" in parts or "E" in parts:
    def Z(t):
        s = mp.mpc(0.5, t)
        ph = mp.im(mp.loggamma((s + a) / 2)) + (t / 2) * mp.log(q / pi)
        return mp.re(mp.exp(1j * ph) * Lfun(s))
    ts = np.arange(0.05, Tz, 0.05); prev = Z(ts[0])
    for i in range(1, len(ts)):
        cur = Z(ts[i])
        if prev * cur < 0: zeros.append(float(mp.findroot(Z, (ts[i-1], ts[i]), solver="anderson")))
        prev = cur
    Cchi = 2 * mp.re(mp.diff(Lfun, 1) / Lfun(1)) + mp.log(q / pi) + mp.digamma((1 + a) / 2)
    part = sum(2 / (0.25 + g_**2) for g_ in zeros)
    Ttop = Tz
    tail = (1 / pi) * (mp.log(q * Ttop / (2 * pi)) + 1) / Ttop  # N(T) density tail, both signs
    print(f"D: {len(zeros)} zeros in (0,{Tz}); first {zeros[:4]}; sum_(|g|<{Tz}) 1/|rho|^2 = {mp.nstr(part,10)}; +N(T) tail ~{mp.nstr(tail,6)} = {mp.nstr(part+tail,8)}; C_chi = {mp.nstr(Cchi,12)}", flush=True)
    out["zeros"] = zeros; out["Cchi"] = float(Cchi); out["partial"] = float(part)
if "E" in parts:
    rhos = [mp.mpc(0.5, g_) for g_ in zeros] + [mp.mpc(0.5, -g_) for g_ in zeros]
    nR = len(rhos); Us = np.array([complex(U(r)) for r in rhos]); rc = np.array([complex(r) for r in rhos])
    gl_x, gl_w = np.polynomial.legendre.leggauss(10)
    def panels(lo, hi, w):
        k = int(np.ceil((hi - lo) / w)); e = np.linspace(lo, hi, k + 1)
        nd = ((e[:-1, None] + e[1:, None]) / 2 + (e[1:, None] - e[:-1, None]) / 2 * gl_x[None, :]).ravel()
        ww = ((e[1:, None] - e[:-1, None]) / 2 * gl_w[None, :]).ravel()
        return nd, ww
    X = 40.0; nodes, wts = panels(1.0, X, 0.25)
    t0 = time.time()
    Dx = np.array([[complex(phi(r, x) - g(x)) for x in nodes] for r in rhos])
    gx = np.array([float(g(x)) for x in nodes])
    print(f"E: delta on [1,{X}] sampled ({len(nodes)} nodes x {nR}) in {time.time()-t0:.0f}s", flush=True)
    def Cint(om): return mp.cos(om) - om * (pi / 2 - mp.si(om))
    g2 = float((1 / pi**2) * (1.5 - 2 * Cint(2 * pi) + 0.5 * Cint(4 * pi))) if a == 1 else float((1 / pi**2) * 0.5 * (1 - Cint(4 * pi)))
    # (0,1]: phi = U x^-rho + r(x), r(x) = g(x) - (rho/(pi x)) J01(rho, x), numpy series; v = log(1/x) in [0, Vc]
    Vc = 60.0; vn, vw = panels(0.0, Vc, 0.05); xv = np.exp(-vn); y = 2 * np.pi * xv
    gv = np.sin(y) / (np.pi * xv) if a == 0 else (1 - np.cos(y)) / (np.pi * xv)
    from math import factorial
    R = np.zeros((nR, len(xv)), complex)
    for i, r in enumerate(rc):
        Jv = np.zeros(len(xv), complex)
        for j in range(0, 45):
            if a == 0: Jv += (-1)**j * y**(2*j+1) / (factorial(2*j+1) * (r + 2*j))
            elif j >= 1: Jv += (-1)**(j+1) * y**(2*j) / (factorial(2*j) * (r + 2*j - 1))
        R[i] = gv - r / (np.pi * xv) * Jv
    # spot check r vs mpmath at x = 0.7
    chk = abs(complex(phi(rhos[0], 0.7) - U(rhos[0]) * mp.mpf(0.7)**(-rhos[0])) - (np.interp(0.7, xv[::-1], R[0].real[::-1]) + 1j*np.interp(0.7, xv[::-1], R[0].imag[::-1])))
    print(f"E: r-series vs mpmath at x=0.7 (interp): {chk:.2e}", flush=True)
    Mn = np.exp(np.outer(rc, vn))          # x^-rho = e^{v rho}
    gd = (Dx.conj() * gx * wts).sum(1); dg = (Dx * gx * wts).sum(1); dd = (Dx * wts) @ Dx.conj().T
    if a == 1:  # [X, inf) tails: delta_rho ~ rho/(pi (1-rho) x), g ~ 1/(pi x) on average
        kap = rc / (np.pi * (1 - rc))
        dd = dd + np.outer(kap, kap.conj()) / X; gd = gd + kap.conj() / (np.pi * X); dg = dg + kap / (np.pi * X)
    tail1 = g2 + gd[None, :] + dg[:, None] + dd
    res = []
    gam = rc.imag
    for Lval in [4, 8, 12, 16, 24, 40, 100, 1000, 10**5]:
        V = Lval + np.log(q)
        wq = vw * xv * (vn <= V)            # only x in (lam/q, 1]  (Vc = 60 caps the near region)
        cross = ((Us[:, None] * Mn) * wq) @ R.conj().T + (R * wq) @ (Us[:, None] * Mn).conj().T + (R * wq) @ R.conj().T
        bnear = ((gv * wq)[None, :] * (Us[:, None] * Mn + R).conj()).sum(1)
        b = bnear + g2 + gd
        dlt = gam[:, None] - gam[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            I = np.where(np.abs(dlt) < 1e-12, V, (np.exp(1j * dlt * V) - 1) / (1j * dlt))
        main = Us[:, None] * Us.conj()[None, :] * I
        Gm = main + cross + tail1               # Gm[i,j] = int phi_i conj(phi_j) = (Y_i, Y_j)
        c = np.linalg.solve(Gm.T, b); proj = float(np.real(np.vdot(b, c)))
        asym = float(np.sum(1 / np.abs(rc)**2))
        res.append(dict(L=Lval, LP=Lval * proj, asym=asym))
        print(f"E: L={Lval:7.0f} L*||P_R 1||^2 = {Lval*proj:.6f}  (R: {nR} zeros, |g|<{Tz}; limit sum_R 1/|rho|^2 = {asym:.6f}); |b0|^2={abs(b[0])**2:.6f} vs 1/|rho0|^2={1/abs(rc[0])**2:.6f}", flush=True)
    out["E"] = res
json.dump(out, open(f"/Users/adriansager/playground/riemann/.scratch/nbdir/burnol_{name}_{Tz}_{parts}.json", "w"))
