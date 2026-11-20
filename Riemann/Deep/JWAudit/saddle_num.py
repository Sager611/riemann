# JWAudit: high-precision falsification test of the explicit complex saddle bound (HollandK DERIVATION Sec. 2).
# Mt(z) = int_0^oo phi(u) u^{2z} du; claim: Mt = 2pi^2 e^{Psi(z)} sqrt(2pi/(-f''(u0))) (1+delta), |delta| <= deltaA,
# and sup_{|z-x0|<=x0/10} |G(z)-G(x0)| <= etaG, G = log Mt - const - Psi.  We compute the k=1 integral numerically
# (real axis for real z / small x0; horizontal ray through the complex saddle u0(z) otherwise) at mp precision.
# Run: uv run --no-project --with mpmath python3 saddle_num.py
import mpmath as mp, sys, time
pi = None
def setp(d):
    global pi
    mp.mp.dps = d; pi = mp.pi
setp(60)
def Z(u): return u * (pi * mp.exp(2 * u) - mp.mpf(9) / 4)
def Zp(u): return pi * mp.exp(2 * u) * (1 + 2 * u) - mp.mpf(9) / 4
def u0of(z, g):
    u = mp.mpc(g)
    for _ in range(400):
        du = (Z(u) - z) / Zp(u); u -= du
        if abs(du) < mp.mpf(10) ** (-mp.mp.dps + 8): break
    return u
def f(u, z): return 2 * z * mp.log(u) + mp.mpf(9) * u / 2 - pi * mp.exp(2 * u)
def f2(u, z): return -2 * z / u ** 2 - 4 * pi * mp.exp(2 * u)
def f3(u, z): return 4 * z / u ** 3 - 8 * pi * mp.exp(2 * u)
def f4(u, z): return -12 * z / u ** 4 - 16 * pi * mp.exp(2 * u)
def amp(u): return 1 - 3 / (2 * pi) * mp.exp(-2 * u)
def J_ray(z, u0):
    # int over the horizontal ray u = u0 + s of e^{f(u)-f(u0)} a(u) ds, s = sigma/sqrt(|kappa|), sigma in [-40, 40]
    kap = -f2(u0, z); r = mp.sqrt(abs(kap)); F0 = f(u0, z)
    g = lambda sg: mp.exp(f(u0 + sg / r, z) - F0) * amp(u0 + sg / r)
    I = mp.quad(g, [-40, -20, -10, -4, 0, 4, 10, 20, 40]) / r
    tail = max(mp.re(f(u0 + 40 / r, z) - F0), mp.re(f(u0 - 40 / r, z) - F0))
    return I, kap, tail
def J_real(x, nk=1):
    # independent: real-axis integral of the first nk theta terms of phi, divided by e^{f(u0(x))} (k=1 normalization)
    u0 = mp.re(u0of(x, mp.log(x) / 2)); kap = -f2(u0, x); w = 60 / mp.sqrt(kap); F0 = f(u0, x)
    def ph(u):
        s = 0
        for k in range(1, nk + 1):
            s += (2 * pi**2 * k**4 * mp.exp(mp.mpf(9) * u / 2) - 3 * pi * k**2 * mp.exp(mp.mpf(5) * u / 2)) * mp.exp(-pi * k * k * mp.exp(2 * u))
        return s
    g = lambda u: ph(u) * mp.exp(2 * x * mp.log(u) - F0)
    pts = [u0 - w, u0 - w / 4, u0, u0 + w / 4, u0 + w]
    I = mp.quad(g, pts) / (2 * pi**2)
    return I, u0, kap
out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
t0 = time.time()
# ---- 1. validate the ray quadrature against the real-axis integral (small x0, real and complex z) ----
setp(50)
for x in [mp.mpf(2000), mp.mpf(10) ** 6]:
    Ir, u0, kap = J_real(x, nk=3)
    Iy, kap2, tl = J_ray(x, mp.mpc(u0))
    P("validate real x=%s: real-axis(k<=3) %s  ray(k=1) %s  rel diff %s" % (mp.nstr(x, 5), mp.nstr(Ir, 20), mp.nstr(mp.re(Iy), 20), mp.nstr(abs(Ir / Iy - 1), 3)))
# complex z at x0 = 2000 via the real axis (oscillatory but resolvable) vs ray through the complex saddle
x0 = mp.mpf(2000); a0 = mp.re(u0of(x0, 2))
for th in [mp.pi / 2, 3 * mp.pi / 4]:
    z = x0 * (1 + mp.mpf('0.1') * mp.expjpi(th / mp.pi))
    u0 = u0of(z, a0); Iy, kap, tl = J_ray(z, u0)
    ur = mp.re(u0of(mp.re(z), a0)); kr = -f2(ur, mp.re(z)); w = 60 / mp.sqrt(kr)
    g = lambda u: (2 * pi**2 * mp.exp(mp.mpf(9) * u / 2) - 3 * pi * mp.exp(mp.mpf(5) * u / 2)) * mp.exp(-pi * mp.exp(2 * u)) * mp.exp(2 * z * mp.log(u) - f(u0, z))
    Ir = mp.quad(g, mp.linspace(ur - w, ur + w, 41)) / (2 * pi**2)
    P("validate complex z=2000(1+0.1e^{i%s}): real-axis %s  ray %s  rel diff %s" % (mp.nstr(th, 4), mp.nstr(Ir, 15), mp.nstr(Iy, 15), mp.nstr(abs(Ir / Iy - 1), 3)))
# ---- 2. the claimed constants vs the true values on the disk edge ----
CL = {24: ('1.16368e-5', '0.0552241', '0.0551866'), 80: ('7.90418e-18', '0.0533654', '0.0534309'), 90: ('5.32143e-20', '0.0532866', '0.0533475')}
cases = [(24, None), (80, None), (90, None), (90, mp.mpf(10) ** 45), (90, mp.mpf(10) ** 50), (90, mp.mpf(10) ** 60)]
NTH = int(sys.argv[1]) if len(sys.argv) > 1 else 16
for Ls, xv in cases:
    setp(70)
    if xv is None: x0 = Z(mp.mpf(Ls) / 2)
    else: x0 = xv
    setp(int(70 + 2.2 * mp.log10(x0)))
    x0 = Z(mp.mpf(Ls) / 2) if xv is None else mp.mpf(xv)
    a0 = mp.re(u0of(x0, mp.log(x0) / 2))
    dA, eG, bet = [mp.mpf(v) for v in CL[Ls]]
    I0, k0, tl0 = J_ray(x0, mp.mpc(a0)); d0 = I0 / mp.sqrt(2 * pi / k0) - 1
    # first-order Laplace prediction of delta at x0 (sanity of the numerics)
    kp = k0; F3 = f3(a0, x0); F4 = f4(a0, x0); A = amp(a0); A1 = 3 / pi * mp.exp(-2 * a0); A2 = -6 / pi * mp.exp(-2 * a0)
    pred = A * (1 + (F4 / (8 * kp) + 5 * F3**2 / (24 * kp**2)) / kp) + (A2 / 2 + A1 * F3 / (2 * kp)) / kp - 1
    G0 = mp.log(I0); mxd = abs(d0); mxG = mp.mpf(0); mxb = mp.mpf(0); argmx = None; mxtail = tl0
    for j in range(NTH):
        th = 2 * mp.pi * j / NTH
        z = x0 * (1 + mp.mpf('0.1') * mp.expjpi(2 * mp.mpf(j) / NTH))
        u0 = u0of(z, a0)
        I, kap, tl = J_ray(z, u0)
        dl = I / mp.sqrt(2 * pi / kap) - 1
        dG = abs(mp.log(I) - G0)
        mxd = max(mxd, abs(dl)); mxb = max(mxb, abs(u0 - a0)); mxtail = max(mxtail, tl)
        if dG > mxG: mxG = dG; argmx = j
    P("L*=%d x0=%s a0=%s: delta(x0)=%s (1st-order Laplace %s); max|delta| on edge %s <= deltaA %s: %s" % (
        Ls, mp.nstr(x0, 6), mp.nstr(a0, 8), mp.nstr(d0, 4), mp.nstr(pred, 4), mp.nstr(mxd, 4), mp.nstr(dA, 4), mxd <= dA))
    P("   max|G(z)-G(x0)| on edge (%d pts) %s at theta=2pi*%d/%d <= etaG %s: %s (margin %s); max|u0(z)-a0| %s <= beta %s: %s; ray tail exponent <= %s" % (
        NTH, mp.nstr(mxG, 8), argmx, NTH, mp.nstr(eG, 6), mxG <= eG, mp.nstr(eG - mxG, 3), mp.nstr(mxb, 6), mp.nstr(bet, 6), mxb <= bet, mp.nstr(mxtail, 4)))
P("time %.1f s" % (time.time() - t0))
open(__file__.replace('.py', '.log'), 'w').write("\n".join(out) + "\n")
