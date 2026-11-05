# fh2lib.py -- rigorous (ball arithmetic) primitives for the diagonal Fasi-Higham proof
# zeta-plane: z_orig = 2 nu zeta, nu = m + 1/2.
import flint
from flint import arb, acb
flint.ctx.prec = 128
PI = arb.pi()
I = acb(0, 1)
def ab(x):
    return x if isinstance(x, arb) else arb(x)
def ball(lo, hi):
    lo, hi = ab(lo), ab(hi)
    return (lo + hi) / 2 + arb(0, ((hi - lo) / 2).upper())   # arb(mid, rad)
def cball(xlo, xhi, ylo, yhi):
    return acb(ball(xlo, xhi), ball(ylo, yhi))
def polar(rlo, rhi, alo, ahi):
    # box containing {rho e^{i a}: rho in [rlo,rhi], a in [alo,ahi]}
    rr = ball(rlo, rhi); aa = ball(alo, ahi)
    return acb(rr * aa.cos(), rr * aa.sin())
def csqrt(w):
    # principal sqrt of a complex ball without the Im-cancellation of the textbook formula;
    # valid when w stays away from the closed negative real axis (|w| + Re w > 0)
    x = w.real; mod = abs(w)
    sr = ((mod + x) / 2).sqrt()
    if not (sr > 0):
        return w.sqrt()
    return acb(sr, w.imag / (2 * sr))
def S(z):            # s = sqrt(1+z^2), principal
    return csqrt(1 + z * z)
def Psie(z):
    s = S(z); return s - (1 + s).log()
def Psi(z): return z + Psie(z)
def eta(z): return Psie(z) + z.log()
def A(a, z): return z + a * Psie(z)
def DkPsie(z, k):
    # factored forms (less dependency): D = z d/dz
    s = S(z); z2 = z * z
    if k == 1: return z2 / (s + 1)
    if k == 2: return z2 / s
    if k == 3: return z2 * (2 + z2) / s**3
    if k == 4: return z2 * (1 / s + 3 / s**5)
def DkA(a, z, k): return z + a * DkPsie(z, k)
def DkP(z, k):       # P = (1/4) log(1+z^2),  u = 1/(1+z^2), 1-u = z^2 u
    u = 1 / (1 + z * z); v = z * z * u
    if k == 1: return v / 2
    if k == 2: return u * v
    if k == 3: return -2 * u * v * (1 - 2 * u)
    if k == 4: return 4 * u * (1 - 7 * u + 12 * u**2 - 6 * u**3)
def fvar(z):
    # |dU1/dzeta| = |3-15p^2|/24 * |z| |p|^3,  p = 1/s
    s = S(z); p = 1 / s
    return abs(3 - 15 * p * p) * abs(z) * abs(p)**3 / 24
def up(x):          # upper endpoint as arb
    return arb(x.upper()) if isinstance(x, arb) else arb(abs(x).upper())
def angle_chunks(alo, ahi, w=0.04):
    alo, ahi = ab(alo), ab(ahi)
    n = max(1, int(float(((ahi - alo) / w).upper())) + 1)
    return [(alo + (ahi - alo) * j / n, alo + (ahi - alo) * (j + 1) / n) for j in range(n)]
def supf_polar(r0, r1, chunks):
    m = arb(0)
    for (c0, c1) in chunks:
        v = up(fvar(polar(r0, r1, c0, c1)))
        if v > m: m = v
    return m
def VI_sector(T, alo, ahi, K=200):
    # int_0^T sup_{a in [alo,ahi]} f(rho e^{ia}) d rho  (I-variation along rays from 0)
    T = ab(T); tot = arb(0); ch = angle_chunks(alo, ahi)
    for k in range(K):
        r0 = T * k / K; r1 = T * (k + 1) / K
        tot += (r1 - r0) * supf_polar(r0, r1, ch)
    return arb(tot.upper())
def VK_sector(T0, alo, ahi, K=200, Smax=8):
    # int_{T0}^inf sup_a f(rho e^{ia}) d rho  (K-variation along rays to infinity); needs |a| < pi/2
    T0 = ab(T0); Smax = ab(Smax); tot = arb(0); ch = angle_chunks(alo, ahi)
    for k in range(K):
        r0 = T0 + (Smax - T0) * k / K; r1 = T0 + (Smax - T0) * (k + 1) / K
        tot += (r1 - r0) * supf_polar(r0, r1, ch)
    # tail rho >= Smax: |1+rho^2 e^{2ia}| >= rho^2-1 ; f <= (3+15/(rho^2-1)) rho (rho^2-1)^{-3/2}/24
    q = Smax * Smax - 1
    tot += (3 + 15 / q) / (24 * q.sqrt())
    return arb(tot.upper())
def VK_horiz(x0, x1, y, K=100):
    # int_{x0}^{x1} f(u + i y) du
    x0, x1, y = ab(x0), ab(x1), ab(y); tot = arb(0)
    for k in range(K):
        u0 = x0 + (x1 - x0) * k / K; u1 = x0 + (x1 - x0) * (k + 1) / K
        z = acb(ball(u0, u1), y)
        tot += (u1 - u0) * up(fvar(z))
    return arb(tot.upper())
def VI_imag(y1, K=200):
    # I-variation along the imaginary axis from 0 to i*y1 (y1<1): p real, U1 monotone for p>=1:
    # exact: |U1(p(iy1)) - U1(1)|, p = (1-y^2)^{-1/2}
    y1 = ab(y1); p = 1 / (1 - y1 * y1).sqrt()
    return arb(abs((3 * p - 5 * p**3) / 24 - arb(-2) / 24).upper())
def epsLG(V, nu):
    return arb((((2 * V) / nu).exp() - 1).upper())

def sup_polar(func, rlo, rhi, plo, phi, w=0.02):
    # upper bound of |func(z)| over the polar box, angle chunks of width <= w
    m = arb(0)
    for (c0, c1) in angle_chunks(plo, phi, w):
        v = up(abs(func(polar(rlo, rhi, c0, c1))))
        if v > m: m = v
    return m
PHIK = arb('1.2')
def VK_box(rlo, rhi, plo, phi):
    # sup of K-variation over points rho e^{i a}, rho in [rlo,rhi], a in [plo,phi] (0<=plo<=phi<=pi/2)
    rlo, rhi, plo, phi = ab(rlo), ab(rhi), ab(plo), ab(phi)
    V = arb(0)
    if plo < PHIK:
        V = VK_sector(rlo, plo, phi if phi < PHIK else PHIK)
    if phi > PHIK:
        p0 = plo if plo > PHIK else PHIK
        # y-range of the part with angle >= p0 ; x in [0, ...]
        ylo = rlo * (p0.sin() if p0.sin() < phi.sin() else phi.sin())
        yhi = rhi
        xend = yhi * PHIK.cos() / PHIK.sin()
        K = 80; tot = arb(0)
        ych = [(ylo + (yhi - ylo) * j / 6, ylo + (yhi - ylo) * (j + 1) / 6) for j in range(6)]
        for k in range(K):
            u0 = xend * k / K; u1 = xend * (k + 1) / K
            mm = arb(0)
            for (y0, y1) in ych:
                v = up(fvar(acb(ball(u0, u1), ball(y0, y1))))
                if v > mm: mm = v
            tot += (u1 - u0) * mm
        tot += VK_sector(ylo / PHIK.sin(), PHIK, PHIK)
        if tot > V: V = tot
    return arb(V.upper())
def VI_box(rlo, rhi, plo, phi):
    return VI_sector(rhi, plo, phi)
def lg_box(rlo, rhi, plo, phi, nu):
    """upper bounds (epsI, epsK, lamI, lamK) for |h_I|,|h_K| (and |log(1+h)|) on the polar box (angles in [0,pi/2])"""
    eI = epsLG(VI_box(rlo, rhi, plo, phi), nu)
    eK = epsLG(VK_box(rlo, rhi, plo, phi), nu)
    if not (eI < 1 and eK < 1): return None
    lI = arb((-(1 - eI).log()).upper()); lK = arb((-(1 - eK).log()).upper())
    return eI, eK, lI, lK
