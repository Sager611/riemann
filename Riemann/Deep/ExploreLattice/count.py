# count.py (explore48): exact zero counting for E*(z,s) (argument principle, no Turing/S(T) assumption).
# N(T) = #{zeros rho of E*(z,.) (any Re), 0 < Im rho <= T}
#      = 1 + (1/pi)[Im lgamma(sR+iT) - T log pi + T log y + Arg Etilde(sR+iT) + D_hor arg G + pi(1/2 - sR)/2] - N_real/2,
# Etilde = y^{-s} E(z,s) (|Etilde - 1| < 1 on Re s = sR, checked), G = e^{-i pi s/2} E*, D_hor = continuous change of arg G
# from sR+iT to 1/2+iT; N_real = # real zeros of E* (0 or 2: real pair in (0,1) iff y > 7.0555 for these families).
# N0(T) = # sign changes of F on (0,T]; off-line count = N(T) - N0(T).
import math, cmath, sys
sys.path.insert(0, '.')
from flint import acb, arb, ctx
from epstein_cs import Estar, F_line, theta_a

def logEt(s, x, y, G=None):
    G = Estar(s, x, y) if G is None else G
    old = ctx.prec; ctx.prec = 80
    lg = acb(s.real, s.imag).lgamma(); ctx.prec = old
    lgc = complex(float(lg.real.mid()), float(lg.imag.mid()))
    return -s * math.log(y) + s * math.log(math.pi) + 1j * math.pi * s / 2 + cmath.log(G) - lgc

def Etilde_real_excess(x, y, sR):
    # Etilde(sR) - 1 bound: sum over (m,n) != (0,+-1) of (y/|mz+n|^2)^sR / 2 ... computed by direct sum (real, positive)
    tot = 0.0
    M = 40
    for m in range(0, M + 1):
        for n in range(-200, 201):
            if m == 0 and n <= 0: continue
            if m == 0 and n == 1: continue
            q = (m * x + n) ** 2 + (m * y) ** 2
            tot += (1.0 / q) ** sR  # (y^s/|mz+n|^{2s}) / y^s summed over half-plane representatives (1/2 sum' -> one of +-)
    return tot

def NE(T, x, y, sR=3.0, maxjump=0.6, Nreal=None):
    if Nreal is None: Nreal = 2 if y > 7.0555 else 0
    # horizontal arg tracking of G from sR+iT to 1/2+iT
    def g(sig): return Estar(complex(sig, T), x, y)
    pts = [(sR, g(sR))]
    total = 0.0
    stack = [(sR, pts[0][1], 0.5, g(0.5))]
    # adaptive subdivision
    def seg(a, ga, b, gb, depth=0):
        d = cmath.phase(gb / ga)
        if abs(d) < maxjump or depth > 18: return d
        m = 0.5 * (a + b); gm = g(m)
        return seg(a, ga, m, gm, depth + 1) + seg(m, gm, b, gb, depth + 1)
    Ga = pts[0][1]; Gb = g(0.5)
    dh = seg(sR, Ga, 0.5, Gb)
    lEt = logEt(complex(sR, T), x, y, G=Ga)
    argEt = math.remainder(lEt.imag, 2 * math.pi)
    old = ctx.prec; ctx.prec = 80
    lg = acb(sR, T).lgamma(); ctx.prec = old
    val = 1 + (float(lg.imag.mid()) - T * math.log(math.pi) + T * math.log(y) + argEt + dh + math.pi * (0.5 - sR) / 2) / math.pi - Nreal / 2
    return val, argEt

def N0_scan(T, x, y, t_start=0.05, step_scale=1.0):
    # sign changes of F on (0, T] with step pi/(8 log(t y/pi)); returns count and list of sign-change positions
    t = t_start; Fp = F_line(t, x, y); n = 0; pos = []
    while t < T:
        h = step_scale * math.pi / (8 * math.log(max(t * y / math.pi, 3.0)))
        tn = min(t + h, T)
        Fn = F_line(tn, x, y)
        if (Fn > 0) != (Fp > 0): n += 1; pos.append(tn)
        t, Fp = tn, Fn
    return n, pos
