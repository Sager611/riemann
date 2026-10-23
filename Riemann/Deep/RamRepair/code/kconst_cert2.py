# kconst_cert2.py -- rigorous enclosures of k1 = H_c(0) and k0 = H_c'(0) + gamma H_c(0) (Padhi App. A, end of the proof of Thm 6.5).
#   H_c(0) = (1/2 pi i) int_P g,  H_c'(0) = (1/2 pi i) int_P g D,   g(s) = A(0; -s, s) e^{eta^2 s^2} / (zeta(1-s) zeta(1+s) s^4),
#   D(s) = d_w log A + eta^2 s - zeta'/zeta(1+s) - 2/s,  P: -c' - iT2 -> k1 - iT2 -> k1 + iT2 -> -c' + iT2, k1 = 0.38, T2 = 400.
# Vertical part: (1/pi) Re int_0^{T2} g(k1 + it) dt (conjugate symmetry).  [0, 80]: acb.integral (rigorous; g is meromorphic and every
# pole inside an evaluation ball makes the ball evaluation non-finite) with the Euler product truncated at P0 = 400.
# Fixes against kconst_cert.py: (i) int |g|, int |g D| by ball subdivision with dyadic steps (old loop gave nan);
# (ii) tau = sum_{p > P0} sup |A_p - 1| summed exactly to 1e6 plus an explicit tail; (iii) |d_w log A| bounded by an exact prime sum
# (the old constant 0.3 was wrong: p = 2 alone contributes 1.07); (iv) horizontal pieces with a computed sup |A| and a Cauchy bound for zeta'/zeta.
import math, json, time, sys
from flint import arb, acb, acb_series, ctx
ctx.prec = 64
eta, kap1, T2, cp = arb('0.02'), arb('0.38'), arb(400), arb('0.01')
P0 = 400
def primes(n):
    sv = bytearray([1])*(n + 1); sv[0:2] = b'\x00\x00'
    for i in range(2, int(n**0.5) + 1):
        if sv[i]: sv[i*i::i] = bytearray(len(sv[i*i::i]))
    return [i for i in range(n + 1) if sv[i]]
PALL = primes(10**6)
PR = [p for p in PALL if p <= P0]; LP = [arb(p).log() for p in PR]
def gD(s):
    """(g, g*D) at complex s (acb ball), Euler product truncated at P0."""
    A = acb(1); dl = acb(0)
    for p, lp in zip(PR, LP):
        x1 = ((s - 1)*lp).exp(); x2 = ((-1 - s)*lp).exp(); zz = acb(1)/p
        A *= (1 - x1 - x2 + zz)*(1 - zz)/((1 - x1)*(1 - x2))
        dl += -lp*x2*(1/(1 - x2) - 1/(1 - x1 - x2 + zz))
    zm = (1 - s).zeta()
    ser = acb_series([1 + s, 1], prec=2).zeta(); zp, dzp = ser[0], ser[1]
    g = A*(eta*eta*s*s).exp()/(zm*zp*s**4)
    D = dl + eta*eta*s - dzp/zp - 2/s
    return g, g*D
t0 = time.time()
# ---- (i) int_0^80 |g|, |g D| : balls t in [j h, (j+1) h], h = 2^-5, halved where the evaluation is not finite ----
def absint(lo, hi, depth=0):
    """rigorous upper bounds of int_lo^hi |g| and int_lo^hi |g D| on Re s = k1: ball evaluation on [lo, hi] (sup bound times length);
    the interval is bisected while the ball bound exceeds 1.25 x the value at the midpoint (and the piece is not < 1e-9) (tightness only; every accepted piece is rigorous)."""
    tb = arb((lo + hi)/2, (hi - lo)/2)
    g, gd = gD(acb(kap1, tb))
    ag, agd = abs(g), abs(gd)
    if ag.is_finite() and agd.is_finite():
        gm, gdm = gD(acb(kap1, arb((lo + hi)/2)))
        w_ = hi - lo
        ok = ((ag.upper() <= 1.25*abs(gm).upper() + 1e-12) and (agd.upper() <= 1.25*abs(gdm).upper() + 1e-12)) \
            or (ag.upper()*w_ <= 1e-9 and agd.upper()*w_ <= 1e-9)
        if ok or depth >= 16:
            w = arb(hi - lo); return arb(ag.upper())*w, arb(agd.upper())*w
    elif depth >= 16: raise RuntimeError('not finite at %r' % lo)
    m = (lo + hi)/2; a1, b1 = absint(lo, m, depth + 1); a2, b2 = absint(m, hi, depth + 1); return a1 + a2, b1 + b2
h = 2.0**-5; absg = arb(0); absgD = arb(0)
for j in range(int(80/h)):
    a_, b_ = absint(j*h, (j + 1)*h); absg += a_; absgD += b_
print('int_0^80 |g| <=', absg.str(8), ' |gD| <=', absgD.str(8), '%.0fs' % (time.time() - t0), flush=True)
# ---- (ii) tau = sum_{p > P0} sup_{Re s = k1} |A_p - 1|,  A_p - 1 = (x1 z + x2 z - z^2 - x1 x2)/((1 - x1)(1 - x2)) ----
# and dtail = sum_{p > P0} sup |d_w log A_p|, d_w log A_p = -log p x2 (z - x1)/((1 - x2)(1 - x1 - x2 + z)), |1 + z - x1 - x2| >= 1 + z - r1 - r2 (z = 1/p real)
def ap_terms(p, e1, e2):
    pa = arb(p); r1, r2, r3 = pa**(-e1), pa**(-e2), 1/pa
    am1 = (r1*r3 + r2*r3 + r3*r3 + r1*r2)/((1 - r1)*(1 - r2))
    den = 1 + r3 - r1 - r2
    dla = pa.log()*r2*(r1 + r3)/((1 - r2)*den) if den > 0 else arb('inf')
    return am1, dla
E1, E2 = 1 - kap1, 1 + kap1
tau = arb(0); dtail = arb(0); dlall = arb(0)
for p in PALL:
    am1, dla = ap_terms(p, E1, E2)
    dlall += dla
    if p > P0: tau += am1; dtail += dla
# p > 1e6: sum_{p > P} p^{-m} <= 1.25506 m P^{1-m}/((m-1) log P) (partial summation, pi(t) < 1.25506 t/log t),
# |A_p - 1| <= 4 p^{-m}/((1 - r1)(1 - r2)) with m = 1 + min(E1, 1) = 1.62; |d_w log A_p| <= 2 log p p^{-2}/((1-r2)(1 + r3 - r1 - r2)),
# sum_{p > P} log p p^{-2} <= 2 * 1.01624/P (theta(t) <= 1.01624 t, Dusart).
Pa = arb(10**6); m = 1 + E1
r1P, r2P = Pa**(-E1), Pa**(-E2)
tl = 4*arb('1.25506')*m*Pa**(1 - m)/((m - 1)*Pa.log())/((1 - r1P)*(1 - r2P))
dtl = 2*2*arb('1.01624')/Pa/((1 - r2P)*(1 - r1P - r2P))
tau += tl; dtail += dtl; dlall += dtl
print('tau', tau.str(8), ' dtail', dtail.str(8), ' sup|d_w log A| on Re s = k1 <=', dlall.str(8), flush=True)
# ---- [0, 80]: acb.integral ----
def f_g(t, analytic): return gD(acb(kap1) + acb(0, 1)*t)[0]
def f_gD(t, analytic): return gD(acb(kap1) + acb(0, 1)*t)[1]
TT = 80
I0 = acb.integral(f_g, 0, TT, rel_tol=arb('1e-6'), eval_limit=200000)
print('int_0^80 g =', I0, '%.0fs' % (time.time() - t0), flush=True)
I1 = acb.integral(f_gD, 0, TT, rel_tol=arb('1e-6'), eval_limit=200000)
print('int_0^80 gD =', I1, '%.0fs' % (time.time() - t0), flush=True)
trunc_g = (tau.exp() - 1)*absg
trunc_gD = (tau.exp() - 1)*absgD + tau.exp()*dtail*absg          # g D = g_P0 T (D_P0 + dD),  |T - 1| <= e^tau - 1, |dD| <= dtail
# ---- [80, T2] on the vertical line and beyond: |g| <= Ak MI iota(1+k1) e^{eta^2 (k1^2 - t^2)}/t^4 ----
import d4_consts as DC
MI, MZ = arb('14.57'), arb('11.89')                 # B1 (b1_cert.py): |1/zeta| <= 14.57, |zeta| <= 11.89 on 0.6 <= sigma <= 1.3, |t| <= 450
iot = arb('1.38').zeta()/arb('2.76').zeta()        # |1/zeta(sigma + it)| <= zeta(sigma)/zeta(2 sigma)
zl = acb_series([acb('1.38'), 1], prec=2).zeta(); mzz = -(zl[1]/zl[0]).real     # |zeta'/zeta(1.38 + it)| <= -zeta'/zeta(1.38)
Ak = DC.Astar(E1, E2, arb(1), verbose=True)      # sup |A(0; -s, s)| on Re s = k1 (branch and bound over independent phases)
a80 = arb(80)
tailv = Ak*MI*iot*(eta*eta*kap1*kap1).exp()*(-(eta*eta*a80*a80)).exp()/(3*a80**3)
Dv = dlall + eta*eta*(kap1*kap1 + T2*T2).sqrt() + mzz + 2/a80
# ---- horizontal pieces s = sigma +- i T2, -c' <= sigma <= k1: real parts 1 - sigma in [0.62, 1.01], 1 + sigma in [0.99, 1.38] ----
def iv(lo, hi): lo, hi = arb(lo), arb(hi); return (lo + hi)/2 + arb(0, ((hi - lo)/2).upper())
AkH = DC.Astar(iv(1 - kap1, 1 + cp), iv(1 - cp, 1 + kap1), arb(1), crude=True)
hor = 2*(kap1 + cp)*AkH*MI*MI*(eta*eta*kap1*kap1).exp()*(-(eta*eta*T2*T2)).exp()/T2**4
# |D| on the horizontal pieces: |d_w log A| <= sum_p log p r2 (r1 + r3)/((1 - r2)(1 + r3 - r1 - r2)) with the largest moduli
# (e1 = 0.62, e2 = 0.99), Cauchy |zeta'(s)| <= max_{|z - s| = 0.3} |zeta(z)|/0.3 <= 11.89/0.3 (B1 box; |zeta| <= zeta(1.3) <= 11.89 for sigma > 1.3)
dlH = arb(0)
for p in PALL:
    dlH += ap_terms(p, E1, 1 - cp)[1]
r1P, r2P = Pa**(-E1), Pa**(-(1 - cp))
dlH += 2*2*arb('1.01624')*Pa**(-(2 - cp - E1) + 1)/((1 - r2P)*(1 - r1P - r2P))/(1 - E1 - (1 - cp) + 1)   # crude tail, exponent 1.61
Dh = dlH + eta*eta*(kap1*kap1 + T2*T2).sqrt() + MZ/arb('0.3')*MI + 2/T2
print('Ak', Ak.str(6), 'AkH', AkH.str(6), 'Dv', Dv.str(6), 'Dh', Dh.str(6), 'tailv', tailv.str(4), 'hor', hor.str(4), flush=True)
k1 = (I0.real + arb(0, (trunc_g + tailv).upper()))/arb.pi() + arb(0, (hor/(2*arb.pi())).upper())
H1 = (I1.real + arb(0, (trunc_gD + tailv*Dv).upper()))/arb.pi() + arb(0, (hor*Dh/(2*arb.pi())).upper())
k0 = H1 + arb.const_euler()*k1
print('k1 =', k1, '  k0 =', k0, flush=True)
res = {'P0': P0, 'I0_0_80': I0.str(20, radius=True), 'I1_0_80': I1.str(20, radius=True), 'absg': absg.str(10), 'absgD': absgD.str(10),
       'tau': tau.str(10), 'dtail': dtail.str(10), 'trunc_g': trunc_g.str(10), 'trunc_gD': trunc_gD.str(10), 'tailv': tailv.str(5), 'hor': hor.str(5),
       'k1': k1.str(20, radius=True), 'k0': k0.str(20, radius=True), 'k1_upper': repr(float(k1.upper())), 'k0_upper': repr(float(k0.upper()))}
json.dump(res, open('out/kconst_cert.json', 'w'), indent=1)
print('done %.0fs' % (time.time() - t0))

