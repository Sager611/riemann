
# Rigorous check (python-flint arb/acb) of a lattice dual certificate for Lamzouri's Thm Main2 LP.
# Certificate: sigma = sum_{k != 0} h(1+eta_|k|) delta_{kh}, eta_k = -1 for k < k0 (kh < lam), eta_k = 0 for k > K.
# Claim: G(al) = sum_{m!=0} delta_{m/h} + P(al) - al*1[|al|<=theta] >= 0, P(al) = (1-h) + 2h sum_k eta_k cos(2 pi al k h).
import json, sys, time, resource
from fractions import Fraction as Fr
from flint import arb, acb, acb_poly, ctx
ctx.prec = 360  # rectangular acb balls grow like (|Re z|+|Im z|)^k in Horner; 360 bits absorb K <= 400
d = json.load(open(sys.argv[1]))
lam = Fr(str(d['lam'])); k0 = d['k0']; h = lam / k0
theta = Fr(399993, 700000) if abs(d['theta'] - 399993/700000) < 1e-12 else Fr(str(d['theta']))
if len(sys.argv) > 3: theta = Fr(sys.argv[3])   # override, e.g. 4/7
eta_f = d['eta']; K = len(eta_f)
D = 10**12
eta = []
for k, e in enumerate(eta_f, start=1):
    if k < k0: q = Fr(-1)
    else:
        q = Fr(round(e * D), D)
        if q < -1: q = Fr(-1)
    eta.append(q)
assert all(k*h >= lam for k in range(1, K+1) if eta[k-1] != -1), "atom inside hard core"
assert all(1 + q >= 0 for q in eta)
assert all(eta[k-1] == -1 for k in range(1, k0))
ha = arb(h.numerator) / h.denominator
tha = arb(theta.numerator) / theta.denominator
co = [acb(0)] + [acb(arb(q.numerator) / q.denominator) for q in eta]
Q = acb_poly(co); dQ = Q.derivative()
pi = arb.pi()
# |P''| <= M2 = 2h sum |eta_k| (2 pi k h)^2
M2 = arb(0)
for k, q in enumerate(eta, start=1):
    M2 += abs(arb(q.numerator)/q.denominator) * (2*pi*k*ha)**2
M2 = 2*ha*M2
amax = Fr(1, 2) / h
delta = Fr(1, int(sys.argv[2]) if len(sys.argv) > 2 else 5000)
# grid: multiples of delta in [0, amax], plus theta and amax
N = int(amax / delta) + 1
pts = sorted(set([Fr(i) * delta for i in range(N + 1) if Fr(i) * delta <= amax] + [theta, amax]))
def PdP(a):
    aa = arb(a.numerator) / a.denominator
    z = acb(0, 2 * pi * aa * ha).exp()
    S = Q(z); dS = dQ(z)
    P = (1 - ha) + 2 * ha * S.real
    dP = 2 * ha * (-(2 * pi * ha) * (z * dS).imag)
    return P, dP
t0 = time.time()
vals = [PdP(a) for a in pts]
worst = None; worst_at = None
for i in range(len(pts) - 1):
    a0, a1 = pts[i], pts[i + 1]
    dl = arb(( a1 - a0).numerator) / (a1 - a0).denominator
    P, dP = vals[i]
    if a1 <= theta:   # g = P - al
        g0 = P - arb(a0.numerator) / a0.denominator; dg = dP - 1
    else:
        g0 = P; dg = dP
    # concave Taylor minorant on [a0,a1]: g >= min(g0, g0 + dg*dl) - M2 dl^2/2
    lo = g0 + dg * dl
    m = (g0 if g0.lower() < lo.lower() else lo) - M2 * dl * dl / 2
    if worst is None or m.lower() < worst.lower():
        worst = m; worst_at = (float(a0), float(a1))
# endpoint amax itself covered by last cell; left-endpoint lower bounds only (cells are closed)
print(f"lam={float(lam)} h={h} k0={k0} K={K} theta={theta} cells={len(pts)-1} M2={float(M2.mid()):.1f} delta={float(delta)}")
print(f"min lower bound of G on [0,1/(2h)] = {worst.lower()} at cell {worst_at}  ({time.time()-t0:.0f}s, maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB)")
print("PASS" if worst.lower() > 0 else "FAIL")
