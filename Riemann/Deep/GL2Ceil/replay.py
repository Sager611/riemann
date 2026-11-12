#!/usr/bin/env python3
"""gl2ceil INDEPENDENT REPLAY of a certificate cell list (re-implementation; shares only the untrusted float proposer ceil.build).
Integrals: exact rationals (flint fmpq_poly) with the finite formula int_0^1 e^{ax} G = e^a X - Y (Degree2Mol cert.py's method, not the 1F1
moments of ceil.py); one Arb exp per factor. Bounds, Taylor assembly, lambda subcells, crude cells and tails re-written from DERIVATION.md.
usage: uv run --no-project --with numpy --with scipy --with python-flint python3 replay.py logs/ceil_<tag>.json [R0]"""
import sys, os, json, time, math
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
cert = json.load(open(sys.argv[1]))
xn, xd = cert['xi'].split('/')
os.environ['HLBITS'] = str(round(-math.log2(cert['HL'])))
sys.argv = ['x', xn, xd, cert['kappa_up']]; sys.path.insert(0, HERE)
import ceil as PROPOSER                        # used ONLY for build() (float proposals of sigma, mu)
from flint import fmpq, fmpq_poly, arb, ctx
ctx.prec = 320
def Q(fr): fr = Fraction(fr); return fmpq(fr.numerator, fr.denominator)
XI = fmpq(int(xn), int(xd)); KUP = Q(Fraction(cert['kappa_up'])); HL = Q(Fraction(cert['HL']))
KP, DEG = cert['KP'], cert['DEG']; RMAX = Q(Fraction(cert['RMAX'])); LAM = Q(Fraction(cert['LAM']))
ONE_K = arb(1) - arb(KUP)
def tgt(R): return (arb(R)*ONE_K).exp()
def CB(i):
    if i == 0: return arb(1), arb(1)/3
    l = arb(fmpq(i)*HL); sh = l.sinh(); t = (2*l).sinh()/(4*l)
    return l*l*(arb(1)/2 + t)/(sh*sh), (t - arb(1)/2)/(sh*sh)
CBc = {}
def CBm(i):
    if i not in CBc: CBc[i] = CB(i)
    return CBc[i]
def iexp(G, a):                                 # exact (Xc, Yc) with int_0^1 e^{ax} G dx = e^a Xc - Yc, a != 0
    X = fmpq(0); Y = fmpq(0); D = G; k = 0; ap = a
    while not D.is_zero():
        s = 1 if k % 2 == 0 else -1
        X += s*D(fmpq(1))/ap; Y += s*D(fmpq(0))/ap; D = D.derivative(); ap *= a; k += 1
    return X, Y
def ipoly(G): I = G.integral(); return I(fmpq(1)) - I(fmpq(0))
x_ = fmpq_poly([0, 1]); omx = fmpq_poly([1, -1])
def compose_1mx(p):                             # p(1-x) by Horner
    c = p.coeffs(); r = fmpq_poly([0])
    for cj in reversed(c): r = r*omx + fmpq_poly([cj])
    return r
def quantities(sc, mc, Rc, Rlo):
    s = []
    for k in range(KP):
        c = [Q(float(t)) for t in sc[k]]
        if k: c[0] = s[-1](fmpq(1))             # exact continuity
        s.append(fmpq_poly(c))
    m = [fmpq_poly([Q(float(t)) for t in mc[j]]) for j in range(KP//2)]
    mt = [m[k - KP//2] if k >= KP//2 else compose_1mx(m[KP//2 - 1 - k]) for k in range(KP)]
    a1 = -Rc/KP; a2 = -2*Rc/KP; a3 = -2*Rlo/KP
    E1 = [(-arb(Rc)*j/KP).exp() for j in range(KP+1)]
    E2 = [(-2*arb(Rc)*j/KP).exp() for j in range(KP+1)]
    E3 = [(-2*arb(Rlo)*j/KP).exp() for j in range(KP+1)]
    K2 = S1 = M = fmpq(0); X = Xp = Xpp = Y = Yp = Ypp = Ylo = arb(0)
    for k in range(KP):
        sk = s[k]; dsk = sk.derivative(); mk = mt[k]; yk = (fmpq(k) + x_)*fmpq(1, KP); yk2 = yk*yk
        K2 += ipoly(sk*sk)/KP; S1 += ipoly(dsk*dsk)*KP; M += ipoly(mk)/KP
        g = dsk*mk; h = mk*mk
        for G, sgn, acc in ((g, 1, 'X'), (yk*g, -1, 'Xp'), (yk2*g, 1, 'Xpp')):
            Xc, Yc = iexp(G, a1); v = sgn*(E1[k+1]*arb(Xc) - E1[k]*arb(Yc))
            if acc == 'X': X += v
            elif acc == 'Xp': Xp += v
            else: Xpp += v
        for G, fac, acc in ((h, 1, 'Y'), (yk*h, -2, 'Yp'), (yk2*h, 4, 'Ypp')):
            Xc, Yc = iexp(G, a2); v = fac*(E2[k+1]*arb(Xc) - E2[k]*arb(Yc))/KP
            if acc == 'Y': Y += v
            elif acc == 'Yp': Yp += v
            else: Ypp += v
        Xc, Yc = iexp(h, a3); Ylo += (E3[k+1]*arb(Xc) - E3[k]*arb(Yc))/KP
    s0 = arb(s[0](fmpq(0))); s1 = arb(s[-1](fmpq(1))); Ma = arb(M)
    eR = E1[KP]; eL = (-arb(Rlo)).exp()
    t = 2*s1 + Ma*eR
    K0 = arb(1)/2 + Ma - 2*s0 - t*t/2; dK0 = t*Ma*eR
    F0 = Ma*Ma*eL*eL + (abs(2*s1) + abs(Ma)*eL)*abs(Ma)*eL
    K1 = arb(S1) + 2*X + Y; dK1 = 2*Xp + Yp; ddK1 = 2*Xpp + Ypp
    G1 = 2*arb(S1).sqrt()*Ylo.sqrt() + 8*Ylo
    return K0, dK0, F0, K1, dK1, ddK1, G1, arb(K2)
def check(qq, Rc, d, ia, ib):
    K0, dK0, F0, K1, dK1, ddK1, G1, K2 = qq; da = arb(d); tA, tB = tgt(Rc - d), tgt(Rc + d)
    stack = [(ia, ib)]
    while stack:
        i, j = stack.pop()
        Ci, _ = CBm(i); _, Bj = CBm(j); x = arb(XI)/Ci; z = 1/(arb(XI)*Bj)
        q2 = (-x*ddK1/2).lower(); q2 = q2 if q2 < 0 else arb(0)
        base = K0 - x*K1 - z*K2 + q2*da*da - F0*da*da/2 - x*G1*da*da*da/6
        lo_m = base - dK0*da + x*dK1*da; lo_p = base + dK0*da - x*dK1*da
        if lo_m > tA and lo_p > tB: continue
        if j - i <= 1: return False
        h = (i + j)//2; stack += [(i, h), (h, j)]
    return True
def crude(R, C):
    R = arb(R); gd = 2*((R/2).tanh()).atan(); K = C/arb(XI)*R.exp()*gd/(4*R); e = (2*R).exp()
    return arb(1)/2 + e*K/(e + 2*K)
t0 = time.time()
cells = cert['cells']; R0 = min(Q(Fraction(c[0])) - Q(Fraction(c[1])) for c in cells)
# tails, recomputed
xi = arb(XI); m0 = arb(2) - 4*xi/(3*xi + 1); s = xi/arb(3).sqrt(); Rt = m0.log()/(1 + s - arb(KUP))
iL = int(Fraction(cert['LAM'])/Fraction(cert['HL'])); CL, BL = CBm(iL); m0L = arb(2) - 1/(arb(3)/4 + CL/(4*xi)); RtL = m0L.log()/(1 + xi*(BL/CL).sqrt() - arb(KUP))
gdM = 2*((arb(RMAX)/2).tanh()).atan(); hM = ((4*arb(RMAX)*xi/gdM).log() + 2*(-arb(RMAX)).exp())/arb(RMAX)
ok_tails = bool(Rt > arb(R0)) and bool(hM < arb(KUP)) and bool(arb(RMAX) > arb(1).exp()*gdM/(4*xi))
Rs = fmpq(int(max(float(Rt.lower()), float(RtL.lower()))*64), 64)
st = [(Rs, RMAX + fmpq(1, 64))]; ok_strip = True; nst = 0
while st:
    a, b = st.pop()
    if crude(a, CL) > tgt(b): nst += 1; continue
    if b - a < fmpq(1, 2**20): ok_strip = False; break
    mid = (a + b)/2; st += [(a, mid), (mid, b)]
ok_strip = ok_strip and bool(arb(Rs) <= RtL) or (ok_strip and bool(arb(Rs) <= Rt))
print(f"tails: R_tail={Rt.str(10)} (core starts {float(R0)}), R_tail(lam>={float(LAM)})={RtL.str(10)}, h(RMAX)={hM.str(10)}: {'OK' if ok_tails else 'FAIL'}; strip from {float(Rs)}: {'OK' if ok_strip else 'FAIL'} ({nst} intervals)", flush=True)
nd = nc = bad = 0
for n, (Rcs, ds, ia, ib, how) in enumerate(cells):
    Rc = Q(Fraction(Rcs)); d = Q(Fraction(ds))
    if how == 'crude':
        C, _ = CBm(ia); ok = bool(crude(Rc - d, C) > tgt(Rc + d)); nc += 1
    else:
        Rf = float(Fraction(Rcs)); lc = (ia + ib)/2*cert['HL']          # same float inputs as the original run
        sc, mc, beta, J = PROPOSER.build(Rf, lc)
        ok = check(quantities(sc, mc, Rc, Rc - d), Rc, d, ia, ib); nd += 1
    if not ok: bad += 1; print("REPLAY FAIL cell", Rcs, ds, ia, ib, how, flush=True)
    if n % 200 == 0: print(f"  {n}/{len(cells)} replayed, {bad} failures, {time.time()-t0:.0f}s", flush=True)
verdict = ok_tails and ok_strip and bad == 0
print(f"REPLAY {sys.argv[3] if len(sys.argv) > 3 else cert['kappa_up']} xi={cert['xi']}: {nd} dual + {nc} crude cells, {bad} failures; tails {ok_tails}, strip {ok_strip} -> {'REPLAY OK' if verdict else 'REPLAY FAIL'} ({time.time()-t0:.0f}s)")
