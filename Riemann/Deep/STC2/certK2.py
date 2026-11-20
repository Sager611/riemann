# certK2.py -- own certificate (stc2) of Casana-Giner's kernel inequality (K): f(b,t) <= U(b,t) for 0 <= b <= 1/2, t != 0
# (Zenodo 23866074, (K), (4.17), Lemma 4.10).  python-flint arb/acb, 80 bits.  Usage: certK2.py NAME edge|box [T1] [TLIM]
# H_s(b,tau) = U_s(b,tau) - s f(b,tau), tau >= 0: s = +1 is t = tau > 0; s = -1 is t = -tau < 0 (f odd, G_R even, G_I odd).
# Structure used: H_s = Re[ sum_k Psi_k(w_k) ] - s pi with w_k = y_k + i tau, y_k in {d+b, d-b, 2d+b, 2d-b},
#   Psi(w) = sum_j A_j w^-j - 2 i s log w (k = 1,2), Psi(w) = i s log w (k = 3,4), A_j = c_Rj + i c_Ij (eps_X c for s = -1).
#   Hence H_bb = -H_tt = Re sum Psi_k'', H_bt = -Im sum sigma_k Psi_k'', H_b = Re sum sigma_k Psi_k', H_t = -Im sum Psi_k'.
# edge: H_s is harmonic in (b,tau) near [-1/2,1/2] x [0,oo), even in b, and -> 0 uniformly as tau -> oo, so by the minimum
#   principle H_s >= 0 there once H_s >= 0 on b = 1/2 (tau >= 0) and on tau = 0. Certified: (a) U(b,0) - pi > 0 on [0,1/2]
#   (this is H_+(b,0); H_-(b,0) = U(b,0) + pi); (b) H_s(1/2,tau) > 0 on [0,T1]; (c) Q_s(1/2,u) = t^2 U - 4d^3 u - K5|u|^3 > 0
#   on 0 < s u <= 1/T1 (Lemma 4.10: t^2 f <= 4d^3 u + K5|u|^3 for |t| >= T1 >= 2d+1/2).  1-D second-order forms
#   v(m) - |v'(m)| r - |v''(box)| r^2/2; accepted only when also within EPS of the centre value (so the printed minima are
#   certified lower bounds within ~1% of the true minima).
# box: no harmonicity used. 2-D bisection of [0,1/2] x [0,T1] with H(m) - |H_b(m)| r_b - |H_t(m)| r_t
#   - |H_bb(box)| (r_b^2 + r_t^2)/2 - |H_bt(box)| r_b r_t > 0; region B: Q on [0,1/2] x [-1/T1,1/T1] with the second-order
#   form in b (Q(b_m,U) - |Q_b(b_m,U)| r_b - |Q_bb(box)| r_b^2/2, natural extension in u).
import json, sys, time
from fractions import Fraction as Fr
from flint import arb, acb, ctx
ctx.prec = 80
KJ = json.load(open('Riemann/Deep/STC2/kernels.json'))
name = sys.argv[1]; MODE = sys.argv[2]
T1 = int(sys.argv[3]) if len(sys.argv) > 3 else 100
TLIM = float(sys.argv[4]) if len(sys.argv) > 4 else 1500
EPS = 0.01
k = KJ[name]
dq = Fr(k['d']).limit_denominator(10**6); d = arb(dq.numerator)/dq.denominator
B = [(X, int(j)) for X, j in k['B']]; C = [arb(c) for c in k['c']]
assert ('I', 1) not in B and B[0] == ('R', 1)
m = max(j for X, j in B); PI = arb.pi(); half = arb(1)/2
cR1 = C[B.index(('R', 1))]; cR2 = C[B.index(('R', 2))]
K5 = (2*(d + half)**5 + 2*d**5 + (2*d + half)**5 + (2*d)**5)/5
assert arb(T1) >= 2*d + half
def Acoef(s):
    A = [acb(0) for _ in range(m+1)]
    for (X, j), c in zip(B, C):
        cs = c if (s == 1 or X == 'R') else -c
        A[j] += acb(cs) if X == 'R' else acb(0, cs)
    return A
AC = {1: Acoef(1), -1: Acoef(-1)}; A0 = Acoef(1)
FT = [(2, d, 1), (2, d, -1), (-1, 2*d, 1), (-1, 2*d, -1)]
def vg(s, b, t):
    A = AC[s]; U = arb(0); S1 = acb(0); S0 = acb(0)
    for sg in (1, -1):
        w = acb(d + sg*b, t); wi = 1/w; p = wi; Fv = acb(0); Fd = acb(0)
        for j in range(1, m+1):
            Fv += A[j]*p; p = p*wi; Fd -= j*A[j]*p
        U += Fv.real; S1 += sg*Fd; S0 += Fd
    f = PI
    for al, base, sg in FT:
        y = base + sg*b; f -= al*(t/y).atan()
        dd = acb(0, -s*al)/acb(y, t); S1 += sg*dd; S0 += dd
    return U - s*f, S1.real, -S0.imag
def hess(s, b, t):
    A = AC[s]; T0 = acb(0); T1_ = acb(0)
    for sg in (1, -1):
        w = acb(d + sg*b, t); wi = 1/w; p = wi*wi; F2 = acb(0)
        for j in range(1, m+1):
            p = p*wi; F2 += (j*(j+1))*A[j]*p
        T0 += F2; T1_ += sg*F2
    for al, base, sg in FT:
        gi = 1/acb(base + sg*b, t); dd = acb(0, s*al)*gi*gi; T0 += dd; T1_ += sg*dd
    return T0.real, -T1_.imag
def ball(lo, hi): return arb((lo + hi)/2, ((hi - lo)/2).upper())
def mag(x): return max(abs(float(x.lower())), abs(float(x.upper())))
def Qv(b, u, sgn):   # t^2 U - 4 d^3 u - K5 |u|^3, natural extension (old certify_K.py form (4.19))
    tot = arb(0)
    for x in (d + b, d - b):
        w = acb(x*u, 1)
        for (X, j), c in zip(B, C):
            if j == 1: tot += c*x/(1 + x*x*u*u)
            else:
                v = (u**(j-2) if j > 2 else arb(1))*w**(-j)
                tot += c*(v.real if X == 'R' else -v.imag)
    return tot - 4*d**3*u - K5*sgn*u**3
def Qder(b, u):      # (d/db, d^2/db^2) of t^2 U: sum_sigma sigma G_x, sum_sigma G_xx, G = Re sum_j A_j u^(j-2) (xu+i)^-j
    g1 = arb(0); g2 = arb(0)
    for sg in (1, -1):
        x = d + sg*b; wi = 1/acb(x*u, 1); p = wi*wi; Gx = acb(0); Gxx = acb(0); up = arb(1)
        for j in range(1, m+1):
            Gx -= j*A0[j]*up*p; p = p*wi; Gxx += (j*(j+1))*A0[j]*up*u*p; up = up*u
        g1 += sg*Gx.real; g2 += Gxx.real
    return g1, g2
t0 = time.time()
def tl():
    if time.time() - t0 > TLIM: print('TIME LIMIT'); sys.exit(3)
def cert1d(fun, segs, label):
    # fun(lo, hi) -> (lower, centre value, scale); bisect until lower > 0 and lower >= min((1-EPS) centre, 1.5 curmin)
    stack = list(segs); n = 0; cmin = None; curmin = None; worst = None
    while stack:
        lo, hi, dep = stack.pop(); n += 1; tl()
        low, cv, sc, scc = fun(lo, hi)
        if cv.upper() < 0: print('COUNTEREXAMPLE %s at [%s,%s] value %s' % (label, lo, hi, cv)); sys.exit(1)
        cvs = float(cv.lower())*scc; curmin = cvs if curmin is None or cvs < curmin else curmin
        ls = float(low.lower())*sc if low.lower() > 0 else None
        if ls is not None and ls >= min((1 - EPS)*cvs, 1.5*curmin):
            cmin = ls if cmin is None or ls < cmin else cmin; continue
        if dep > 60: print('DEPTH LIMIT %s at [%s,%s]' % (label, lo, hi)); sys.exit(2)
        mid = (lo + hi)/2; stack += [(lo, mid, dep+1), (mid, hi, dep+1)]
    print('%s %s: certified, %d intervals, certified min (scaled) >= %.6e, smallest centre value (scaled) %.6e, %.1f s' % (name, label, n, cmin, curmin, time.time()-t0)); sys.stdout.flush()
    return cmin
def bottom(lo, hi):   # (a): H_+(b,0) = U(b,0) - pi on b in [lo,hi]
    bm = (lo + hi)/2; r = (hi - lo)/2
    H, Hb, Ht = vg(1, bm, arb(0)); hbb, hbt = hess(1, ball(lo, hi), arb(0))
    return H - abs(Hb)*r - abs(hbb).upper()*r*r/2, H, 1.0, 1.0
def side(s):          # (b): H_s(1/2,tau) on tau in [lo,hi], scaled by max(1,tau^2)
    def fn(lo, hi):
        tm = (lo + hi)/2; r = (hi - lo)/2
        H, Hb, Ht = vg(s, half, tm); hbb, hbt = hess(s, half, ball(lo, hi))
        return H - abs(Ht)*r - abs(hbb).upper()*r*r/2, H, max(1.0, float(lo.lower())**2), max(1.0, float(tm.lower())**2)
    return fn
def tailQ(sgn):       # (c): Q at b = 1/2 on u in [lo,hi] (natural extension in u)
    def fn(lo, hi):
        q = Qv(half, ball(lo, hi), sgn); qm = Qv(half, (lo + hi)/2, sgn)
        return q, qm, 1.0, 1.0
    return fn
def tsegs():
    tb = [arb(0), arb(1)/8, arb(1)/4, arb(1)/2]; v = 1
    while v < T1: tb.append(arb(v)); v *= 2
    tb[-1] = arb(T1) if tb[-1] > T1 else tb[-1]
    if tb[-1] < T1: tb.append(arb(T1))
    return [(tb[i], tb[i+1], 0) for i in range(len(tb)-1)]
if MODE == 'edge':
    res = {}
    res['a'] = cert1d(bottom, [(arb(i)/16, arb(i+1)/16, 0) for i in range(8)], '(a) tau=0: U(b,0)-pi, b in [0,1/2]')
    for s in (1, -1):
        res['b%+d' % s] = cert1d(side(s), tsegs(), '(b) b=1/2, s=%+d: max(1,tau^2) H_s, tau in [0,%d]' % (s, T1))
        uu = [(arb(0), arb(1)/T1, 0)] if s == 1 else [(-arb(1)/T1, arb(0), 0)]
        res['c%+d' % s] = cert1d(tailQ(s), uu, '(c) b=1/2, s=%+d: Q (Lemma 4.10 tail), |t|>=%d' % (s, T1))
    print('%s: (K) CERTIFIED on 0<=b<=1/2, t != 0 (harmonic reduction to the boundary) %s' % (name, json.dumps(res)))
else:
    # region A: 2-D second-order forms
    for s in (1, -1):
        bs = [arb(i)/16 for i in range(9)]; stack = [(bs[i], bs[i+1], a, b_, 0) for i in range(8) for (a, b_, _) in tsegs()]
        n = 0; minlo = None
        while stack:
            b0, b1, ta, tb_, dep = stack.pop(); n += 1; tl()
            bm = (b0 + b1)/2; tm = (ta + tb_)/2; rb = (b1 - b0)/2; rt = (tb_ - ta)/2
            H, Hb, Ht = vg(s, bm, tm)
            if H.upper() < 0: print('COUNTEREXAMPLE s=%d b=%s t=%s H=%s' % (s, bm, tm, H)); sys.exit(1)
            hbb, hbt = hess(s, ball(b0, b1), ball(ta, tb_)); qa = abs(hbb).upper(); qb = abs(hbt).upper()
            low = H - abs(Hb)*rb - abs(Ht)*rt - qa*(rb*rb + rt*rt)/2 - qb*rb*rt
            if low.lower() > 0:
                lo = low.lower(); minlo = lo if minlo is None or lo < minlo else minlo; continue
            if dep > 60: print('DEPTH LIMIT'); sys.exit(2)
            cb = mag(Hb)*float(rb) + float(qa)*float(rb)**2/2 + float(qb)*float(rb*rt)/2
            ct = mag(Ht)*float(rt) + float(qa)*float(rt)**2/2 + float(qb)*float(rb*rt)/2
            if cb >= ct: stack += [(b0, bm, ta, tb_, dep+1), (bm, b1, ta, tb_, dep+1)]
            else: stack += [(b0, b1, ta, tm, dep+1), (b0, b1, tm, tb_, dep+1)]
        print('%s box region A s=%+d (|t|<=%d): certified, %d boxes, min lower bound %s, %.0f s' % (name, s, T1, n, minlo, time.time()-t0)); sys.stdout.flush()
    for sgn in (1, -1):
        U0 = (arb(0), arb(1)/T1) if sgn > 0 else (-arb(1)/T1, arb(0))
        stack = [(arb(0), half, U0[0], U0[1], 0)]; n = 0; minlo = None
        while stack:
            b0, b1, u0, u1, dep = stack.pop(); n += 1; tl()
            bm = (b0 + b1)/2; rb = (b1 - b0)/2; ub = ball(u0, u1)
            q = Qv(bm, ub, sgn); g1, g2c = Qder(bm, ub); g1x, g2 = Qder(ball(b0, b1), ub)
            low = q - abs(g1)*rb - abs(g2).upper()*rb*rb/2
            if low.lower() > 0:
                lo = low.lower(); minlo = lo if minlo is None or lo < minlo else minlo; continue
            if dep > 50: print('FAIL region B box b=[%s,%s] u=[%s,%s] low=%s' % (b0, b1, u0, u1, low)); sys.exit(4)
            wb = mag(g1)*float(rb) + mag(g2)*float(rb)**2/2; wu = float(q.rad())*2
            if wb >= wu: stack += [(b0, bm, u0, u1, dep+1), (bm, b1, u0, u1, dep+1)]
            else: um = (u0 + u1)/2; stack += [(b0, b1, u0, um, dep+1), (b0, b1, um, u1, dep+1)]
        print('%s box region B u %s 0 (|t|>=%d): certified, %d boxes, min lower bound %s, %.0f s' % (name, '>' if sgn > 0 else '<', T1, n, minlo, time.time()-t0)); sys.stdout.flush()
    print('%s: (K) CERTIFIED on 0<=b<=1/2, t != 0 (2-D boxes, no harmonicity)' % name)

