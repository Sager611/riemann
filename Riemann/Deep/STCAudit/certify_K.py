# Interval-arithmetic certificate of Casana-Giner's kernel inequality (K) for the line l20 (d = 1.2, m = 5; the kernel
# with the smallest sampled margin, 1.01e-4 at b = 1/2, t = 1.006), own code, python-flint arb (80 bits).
# Region A (|t| <= T1): H_s(b,tau) = sum_i eps_s(X) c_i G_i(b,tau) - s f(b,tau) > 0 on [0,1/2] x [0,T1], s = +1 (t = tau > 0)
#   and s = -1 (t = -tau < 0; uses f odd, G_R even, G_I odd), f in the form (4.17) (analytic at tau = 0, the t -> 0 limit).
#   Centred form H(mid) + dH/db(box)[-rb,rb] + dH/dtau(box)[-rt,rt], adaptive bisection.
# Region B (|t| >= T1, the t -> infinity tail, Lemma 4.10): Q(b,u) = t^2 U(b,t) - 4 d^3 u - K5 |u|^3 > 0 for u = 1/t in
#   [-1/T1, 1/T1] via (4.19), natural interval extension, bisection; Lemma 4.10 gives t^2 f <= 4 d^3 u + K5 |u|^3.
import json, sys, time
from flint import arb, acb, ctx
ctx.prec = 80
K = json.load(open('Riemann/Deep/STCAudit/kernels.json'))
name = sys.argv[1] if len(sys.argv) > 1 else 'L20'
T1 = int(sys.argv[2]) if len(sys.argv) > 2 else 100
TLIM = float(sys.argv[3]) if len(sys.argv) > 3 else 1500
k = K[name]; d = arb(k['d']); B = [(X, int(j)) for X, j in k['B']]; C = [arb(c) for c in k['c']]
m = max(j for X, j in B); PI = arb.pi(); I = acb(0, 1)
half = arb(1)/2
def phi(X, w): return w.real if X == 'R' else -w.imag
def H(s, b, t, grad):
    # value (grad=False) or gradient (grad=True) of H_s at arb balls b, t
    val = arb(0); gb = arb(0); gt = arb(0)
    for x, sx in ((d + b, 1), (d - b, -1)):
        z = acb(x, t); zi = 1/z
        p = [acb(1)]
        for j in range(1, m + 2): p.append(p[-1]*zi)
        for (X, j), c in zip(B, C):
            cc = c if (s == 1 or X == 'R') else -c
            if grad:
                gt += cc*phi(X, -I*j*p[j+1]); gb += cc*phi(X, -j*p[j+1])*sx
            else:
                val += cc*phi(X, p[j])
    ys = ((2, d + b, 1), (2, d - b, -1), (-1, 2*d + b, 1), (-1, 2*d - b, -1))
    if grad:
        # f = pi - sum a*atan(t/y): df/dt = -sum a*y/(y^2+t^2); df/db = sum a*t/(y^2+t^2)*dy/db
        dft = -sum((a*y/(y*y + t*t) for a, y, sy in ys), arb(0)); dfb = sum((a*t/(y*y + t*t)*sy for a, y, sy in ys), arb(0))
        return gb - s*dfb, gt - s*dft
    f = PI - sum((a*(t/y).atan() for a, y, sy in ys), arb(0))
    return val - s*f
def certA(s, t_hi, t0):
    stack = []
    bs = [arb(i)/16 for i in range(9)]
    tb = [arb(0), arb(1)/8, arb(1)/4, arb(1)/2]
    v = 1
    while v < t_hi: tb.append(arb(v)); v *= 2
    tb.append(arb(t_hi))
    for i in range(8):
        for jj in range(len(tb)-1): stack.append((bs[i], bs[i+1], tb[jj], tb[jj+1], 0))
    n = 0; minlo = None; maxdepth = 0
    while stack:
        b0, b1, ta, tb_, dep = stack.pop(); n += 1
        bm = (b0 + b1)/2; tm = (ta + tb_)/2; rb = (b1 - b0)/2; rt = (tb_ - ta)/2
        hm = H(s, bm, tm, False)
        if hm.upper() < 0:
            print('COUNTEREXAMPLE s=%d b=%s t=%s H=%s' % (s, bm, tm, hm)); sys.exit(1)
        gb, gt = H(s, arb(bm, rb.upper()), arb(tm, rt.upper()), True)
        enc = hm + gb*arb(0, rb.upper()) + gt*arb(0, rt.upper())
        if enc.lower() > 0:
            lo = enc.lower(); minlo = lo if minlo is None or lo < minlo else minlo; maxdepth = max(maxdepth, dep); continue
        if dep > 60: print('DEPTH LIMIT s=%d box b=[%s,%s] t=[%s,%s]' % (s, b0, b1, ta, tb_)); sys.exit(2)
        if time.time() - t0 > TLIM: print('TIME LIMIT in region A s=%d after %d boxes' % (s, n)); sys.exit(3)
        # bisect the side with the larger contribution to the enclosure width
        if max(abs(float(gb.lower())), abs(float(gb.upper())))*float(rb.upper()) >= max(abs(float(gt.lower())), abs(float(gt.upper())))*float(rt.upper()):
            stack.append((b0, bm, ta, tb_, dep+1)); stack.append((bm, b1, ta, tb_, dep+1))
        else:
            stack.append((b0, b1, ta, tm, dep+1)); stack.append((b0, b1, tm, tb_, dep+1))
    return n, minlo, maxdepth
def Q(b, u, sgn):
    tot = arb(0)
    K5 = (2*(d + half)**5 + 2*d**5 + (2*d + half)**5 + (2*d)**5)/5
    for x in (d + b, d - b):
        w = acb(x*u, 1)
        for (X, j), c in zip(B, C):
            if j == 1: tot += c*x/(1 + x*x*u*u)
            else: tot += c*phi(X, (u**(j-2) if j > 2 else arb(1))*w**(-j))
    return tot - 4*d**3*u - K5*sgn*u**3
def certB(sgn, t0):
    lo_u, hi_u = (arb(0), arb(1)/T1) if sgn > 0 else (-arb(1)/T1, arb(0))
    stack = [(arb(0), half, lo_u, hi_u, 0)]; n = 0; minlo = None
    while stack:
        b0, b1, u0, u1, dep = stack.pop(); n += 1
        e = Q(arb((b0+b1)/2, ((b1-b0)/2).upper()), arb((u0+u1)/2, ((u1-u0)/2).upper()), sgn)
        if e.lower() > 0:
            minlo = e.lower() if minlo is None or e.lower() < minlo else minlo; continue
        if dep > 40 or time.time() - t0 > TLIM: print('FAIL region B sgn=%d box b=[%s,%s] u=[%s,%s] enc=%s' % (sgn, b0, b1, u0, u1, e)); sys.exit(4)
        if dep % 2 == 0: bm = (b0+b1)/2; stack += [(b0, bm, u0, u1, dep+1), (bm, b1, u0, u1, dep+1)]
        else: um = (u0+u1)/2; stack += [(b0, b1, u0, um, dep+1), (b0, b1, um, u1, dep+1)]
    return n, minlo
t0 = time.time()
assert arb(T1) >= 2*d + half
for s in (1, -1):
    n, lo, dep = certA(s, T1, t0)
    print('%s region A s=%+d (t %s 0, |t|<=%d): certified, %d boxes, min lower bound %s, max depth %d, %.0f s' % (name, s, '>' if s > 0 else '<', T1, n, lo, dep, time.time()-t0)); sys.stdout.flush()
for sgn in (1, -1):
    n, lo = certB(sgn, t0)
    print('%s region B u %s 0 (|t|>=%d): certified, %d boxes, min lower bound %s, %.0f s' % (name, '>' if sgn > 0 else '<', T1, n, lo, time.time()-t0)); sys.stdout.flush()
print('%s: (K) CERTIFIED on 0<=b<=1/2, t != 0' % name)
