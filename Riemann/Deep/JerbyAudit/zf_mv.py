# JerbyAudit: independent rigorous re-certification of the zero-free Gram blocks by the MEAN-VALUE form (different algorithm from
# certify_zerofree.py's naive ball evaluation): on [a, b] with centre c, radius r: Z([a,b]) in Z(c) + [-r, r] * Z'([a,b]).
# Gram endpoints by own interval Newton in Arb (start: own asymptotic + float Newton, no mpmath.grampoint); cover starts 1e-12 below
# the left Gram ball and ends 1e-12 above the right one. Also: theta'(7) > 0 certified.  usage: zf_mv.py n [n ...]
import sys, time, math
from flint import arb, acb, ctx, acb_series
ctx.prec = 160
PI = arb.pi()
def theta(t): return acb(arb(0.25), t/2).lgamma().imag - t/2*PI.log()
def dtheta(t): return acb(arb(0.25), t/2).digamma().real/2 - PI.log()/2
def gram_ball(k):
    # asymptotic start t log(t/(2 pi e)) = 2 pi (k + 1/8), float Newton, then Arb Newton and an interval-Newton check
    t = 2*math.pi*(k + 0.125)/math.log(k + 1.0)
    for _ in range(60): t = t - (t*math.log(t/(2*math.pi*math.e)) - 2*math.pi*(k + 0.125))/math.log(t/(2*math.pi))
    m = arb(t)
    for _ in range(10): m = arb((m - (theta(m) - k*PI)/dtheta(m)).mid())
    X = arb(m.mid(), arb('1e-30')); N = arb(m.mid()) - (theta(arb(m.mid())) - k*PI)/dtheta(X)
    assert X.contains(N), 'interval Newton failed'
    return N
def ZdZ(t):
    s = acb_series([acb(arb(0.5), t), acb(0, 1)], prec=2); z = s.zeta().coeffs()
    th = theta(t); e = acb(th.cos(), th.sin())
    return (e*z[0]).real, (e*(acb(0, dtheta(t))*z[0] + z[1])).real
print(f"theta'(7) in {dtheta(arb(7)).str(10)} (> 0: {bool(dtheta(arb(7)) > 0)})", flush=True)
for n in [int(x) for x in sys.argv[1:]]:
    t0 = time.time()
    gl, gn, gr = gram_ball(n - 1), gram_ball(n), gram_ball(n + 1)
    a0 = (gl.lower() - arb('1e-12')).mid(); b0 = (gr.upper() + arb('1e-12')).mid()
    stack = [(arb(b0), arb(a0))]; pieces = 0; sign = None; covered_to = arb(a0); evals = 0
    # process left to right so that coverage is contiguous: stack holds (b, a), pop the leftmost
    todo = [(arb(a0), arb(b0))]
    while todo:
        a, b = todo.pop()
        c = arb(((a + b)/2).mid()); r = arb(max(abs(float((c - a).upper())), abs(float((b - c).upper())))*1.0000001)
        Zc, _ = ZdZ(c); _, dZI = ZdZ(arb(c, r)); evals += 2
        enc = Zc + arb(0, r)*dZI
        if enc > 0 or enc < 0:
            s = 1 if enc > 0 else -1
            assert sign is None or s == sign, 'sign change'
            assert abs(float((a - covered_to).mid())) == 0.0, 'gap in cover'
            sign = s; pieces += 1; covered_to = b
        else:
            if float((b - a).mid()) < 1e-9: raise SystemExit(f'cannot exclude a zero near {c.str(20)}')
            todo.append((c, b)); todo.append((a, c))   # left half popped first
    assert covered_to == arb(b0) or abs(float((covered_to - arb(b0)).mid())) == 0.0
    ok = bool(arb(a0) < gl) and bool(arb(b0) > gr)
    print(f'n={n}: g_(n-1) in {gl.str(25)}, g_(n+1) in {gr.str(25)}; mean-value cover of [{arb(a0).str(22)}, {arb(b0).str(22)}] '
          f'(contains the block: {ok}), {pieces} contiguous pieces, {evals} Arb evaluations, constant sign {sign:+d}: zero-free (certified)  [{time.time() - t0:.1f} s]', flush=True)
