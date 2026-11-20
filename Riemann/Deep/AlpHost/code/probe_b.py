# probe (b): certified zero-free boxes Re s > 1 - 1/400, |Im s| <= T for small real primitive characters
# argument principle on [s0, 3] x [-T, T] with arb ball enclosures (python-flint); counts zeros (expect 0)
import sys, time, math
from flint import acb, arb, dirichlet_char, ctx
ctx.prec = 64
s0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1 - 1/400.0   # left edge; argv[3] = 0.4 for the positive control
S1 = 3.0
def kron(D, n):
    # Kronecker symbol (D/n) for n >= 1
    if math.gcd(D, n) != 1: return 0
    r = 1
    while n % 2 == 0:
        n //= 2
        if D % 8 in (3, 5): r = -r
    # Jacobi (D/n), n odd
    a, m = D % n, n
    while a:
        while a % 2 == 0:
            a //= 2
            if m % 8 in (3, 5): r = -r
        a, m = m, a
        if a % 4 == 3 and m % 4 == 3: r = -r
        a %= m
    return r if m == 1 else 0
def find_char(D):
    q = abs(D)
    for l in range(1, q + 1):
        if math.gcd(l, q) != 1: continue
        c = dirichlet_char(q, l)
        if not (c.is_real() and c.is_primitive()): continue
        if all(abs(complex(c(n)) - kron(D, n)) < 1e-9 for n in range(1, 4 * q + 1)):
            return c
    return None
def ball(a, b):
    m = (a + b) / 2; r = abs(b - a) / 2
    if a.real == b.real:  # vertical segment
        return acb(arb(m.real), arb(m.imag, r))
    return acb(arb(m.real, r), arb(m.imag))
def L(chi, z):
    return acb.dirichlet_l(z, chi)
def seg_ok(E):
    rm, rr = float(E.real.mid()), float(E.real.rad()); im_, ir = float(E.imag.mid()), float(E.imag.rad())
    mod = math.hypot(rm, im_); rad = math.hypot(rr, ir)
    return mod > 0 and rad < mod * math.sin(math.pi / 4) * 0.999
def edge(chi, a, b, h, stats):
    # returns total arg increment along a->b (rigorous up to float rounding of principal args; each |inc| < pi/2 certified)
    n = max(1, int(math.ceil(abs(b - a) / h)))
    pts = [a + (b - a) * k / n for k in range(n + 1)]
    tot = 0.0
    stack = [(pts[k], pts[k + 1], 0) for k in range(n)][::-1]
    val = {}
    def pv(z):
        if z not in val: val[z] = complex(L(chi, acb(z.real, z.imag)).mid())
        return val[z]
    while stack:
        u, v, dep = stack.pop()
        E = L(chi, ball(u, v)); stats['ev'] += 1
        if seg_ok(E):
            inc = math.atan2((pv(v) / pv(u)).imag, (pv(v) / pv(u)).real)
            tot += inc
        else:
            if dep > 30: raise RuntimeError('cannot certify segment %s %s' % (u, v))
            w = (u + v) / 2; stats['bis'] += 1
            stack.append((w, v, dep + 1)); stack.append((u, w, dep + 1))
    return tot
Ds = [int(x) for x in sys.argv[1].split(',')]
Tcap = float(sys.argv[2]) if len(sys.argv) > 2 else 200.0
for D in Ds:
    t0 = time.time()
    chi = find_char(D)
    if chi is None: print('D', D, 'no char found'); continue
    T = min(8.0 * abs(D), Tcap)
    st = {'ev': 0, 'bis': 0}
    c1, c2, c3, c4 = complex(s0, -T), complex(S1, -T), complex(S1, T), complex(s0, T)
    tot = edge(chi, c1, c2, 0.05, st) + edge(chi, c2, c3, 1.0, st) + edge(chi, c3, c4, 0.05, st) + edge(chi, c4, c1, 0.05, st)
    w = tot / (2 * math.pi)
    print('D=%d conductor=%d parity=%d T=%.0f s0=%.4f: zeros in (s0, 3] x [-T, T] = %d (winding %.6f) evals=%d bisections=%d %.1fs'
          % (D, chi.conductor(), chi.parity(), T, s0, round(w), w, st['ev'], st['bis'], time.time() - t0), flush=True)
