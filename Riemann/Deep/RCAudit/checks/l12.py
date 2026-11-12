
# Lemma 1 / Lemma 2 audit (Combos NOTE §1): exact sharp constant of the weighted large sieve.
# For D<0, classes C <-> reduced forms Q. Primitive O-ideals prime to f in class [a_Q]^{-1} <-> primitive (x,y), Q(x,y)=N,
# gcd(N,f)=1, w-to-1 (w=|O^*|). All ideals m*b, (m,f)=1, same class. Sharp constant of
#   sum_chi |sum c_a chi(a)|^2 <= K * sum m(a)^{1+2e}|c_a|^2   is   K = h * max_C sum_{a in C, Na<=N} m(a)^{-1-2e}.
import numpy as np, math, sys, time, json, resource
OUT = sys.argv[1]; DMAX = int(sys.argv[2]); MODE = sys.argv[3]   # MODE: 'all' (N=|D|^.5,|D|,4|D|) or 'big' (N up to |D|^2, sample)
EPS = [0.5, 0.1, 0.01]
def fund(D):
    n = -D; best = 1
    for g in range(1, math.isqrt(n) + 1):
        if n % (g*g) == 0 and (D // (g*g)) % 4 in (0, 1): best = g
    return D // (best*best), best
def red_forms(D):
    n = -D; F = []; a = 1
    while 3*a*a <= n:
        for b in range(-a+1, a+1):
            if (b*b - D) % (4*a): continue
            c = (b*b - D) // (4*a)
            if c < a or (b < 0 and c == a): continue
            if math.gcd(math.gcd(a, abs(b)), c) != 1: continue
            F.append((a, b, c))
        a += 1
    return F
def prim_norms(a, b, c, N, n):
    V = math.isqrt(4*a*N // n) + 1
    ys = np.arange(-V, V+1, dtype=np.int64)
    disc = N - n*ys.astype(np.float64)**2/(4*a)
    ok = disc >= 0; ys = ys[ok]; disc = disc[ok]
    half = np.sqrt(disc/a); cen = -b*ys/(2.0*a)
    xlo = np.floor(cen - half).astype(np.int64) - 1; xhi = np.ceil(cen + half).astype(np.int64) + 1
    cnt = xhi - xlo + 1
    Y = np.repeat(ys, cnt); st = np.repeat(xlo, cnt)
    off = np.arange(cnt.sum(), dtype=np.int64) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    X = st + off
    Qv = a*X*X + b*X*Y + c*Y*Y
    k = (Qv <= N) & (Qv > 0); X = X[k]; Y = Y[k]; Qv = Qv[k]
    return np.sort(Qv[np.gcd(X, Y) == 1])
t0 = time.time(); rows = []
Ds = [D for D in range(-3, -DMAX-1, -1) if D % 4 in (0, 1)]
if MODE == 'big':
    Ds = [D for D in Ds if -D <= 400] + [D for D in Ds if 400 < -D and (-D) % 97 in (0, 1, 4)]
for D in Ds:
    n = -D; dK, f = fund(D); F = red_forms(D); h = len(F); w = 6 if D == -3 else (4 if D == -4 else 2)
    if MODE == 'all': Ns = sorted({max(1, int(n**0.5)), n, 4*n})
    else: Ns = sorted({max(1, int(n**0.5)), n, int(n**1.5), n*n if n <= 400 else int(n**1.5)})
    Nmax = Ns[-1]
    data = []
    for (a, b, c) in F:
        q = prim_norms(a, b, c, Nmax, n)
        qf = q[np.gcd(q, f) == 1] if f > 1 else q
        data.append((a, q, qf))
    for N in Ns:
        bL1 = 2 + 8*N/n**0.5 + 4*N**0.5/n**0.25
        P = max(int(np.searchsorted(q, N, 'right')) for (_, q, _) in data)
        ms = np.array([m for m in range(1, math.isqrt(N) + 1) if math.gcd(m, f) == 1], dtype=np.int64)
        lims = N // (ms*ms)
        best = {e: 0.0 for e in EPS}; unw = 0.0; piC = 0.0
        for (_, q, qf) in data:
            cnts = np.searchsorted(qf, lims, 'right') / w
            piC = max(piC, cnts[0]); unw = max(unw, cnts.sum())
            for e in EPS: best[e] = max(best[e], float((cnts * ms.astype(float)**(-1-2*e)).sum()))
        base = N + n**0.5; L2 = math.log(3*n)**2
        rows.append(dict(D=D, dK=dK, f=f, h=h, N=N, L1=P/bL1, hpi=h*piC/base, unw=h*unw/base,
                         **{f'K{e}': h*best[e]/base for e in EPS}, **{f'R{e}': h*best[e]/((30/e)*L2*base) for e in EPS}))
    if time.time() - t0 > 1500: print('TIME CAP at D', D); break
json.dump(rows, open(OUT, 'w'))
print('done', len(rows), 'rows', round(time.time()-t0, 1), 's; maxrss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20)
