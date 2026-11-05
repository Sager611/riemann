# t5_family.py — character decomposition of the periodic functions in Hikami's strange identities / Lawrence–Zagier:
# torus knots T(s,t) (chi mod 2st: +1 at +-(st-s-t), -1 at +-(st+s-t)) and Brieskorn spheres Sigma(p1,p2,p3) (chi mod 2P,
# n = P(1 + e1/p1 + e2/p2 + e3/p3) -> -e1 e2 e3). A single character <=> Euler product; else Saias–Weingartner => zeros in Re s > 1.
import numpy as np, itertools, math
def factor(n):
    f, p = {}, 2
    while p*p <= n:
        while n % p == 0: f[p] = f.get(p, 0)+1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f
def char_table(N):
    units = [n for n in range(N) if math.gcd(n, N) == 1]
    # generators of (Z/N)^* via CRT components
    gens = []
    for p, e in factor(N).items():
        q = p**e
        gens.append((p, e, q))
    # brute force: build group as product of cyclic factors by searching elements; simpler: characters via discrete logs per prime power
    tables = []
    for p, e, q in gens:
        U = [n for n in range(q) if math.gcd(n, q) == 1]
        if p == 2 and e >= 3:
            # (Z/2^e)^* = <-1> x <5>
            order5 = q//4
            log = {}
            for a in range(2):
                for b in range(order5):
                    log[((-1)**a*pow(5, b, q)) % q] = (a, b)
            tables.append(('2', q, [2, order5], log))
        else:
            phi = len(U)
            g = next(g for g in U if len({pow(g, k, q) for k in range(phi)}) == phi)
            log = {pow(g, k, q): (k,) for k in range(phi)}
            tables.append(('c', q, [phi], log))
    orders = sum((t[2] for t in tables), [])
    chars = []
    for exps in itertools.product(*[range(o) for o in orders]):
        vals = np.zeros(N, complex)
        for n in units:
            ph = 0.0; i = 0
            for t in tables:
                lg = t[3][n % t[1]]
                for j, l in enumerate(lg):
                    ph += exps[i+j]*l/t[2][j]
                i += len(t[2])
            vals[n] = np.exp(2j*np.pi*ph)
        chars.append(vals)
    return units, chars
def conductor(chi, N):
    for d in sorted(d for d in range(1, N+1) if N % d == 0):
        if all(abs(chi[n]-1) < 1e-9 for n in range(N) if math.gcd(n, N) == 1 and n % d == 1 % d): return d
def decompose(fvals, N):
    units, chars = char_table(N)
    out = []
    for c in chars:
        coef = sum(fvals[n]*np.conj(c[n]) for n in units)/len(units)
        if abs(coef) > 1e-9: out.append((round(abs(coef), 6), conductor(c, N), 'even' if abs(c[N-1]-1) < 1e-9 else 'odd'))
    return out
def torus(s, t):
    N = 2*s*t; f = np.zeros(N)
    for r, v in [(s*t-s-t, 1), (s*t+s-t, -1), (s*t-s+t, -1), (s*t+s+t, 1)]: f[r % N] += v
    return N, f
def brieskorn(p):
    P = p[0]*p[1]*p[2]; N = 2*P; f = np.zeros(N)
    for e in itertools.product([1, -1], repeat=3):
        n = P + e[0]*P//p[0] + e[1]*P//p[1] + e[2]*P//p[2]
        f[n % N] += -e[0]*e[1]*e[2]
    return N, f
single = []
for s in range(2, 8):
    for t in range(s+1, 13):
        if math.gcd(s, t) != 1: continue
        N, f = torus(s, t)
        d = decompose(f, N)
        if len(d) == 1: single.append((s, t))
        if (s, t) in [(2, 3), (2, 5), (3, 4), (2, 7), (3, 5), (2, 9), (4, 5)]:
            print(f'T({s},{t}) mod {N}: {len(d)} characters [(|coef|, conductor, parity)]:', d)
print('torus knots T(s,t), 2<=s<8, s<t<13, whose strange-identity function is a single character:', single)
for p in [(2, 3, 5), (2, 3, 7), (2, 5, 7), (3, 4, 5)]:
    N, f = brieskorn(p)
    d = decompose(f, N)
    print(f'Sigma{p} mod {N}: support {sorted(int(n) for n in np.nonzero(f)[0])[:10]}..., {len(d)} characters:', d[:8])
