#!/usr/bin/env python3
"""Independent checks for the Kopp-Lagarias Problem 1.7 certificate (cert.py).  Run:
  uv run --no-project --with numpy --with python-flint python3 checks.py
C1  Lemma A brute force: no primitive alpha with 1 < |N alpha| < sqrt(Delta)/2 in the normalized window (all n <= 3000, 40 random n <= 10^9 per family)
C2  principal cycle = {(1,n-2,2-n),(2-n,n-2,1)} (N=+1) or {(1,n,-1),(-1,n,1)} (N=-1), all n <= 20000
C3  Cor. B mechanism: for split p with 4p^4 < Delta, the reduced cycle of p^2 (composition) contains no form with |a| = 1 (wide non-principal); n <= 2000 and 30 random n <= 10^9 per family
C4  genus theory: #{C in Cl+ : C^2 = 1} = 2^(mu(Delta)-1) (Halter-Koch Thm 6.5.2) for every n <= 1000
C5  sieve implementation: brute-force split-prime test (Jacobi-symbol code, independent of cert.py) over all n in [1001, 3*10^6] = cert.py survivors in that range
C6  every 2-torsion order found has no split prime p with 4p^4 < Delta (consistency with Cor. B), Delta >= 100
C7  writes pari_checks.gp (class number formula with L(1,chi), quadclassunit, quadunit) for gp
C8  Lemma A analytic inequalities on a grid of Delta in [100, 10^40] (arb)
"""
import sys, os, random, time
from math import isqrt, gcd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from forms import class_data, compose, reduce_form, rho, is_reduced, factor
import flint
from flint import arb
random.seed(20261106)
T0 = time.time()
LOG = []
def out(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.append(s)

def jacobi(a, n):
    """Jacobi symbol (a/n), n odd positive (binary algorithm, independent of Euler's criterion)."""
    assert n > 0 and n % 2 == 1
    a %= n; t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: t = -t
        a %= n
    return t if n == 1 else 0
def primes_upto(m):
    s = bytearray([1]) * (m + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(m) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(m + 1) if s[i]]
PR = primes_upto(40000)
def splits(D, p):
    if p == 2: return D % 8 == 1
    return jacobi(D, p) == 1
def passes(n, N):
    D = n * n - 4 * N
    for p in PR:
        if 4 * p ** 4 >= D: return True
        if splits(D, p): return False
    raise RuntimeError

def fams(nmax):
    for N in (1, -1):
        for n in range(1, nmax + 1):
            if N == 1 and n < 3: continue
            yield n, N

# C1 -------------------------------------------------------------------
def lemmaA_brute(n, N):
    """search the window of the proof of Lemma A: y >= 1, y < (eps+1) sqrt(m)/sqrt(D) <= (eps+1) (D/4)^(1/4)/sqrt(D) (+1 slack), |x - y sqrt D| < 3."""
    D = n * n - 4 * N
    s = isqrt(D)
    # eps + 1 < sqrt(D) + 2 ; sqrt(m) < (D/4)^(1/4)
    ymax = int(((s + 3) * (isqrt(isqrt(D // 4)) + 2)) // s) + 2
    bad = []
    for y in range(1, ymax + 1):
        xc = isqrt(D * y * y)
        for x in range(xc - 3, xc + 5):
            if x < 0 or (x - y * D) % 2: continue
            v = x * x - D * y * y
            if v % 4: continue
            m = abs(v) // 4
            if m <= 1 or 4 * m * m >= D: continue
            u = (x - y * D) // 2
            if gcd(u, y) == 1:
                bad.append((x, y, m))
    return bad
cnt = 0; bad = []
for n, N in fams(3000):
    if n * n - 4 * N < 100: continue
    b = lemmaA_brute(n, N); cnt += 1
    if b: bad.append((n, N, b[:3]))
for N in (1, -1):
    for _ in range(40):
        n = random.randint(3001, 10 ** 9)
        b = lemmaA_brute(n, N); cnt += 1
        if b: bad.append((n, N, b[:3]))
out(f'C1 Lemma A brute force: {cnt} orders, counterexamples: {bad[:5]}  [{time.time()-T0:.1f}s]')
assert not bad

# C2 -------------------------------------------------------------------
badc = []
for n, N in fams(20000):
    D = n * n - 4 * N
    if D <= 8: continue
    s = isqrt(D)
    b0 = s if (s - D) % 2 == 0 else s - 1
    f = (1, b0, (b0 * b0 - D) // 4)
    assert is_reduced(f, s)
    cyc = [f]; g = rho(f, D, s)
    while g != f:
        cyc.append(g); g = rho(g, D, s)
        if len(cyc) > 10: break
    want = {(1, n - 2, 2 - n), (2 - n, n - 2, 1)} if N == 1 else {(1, n, -1), (-1, n, 1)}
    if set(cyc) != want and not (N == 1 and n == 3):
        badc.append((n, N, cyc[:4]))
out(f'C2 principal cycles: n <= 20000 both families, mismatches: {badc[:5]}  [{time.time()-T0:.1f}s]')
assert not badc

# C3 -------------------------------------------------------------------
def sqrt_mod_p(a, p):
    """Tonelli-Shanks square root of a quadratic residue a mod odd prime p (result verified)."""
    a %= p
    if p % 4 == 3:
        r = pow(a, (p + 1) // 4, p)
    else:
        q, s = p - 1, 0
        while q % 2 == 0: q //= 2; s += 1
        z = 2
        while pow(z, (p - 1) // 2, p) != p - 1: z += 1
        m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
        while t != 1:
            i, t2 = 0, t
            while t2 != 1: t2 = t2 * t2 % p; i += 1
            b = pow(c, 1 << (m - i - 1), p)
            m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    assert (r * r - a) % p == 0
    return r
def sqrt_mod_4p(D, p):
    """b with b^2 = D (mod 4p), b = D (mod 2), 0 <= b < 2p (p odd prime, (D/p) = 1)."""
    r = sqrt_mod_p(D, p)
    for b in (r, p - r, r + p, 2 * p - r):
        if 0 <= b < 2 * p and (b * b - D) % (4 * p) == 0:
            return b
    raise RuntimeError
def check_corB(n, N, pmax=None):
    D = n * n - 4 * N
    s = isqrt(D)
    res = 0; viol = []
    for p in PR[1:]:
        if 4 * p ** 4 >= D: break
        if not splits(D, p): continue
        b = sqrt_mod_4p(D, p)
        F = (p, b, (b * b - D) // (4 * p))
        S = reduce_form(compose(F, F, D), D, s)
        g = S; seen = 0; principal = False
        while True:
            if abs(g[0]) == 1: principal = True
            g = rho(g, D, s); seen += 1
            if g == S or seen > 10 ** 5: break
        assert seen <= 10 ** 5
        res += 1
        if principal: viol.append(p)
    return res, viol
tot = 0; viols = []
for n, N in fams(2000):
    if n * n - 4 * N < 100: continue
    r, v = check_corB(n, N); tot += r
    if v: viols.append((n, N, v))
for N in (1, -1):
    for _ in range(30):
        n = random.randint(2001, 10 ** 9)
        r, v = check_corB(n, N); tot += r
        if v: viols.append((n, N, v))
out(f'C3 Cor. B mechanism: {tot} (n, p) pairs with p split, 4p^4 < Delta; p^2 wide-principal in: {viols[:5]}  [{time.time()-T0:.1f}s]')
assert not viols
# sharpness probe: split p with Delta/4 <= p^4 < Delta^2 ... count p^2 principal (allowed beyond the bound)
sharp = []
for n, N in fams(400):
    D = n * n - 4 * N
    if D < 100: continue
    s = isqrt(D)
    for p in PR[1:]:
        if p * p > s: break
        if 4 * p ** 4 < D or not splits(D, p): continue
        b = sqrt_mod_4p(D, p)
        F = (p, b, (b * b - D) // (4 * p))
        S = reduce_form(compose(F, F, D), D, s)
        g = S; hit = False
        while True:
            if abs(g[0]) == 1: hit = True
            g = rho(g, D, s)
            if g == S: break
        if hit: sharp.append((n, N, p))
out(f'C3b sharpness: split p with sqrt(Delta)/2 <= p^2 < sqrt(Delta) and p^2 principal, n <= 400: {len(sharp)} cases, e.g. {sharp[:6]}')

# C4 -------------------------------------------------------------------
def mu(D):
    odd = D
    while odd % 2 == 0: odd //= 2
    r = len(factor(odd))
    if D % 4 == 1 or D % 16 == 4: return r
    if D % 16 in (8, 12) or D % 32 == 16: return r + 1
    if D % 32 == 0: return r + 2
    raise ValueError(D)
badg = []; n2 = 0
for n, N in fams(1000):
    D = n * n - 4 * N
    d = class_data(D)
    s = isqrt(D)
    amb = sum(1 for f in d['reps'] if d['cyc'][reduce_form(compose(f, f, D), D, s)] == d['pid'])
    if amb != 2 ** (mu(D) - 1): badg.append((n, N, amb, mu(D)))
    if d['wide2']:
        n2 += 1
        assert d['h'] <= 2 ** (mu(D) - 1)
out(f'C4 genus count |Cl+[2]| = 2^(mu-1): {2 * 1000 - 2} orders, mismatches {badg[:5]}; 2-torsion orders satisfy h <= 2^(mu-1): {n2}  [{time.time()-T0:.1f}s]')
assert not badg

# C5 -------------------------------------------------------------------
surv = {1: set(), -1: set()}
NEX = 1000
for line in open(os.path.join(HERE, 'logs', 'survivors.txt')):
    if line.startswith('N='):
        parts = dict(kv.split('=') for kv in line.split()[:2])
        surv[int(parts['N'])].add(int(parts['n']))
LIM = 3 * 10 ** 6
for N in (1, -1):
    brute = {n for n in range(NEX + 1, LIM + 1) if passes(n, N)}
    mine = {n for n in surv[N] if n <= LIM}
    out(f'C5 N={N:+d}: brute-force survivors in [{NEX+1}, {LIM}]: {len(brute)}; cert.py survivors there: {len(mine)}; equal: {brute == mine}  [{time.time()-T0:.1f}s]')
    assert brute == mine

# C6 -------------------------------------------------------------------
listed = []
for line in open(os.path.join(HERE, 'logs', 'cert_output.txt')):
    line = line.strip()
    if line.startswith('n=') and 'Delta=' in line:
        listed.append(int(__import__('re').search(r'Delta=\s*(\d+)', line).group(1)))
fam_of = {}
for n, N in fams(400):
    fam_of.setdefault(n * n - 4 * N, []).append((n, N))
badl = [D for D in listed if D >= 100 and not all(passes(n, N) for n, N in fam_of[D])]
out(f'C6 the {len(listed)} listed entries with Delta >= 100 all pass the split-prime test: {not badl}')
assert not badl and len(listed) == 90

# C7 -------------------------------------------------------------------
recs = []
for n, N in fams(NEX):
    recs.append((n, N))
for N in (1, -1):
    for n in sorted(surv[N]): recs.append((n, N))
def fund(D):
    f = 1
    for p, e in factor(D).items(): f *= p ** (e // 2)
    d0 = D // (f * f)
    if d0 % 4 not in (0, 1): f //= 2; d0 *= 4
    return f, d0
hs = {}
for n, N in recs:
    D = n * n - 4 * N
    if (n, N) not in hs:
        d = class_data(D)
        hs[(n, N)] = (d['h'], 1 if d['wide2'] else 0)
with open(os.path.join(HERE, 'pari_checks.gp'), 'w') as fh:
    fh.write('\\\\ generated by checks.py: rows [n, N, Delta, f, Delta0, h_exact, wide2_exact]\n')
    fh.write('default(parisizemax, 200000000); default(realprecision, 38);\n')
    fh.write('R = [' + ';'.join(f'{n},{N},{n*n-4*N},{fund(n*n-4*N)[0]},{fund(n*n-4*N)[1]},{hs[(n,N)][0]},{hs[(n,N)][1]}' for n, N in recs) + '];\n')
    fh.write(open(os.path.join(HERE, 'pari_checks_body.gp')).read())
out(f'C7 wrote pari_checks.gp with {len(recs)} orders  [{time.time()-T0:.1f}s]')

# C8 -------------------------------------------------------------------
flint.ctx.prec = 128
worst1 = worst2 = None
for e10 in range(200, 4001):
    D = arb(10) ** (arb(e10) / 100)          # Delta = 10^(e10/100), 100 .. 10^40
    sD = D.sqrt(); q = D.sqrt().sqrt()
    phi1 = sD / 2 - 1 - 2 / sD - q / arb(2).sqrt()                       # lower bound for phi(1) (DERIVATION 2.2)
    phi2 = sD - 4 / sD - (q / arb(2).sqrt()) * (1 + 2 / sD) ** 2 - q / arb(2).sqrt() - 1   # lower bound for phi(Y)/sqrt(m)
    worst1 = phi1 if worst1 is None or phi1 < worst1 else worst1
    worst2 = phi2 if worst2 is None or phi2 < worst2 else worst2
out(f'C8 Lemma A endpoint bounds on Delta in [100, 1e40]: min phi(1) lower bound {worst1.str(5)}, min phi(Y) bound {worst2.str(5)} (both must be > 0)')
assert worst1 > 0 and worst2 > 0
out(f'all checks passed  [{time.time()-T0:.1f}s]')
with open(os.path.join(HERE, 'logs', 'checks_output.txt'), 'w') as fh:
    fh.write('\n'.join(LOG) + '\n')
