# Exact finite check of AskMainBalance (colour_balance, StronglyIrreducible.lean:695;
# definitions FirstPrimeNeSecondPrime.lean:61-101, PartnersToFactorsAt.lean:13-29, StatementOC.lean:11-13).
# Values live in Q(i) tensor (formal symbols): log p, gamma, log a, log w_rho, ell_psi. Exact Fractions, no floats.
import random, sys
from fractions import Fraction as Fr
from itertools import product

def cadd(x, y): return (x[0]+y[0], x[1]+y[1])
def cmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def cconjg(x): return (x[0], -x[1])
ZERO = (Fr(0), Fr(0)); ONE = (Fr(1), Fr(0))

class E(dict):  # formal linear combination: key -> Gaussian rational
    def add(self, k, c):
        v = cadd(self.get(k, ZERO), c)
        if v == ZERO: self.pop(k, None)
        else: self[k] = v
        return self
    def addE(self, o, c=ONE):
        for k, v in o.items(): self.add(k, cmul(v, c))
        return self
def num(q): return (Fr(q), Fr(0))

def primes_upto(N): return [p for p in range(2, N+1) if all(p % r for r in range(2, int(p**0.5)+1))]
def factor(n):
    f = {}; p = 2
    while p*p <= n:
        while n % p == 0: f[p] = f.get(p, 0)+1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f

def run(nH, seed, break_host=False, wrong_main=False):
    rnd = random.Random(seed)
    G = [((e1, e2), h) for e1 in (0, 1) for e2 in (0, 1) for h in range(nH)]
    gadd = lambda g, k: (((g[0][0]+k[0][0]) % 2, (g[0][1]+k[0][1]) % 2), (g[1]+k[1]) % nH)
    gsmul = lambda i, g: (((i*g[0][0]) % 2, (i*g[0][1]) % 2), (i*g[1]) % nH)
    gsub = lambda g, k: gadd(g, gsmul(-1 % (2*nH), k)) if False else (((g[0][0]-k[0][0]) % 2, (g[0][1]-k[0][1]) % 2), (g[1]-k[1]) % nH)
    inPhi = lambda g: g[0][0] == 0
    Phi = [g for g in G if inPhi(g)]
    assert nH in (1, 2, 4)
    zeta = {1: [ONE], 2: [ONE, num(-1)], 4: [ONE, (Fr(0), Fr(1)), num(-1), (Fr(0), Fr(-1))]}[nH]
    chars = [(e1, e2, h) for e1 in (0, 1) for e2 in (0, 1) for h in range(nH)]
    def psi(c, g):
        s = num((-1)**((c[0]*g[0][0] + c[1]*g[0][1]) % 2))
        return cmul(s, zeta[(c[2]*g[1]) % nH])
    chC = (1, 0, 0); triv = (0, 0, 0)
    cconj = ((1, 1), 0); sgc0 = ((1, 0), 0)
    partner = lambda t, g: gadd(gadd(cconj, g), ((0, 0), t))
    # ArithData + HostFacts: q1 = 3, d = 7 (odd, distinct primes)
    q1, d = 3, 7
    a = Fr(rnd.randint(40, 90), rnd.choice([1, 2, 3]))
    w = {g: Fr(rnd.randint(10, 30), 10) for g in G}
    N = int(max(a / w[r] for r in Phi)) + 1
    P = primes_upto(N + 5)
    frob = {}; inertia = {}
    for p in P:
        if p in (2, q1):
            frob[p] = rnd.choice(Phi)
            inertia[p] = [((0, 0), 0)] if p == q1 else [((0, 0), 0), ((0, 1), 0)]  # inside Phi
        elif p == d:
            frob[p] = rnd.choice(G); inertia[p] = [((0, 0), 0), sgc0]
        else:
            frob[p] = rnd.choice(G); inertia[p] = [((0, 0), 0)]
    if break_host:  # violate HostFacts.unram at p = 5: nontrivial inertia at an unramified prime
        inertia[5] = [((0, 0), 0), ((1, 1), 0)]
    def cosetMean(I, h, c):
        s = ZERO
        for x in I: s = cadd(s, psi(c, gadd(h, x)))
        return (s[0]/len(I), s[1]/len(I))
    def primVal(c, n):
        v = ONE
        for p, e in factor(n).items(): v = cmul(v, cosetMean(inertia[p], gsmul(e, frob[p]), c))
        return v
    def Lam(n):
        f = factor(n)
        return ('log', list(f)[0]) if len(f) == 1 else None
    def ramp_terms(c, y, scale):  # scale * sum_{n<=y} Lam(n) primVal(c,n) (y/n - 1)
        out = E()
        for n in range(2, int(y) + 1):
            k = Lam(n)
            if k: out.add(k, cmul(cmul(primVal(c, n), num(y / n - 1)), num(scale)))
        return out
    def GPsi_times_w(c, wv):  # w * GPsi(c, a/w), exact; Gz for c = 0
        y = a / wv
        if c == triv:
            out = ramp_terms(c, y, wv)
            # - y log y + (1+gamma) y, times w: -a(log a - log w) + (1+gamma) a
            out.add(('loga',), num(-a)); out.add(('logw', wv), num(a)); out.add(('1',), num(a)); out.add(('gamma',), num(a))
            return out
        out = ramp_terms(c, y, wv)
        out.add(('ell', c), num(a))
        return out
    # ASK
    t = rnd.randrange(nH)
    ASK = E()
    for p in P:
        if p > N: continue
        if p not in (2, q1, d) and not inPhi(frob[p]):
            for i in range(1, N + 1):
                g = gsmul(i, frob[p])
                if not inPhi(g):
                    v = a / Fr(p)**i - w[partner(t, g)]
                    if v > 0: ASK.add(('log', p), num(2 * v))
    for i in range(1, N + 1):
        for g in G:
            if gsub(g, gsmul(i, frob[d])) in inertia[d] and not inPhi(g):
                v = a / Fr(d)**i - w[partner(t, g)]
                if v > 0: ASK.add(('log', d), num(v))
    MAIN = E()
    MAIN.add(('loga',), num(a)); MAIN.add(('1',), num(-a)); MAIN.add(('gamma',), num(-a if not wrong_main else a))
    MAIN.add(('ell', chC), num(a))
    for r in Phi: MAIN.add(('logw', w[r]), num(-a / len(Phi)))
    De = E()
    for g in G:
        if inPhi(g): continue
        wv = w[partner(t, g)]
        for c in chars:
            De.addE(GPsi_times_w(c, wv), cmul(cconjg(psi(c, g)), num(Fr(1, len(Phi)))))
    diff = E().addE(ASK).addE(MAIN, num(-1)).addE(De, num(-1))
    return len(diff), len(ASK), len(De), N, a, t

if __name__ == '__main__':
    tot = 0
    for nH in (1, 2, 4):
        for seed in range(6):
            nd, na, nde, N, a, t = run(nH, seed)
            tot += 1
            print(f"H=Z/{nH} seed={seed} a={a} N={N} t={t}: nonzero keys in ASK-MAIN-De = {nd} (ASK keys {na}, De keys {nde})")
            assert nd == 0
    print(f"PASS: {tot} exact instances, ASK - MAIN = De identically in all formal symbols")
    for nH in (2, 4):
        nd = run(nH, 1, break_host=True)[0]; print(f"control H=Z/{nH}: HostFacts.unram violated at p=5 -> nonzero keys {nd}")
        nd = run(nH, 1, wrong_main=True)[0]; print(f"control H=Z/{nH}: MAIN with +gamma instead of -gamma -> nonzero keys {nd}")
