#!/usr/bin/env python3
"""Falsification tests over F = Q(zeta12) for Liu arXiv 2610.26835v2 (number_fields_v4.tex):
 Lemma 12.9  (l.6435-6472)  full Gauss sums incl. imprimitive powers + CRT
 Lemma 12.8  (l.6355-6363)  local Gauss table |P^{-i/2} G(p^i, h p^6)|^2 = 0, (1-1/P)^2, 1/P at i = 1, 6, 7
 Lemma 12.10 (l.6474-6552)  full correlation F(u,v;j): Fourier form, j=0, closed form via C=(u,v), complete support.
Exact field arithmetic in Z[zeta12] (zeta^4 = zeta^2 - 1); numerics in float64/complex128.
Prime ideals of degree one over distinct rational primes p = 1 mod 12, so O/u = Z/N(u) with integer representatives.
chi_pi(x) = j6(lift of xbar^{(p-1)/6}), j(zeta12) = exp(2 pi i/12)  (l.1884-1886). e(z) = exp(2 pi i Tr z).
The reciprocity factor chi_b(a)/chi_a(b) for the chosen (not necessarily allowed) generators is printed separately:
Lemma 3.4 makes it 1 for allowed generators; the identities are tested with that factor inserted explicitly.
"""
import sys, itertools, cmath, math
import numpy as np

def mul(a, b):
    c = [0]*7
    for i in range(4):
        if a[i]:
            for j in range(4):
                c[i+j] += a[i]*b[j]
    # zeta^6 = -1, zeta^5 = zeta^3 - zeta, zeta^4 = zeta^2 - 1
    for k in (6, 5, 4):
        t = c[k]; c[k] = 0
        if k == 6: c[0] -= t
        elif k == 5: c[3] += t; c[1] -= t
        else: c[2] += t; c[0] -= t
    return tuple(c[:4])

def add(a, b): return tuple(x+y for x, y in zip(a, b))
def neg(a): return tuple(-x for x in a)
def sc(k, a): return tuple(k*x for x in a)
ONE = (1, 0, 0, 0)
def pw(a, e):
    r = ONE
    for _ in range(e): r = mul(r, a)
    return r

def sigma(a, k):
    # zeta -> zeta^k, k in {1,5,7,11}
    z = (0, 1, 0, 0); zk = pw(z, k)
    r = (0, 0, 0, 0); zp = ONE
    for i in range(4):
        r = add(r, sc(a[i], zp)); zp = mul(zp, zk)
    return r

def star(a):  # product of the other three conjugates: a * star(a) = N(a)
    return mul(mul(sigma(a, 5), sigma(a, 7)), sigma(a, 11))

def norm(a):
    n = mul(a, star(a))
    assert n[1] == n[2] == n[3] == 0, n
    return n[0]

def tr(a): return 4*a[0] + 2*a[2]

def ev(a, r, m): return (a[0] + a[1]*r + a[2]*r*r + a[3]*r**3) % m

def phi12(x, m): return (x**4 - x**2 + 1) % m

def roots_mod_p(p): return [r for r in range(p) if phi12(r, p) == 0]

def hensel(r, p, e):
    m = p; x = r
    for _ in range(e):
        m = min(m*p, p**e) if m < p**e else m
    x = r; mod = p
    while mod < p**e:
        mod = min(mod*mod, p**e)
        f = (x**4 - x**2 + 1) % mod; d = (4*x**3 - 2*x) % mod
        x = (x - f*pow(d, -1, mod)) % mod
    assert phi12(x, p**e) == 0
    return x

def find_gen(p, r):
    rng = range(-4, 5)
    best = None
    for a in itertools.product(rng, repeat=4):
        if ev(a, r, p) == 0 and norm(a) == p:
            s = sum(abs(t) for t in a)
            if best is None or s < best[0]: best = (s, a)
    assert best, (p, r)
    return best[1]

ZETA = [cmath.exp(2j*math.pi*k/6) for k in range(6)]

class Prime:
    def __init__(self, p, r):
        self.p, self.r = p, r
        self.pi = find_gen(p, r)
        self.pistar = star(self.pi)
        self.z6 = [pow(r, 2*k, p) for k in range(6)]
        e = (p-1)//6
        tab = np.zeros(p, dtype=complex)
        for x in range(1, p):
            t = pow(x, e, p); tab[x] = ZETA[self.z6.index(t)]
        self.tab = tab  # chi_pi on integers mod p
    def chi(self, a):  # chi_pi(a), a element of O
        x = ev(a, self.r, self.p)
        return complex(self.tab[x])
    def val(self, a):  # v_pi(a) for nonzero element a
        v = 0
        while True:
            if ev(a, self.r, self.p) != 0: return v
            a = div_pi(a, self); v += 1

def div_pi(a, P):
    # exact division a / pi = a * pistar / p
    b = mul(a, P.pistar)
    assert all(t % P.p == 0 for t in b), (a, P.p)
    return tuple(t//P.p for t in b)

class Ideal:
    """u = prod pi_l^{e_l} with generator prod of the chosen pi_l (times an optional unit)."""
    def __init__(self, exps, unit=ONE):
        self.exps = [(P, e) for P, e in exps if e > 0]
        g = unit
        for P, e in self.exps: g = mul(g, pw(P.pi, e))
        self.gen = g
        self.N = 1
        for P, e in self.exps: self.N *= P.p**e
        self.ustar = star(g)
        self.mods = [(P, e, P.p**e, hensel(P.r, P.p, e)) for P, e in self.exps]
    def chi_int(self, x):  # chi_u on integer array x
        out = np.ones(x.shape, dtype=complex)
        for P, e in self.exps:
            out *= P.tab[x % P.p]**e
        return out
    def chi(self, a):
        out = 1+0j
        for P, e in self.exps: out *= P.chi(a)**e
        return out

def gauss_bf(u, h):
    """G(u,h) = q_u^{-1/2} sum_{x mod u} chi_u(x) e(hx/u), brute force over integer reps, chunked."""
    N = u.N
    t = tr(mul(h, u.ustar)) % N   # Tr(h x u*) = x Tr(h u*) for integer x
    tot = 0j; CH = 1 << 21
    for s in range(0, N, CH):
        x = np.arange(s, min(N, s+CH), dtype=np.int64)
        ph = (x * t) % N
        tot += np.sum(u.chi_int(x) * np.exp(2j*np.pi*ph/N))
    return tot / math.sqrt(N)

def main():
    log = []
    def rep(name, ok, info=''):
        log.append((name, ok)); print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''), flush=True)
    P13 = Prime(13, roots_mod_p(13)[0]); P37 = Prime(37, roots_mod_p(37)[0]); P61 = Prime(61, roots_mod_p(61)[0])
    for P in (P13, P37, P61): print('prime', P.p, 'root', P.r, 'gen', P.pi, 'norm', norm(P.pi))
    # sanity: trace formula, chi multiplicative, chi(-1) = 1
    assert tr(ONE) == 4 and tr((0, 0, 1, 0)) == 2 and tr((0, 1, 0, 0)) == 0 and tr((0, 0, 0, 1)) == 0
    for P in (P13, P37):
        assert abs(P.chi(neg(ONE)) - 1) < 1e-12
    tol = 1e-7
    # ---------- Lemma 12.9 single prime powers ----------
    for P, imax in ((P13, 7), (P37, 4)):
        p = P.p
        gam = {}
        for i in range(1, imax+1):
            u = Ideal([(P, i)])
            hs = []
            for m in range(0, i+2):
                for h0 in ((1, 0, 0, 0), (2, 1, 0, 0), (0, 1, 1, 0), (3, 0, -1, 2)):
                    if ev(h0, P.r, p) == 0: continue
                    hs.append((m, h0, mul(pw(P.pi, m), h0)))
            hs.append((10**9, None, (0, 0, 0, 0)))
            for m, h0, h in hs:
                G = gauss_bf(u, h)
                vh = m
                if i % 6:
                    mag = p**((i-1)/2) * (1 if vh == i-1 else 0)
                    ok = abs(abs(G) - mag) < tol*max(1, mag)
                    if ok and vh == i-1:
                        # phase gamma_i(p) * conj(chi_p(h/p^{i-1}))^i ; gamma_i = P^{-1/2} sum chi^i(x0) e(x0/pi)
                        if i not in gam:
                            x0 = np.arange(p, dtype=np.int64)
                            gam[i] = np.sum(P.tab[x0]**i * np.exp(2j*np.pi*((x0*tr(P.pistar)) % p)/p))/math.sqrt(p)
                        pred = mag * gam[i] * np.conj(P.chi(h0))**i
                        ok = abs(G - pred) < tol*mag
                    rep(f'L12.9 p={p} i={i} v(h)={vh if vh < 10**8 else "inf"} h0={h0}', ok, f'G={G:.6g}')
                else:
                    pred = p**(i/2)*(1 if vh >= i else 0) - p**(i/2-1)*(1 if vh >= i-1 else 0)
                    ok = abs(G - pred) < tol*max(1, abs(pred))
                    rep(f'L12.9 p={p} i={i} v(h)={vh if vh < 10**8 else "inf"} h0={h0}', ok, f'G={G:.6g} pred={pred:.6g}')
    # ---------- Lemma 12.8 local table at P=13 ----------
    p = 13; P = P13
    for h0 in ((1, 0, 0, 0), (2, 1, 0, 0)):
        h = mul(pw(P.pi, 6), h0)
        vals = []
        for i in (1, 6, 7):
            G = gauss_bf(Ideal([(P, i)]), h); vals.append(abs(p**(-i/2)*G)**2)
        pred = [0, (1-1/p)**2, 1/p]
        # the table is printed for P^{-i/2} G with G already normalized by q_u^{-1/2}: compare both readings
        G1 = [abs(gauss_bf(Ideal([(P, i)]), h))**2 for i in (1, 6, 7)]
        predG = [0, (p**3-p**2)**2, p**6]
        ok = all(abs(a-b) < 1e-6*max(1, b) for a, b in zip(G1, predG))
        rep(f'L12.8 table p=13 h0={h0}: |G(p^i,h p^6)|^2 at i=1,6,7 = 0,(P^3-P^2)^2,P^6 i.e. |P^-i/2 G|^2 reading {[round(x,9) for x in vals]}', ok,
            'with Lemma 12.9 normalisation the printed table holds for P^{-i/2} * (unnormalised sum) = G')
        # original frequency only i=1 with |.|^2 = P^{-1} in the printed normalisation
    # ---------- Lemma 12.9 CRT ----------
    pairs = [(Ideal([(P13, 1)]), Ideal([(P37, 1)])), (Ideal([(P13, 2)]), Ideal([(P37, 1)])), (Ideal([(P13, 1)]), Ideal([(P37, 2)])),
             (Ideal([(P13, 1)]), Ideal([(P37, 1), (P61, 1)])), (Ideal([(P13, 3)]), Ideal([(P37, 1)])), (Ideal([(P13, 6)]), Ideal([(P37, 1)]))]
    for a, b in pairs:
        ab = Ideal(a.exps + b.exps)
        assert ab.gen == mul(a.gen, b.gen)
        for h in ((1, 0, 0, 0), (2, 1, 0, 0), (0, 0, 0, 0), mul(a.gen, (1, 1, 0, 0)), mul(P13.pi, (0, 1, 0, 1))):
            Gab = gauss_bf(ab, h)
            pred = a.chi(b.gen)*b.chi(a.gen)*gauss_bf(a, h)*gauss_bf(b, h)
            ok = abs(Gab-pred) < 1e-6*max(1, abs(pred))
            rep(f'L12.9 CRT N(a)={a.N} N(b)={b.N} h={h}', ok, f'|G(ab)|={abs(Gab):.5g}')
    print('SUMMARY_GAUSS', sum(ok for _, ok in log), 'pass of', len(log))
    return log

if __name__ == '__main__':
    main()
