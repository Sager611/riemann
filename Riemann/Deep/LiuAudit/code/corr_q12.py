#!/usr/bin/env python3
"""Lemma 12.10 (l.6474-6552 of number_fields_v4.tex) over F = Q(zeta12): full correlation
F(u,v;j) = sum_{x mod u, y mod v, vx-uy = j mod uv} chi_u(x) conj(chi_v(y)).
Checks: (a) Fourier form with the full Gauss sums of Lemma 12.9; (b) j=0: 1_{u=v} phi(u);
(c) closed form chi_{n1}(k) conj(chi_{n2}(-k)) L_C(n1,n2;k) with the reciprocity factor rho = chi_{n2}(n1)/chi_{n1}(n2)
    of our (not necessarily allowed) generators inserted (Lemma 3.4: rho = 1 for allowed generators);
(d) complete-support form F(Da,Eb;j) = F(D,E;j) chi_a(j) conj(chi_b(-j)) times the unit factor of l.6537-6541.
"""
import sys, os, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from gauss_q12 import Prime, Ideal, roots_mod_p, mul, pw, neg, ev, tr, star, ONE

def F_bf(u, v, j):
    uv = Ideal(merge(u.exps, v.exps))
    Nu, Nv = u.N, v.N
    y = np.arange(Nv, dtype=np.int64)
    cy = np.conj(v.chi_int(y))
    tot = 0j
    xs = np.arange(Nu, dtype=np.int64)
    cx = u.chi_int(xs)
    for x in range(Nu):
        if cx[x] == 0: continue
        ok = np.ones(Nv, dtype=bool)
        for P, e, m, r in uv.mods:
            a = ev(v.gen, r, m); b = ev(u.gen, r, m); c = ev(j, r, m)
            ok &= ((a*x - b*y - c) % m) == 0
        tot += cx[x]*np.sum(cy[ok])
    return tot

def merge(e1, e2):
    d = {}
    for P, e in list(e1)+list(e2): d[P] = d.get(P, 0)+e
    return [(P, e) for P, e in d.items()]

def gauss_all(u):
    # Ghat[t] = N^{-1/2} sum_x chi_u(x) e(x t / N), t mod N
    N = u.N
    chi = u.chi_int(np.arange(N, dtype=np.int64))
    return np.fft.ifft(chi)*N/math.sqrt(N)

def F_fourier(u, v, j):
    uv = Ideal(merge(u.exps, v.exps))
    N = uv.N
    Gu = gauss_all(u); Gv = gauss_all(v)
    h = np.arange(N, dtype=np.int64)
    Tu = tr(u.ustar) % u.N; Tv = tr(v.ustar) % v.N
    Tj = tr(mul(j, uv.ustar)) % N
    gu = Gu[(h*Tu) % u.N]; gv = Gv[(h*Tv) % v.N]
    ph = np.exp(-2j*np.pi*((h*Tj) % N)/N)
    return np.sum(gu*np.conj(gv)*ph)/math.sqrt(u.N*v.N)

def ideal_part(exps, keep):
    return [(P, e) for P, e in exps if P in keep]

def div_gen(a, I):
    # exact a / I.gen
    b = mul(a, I.ustar)
    assert all(t % I.N == 0 for t in b), 'not divisible'
    return tuple(t//I.N for t in b)

def divisible(a, I):
    b = mul(a, I.ustar)
    return all(t % I.N == 0 for t in b)

def closed_form(u, v, j):
    du = dict(u.exps); dv = dict(v.exps)
    primes = set(du) | set(dv)
    C = Ideal([(P, min(du.get(P, 0), dv.get(P, 0))) for P in primes])
    n1 = Ideal([(P, du.get(P, 0)-min(du.get(P, 0), dv.get(P, 0))) for P in primes])
    n2 = Ideal([(P, dv.get(P, 0)-min(du.get(P, 0), dv.get(P, 0))) for P in primes])
    assert mul(C.gen, n1.gen) == u.gen and mul(C.gen, n2.gen) == v.gen
    if not divisible(j, C): return 0j, 1+0j
    k = div_gen(j, C)
    rho = n2.chi(n1.gen)/n1.chi(n2.gen) if (n1.N > 1 and n2.N > 1) else 1+0j
    val = n1.chi(k)*np.conj(n2.chi(neg(k)))
    for P, c in C.exps:
        Pn = P.p
        pk = ev(k, P.r, Pn) == 0
        u1 = ev(n1.gen, P.r, Pn) != 0; u2 = ev(n2.gen, P.r, Pn) != 0
        if u1 and u2:
            ph = (P.chi(n1.gen)*np.conj(P.chi(n2.gen)))**c
            if pk: L = Pn**(c-1)*ph*(Pn-1)
            elif c % 6: L = Pn**(c-1)*ph*(-1)
            else: L = Pn**(c-1)*ph*(Pn-2)
        else:
            L = Pn**(c-1)*(Pn-1)*(1 if c % 6 == 0 else 0)*(0 if pk else 1)
        val *= L
    return val, rho

def main():
    P13 = Prime(13, roots_mod_p(13)[0]); P37 = Prime(37, roots_mod_p(37)[0]); P61 = Prime(61, roots_mod_p(61)[0]); P73 = Prime(73, roots_mod_p(73)[0])
    I = lambda *e: Ideal(list(e))
    cases = [
        (I((P13, 1)), I((P13, 1))), (I((P13, 1)), I((P13, 2))), (I((P13, 2)), I((P13, 1))), (I((P13, 2)), I((P13, 2))),
        (I((P13, 1), (P37, 1)), I((P13, 1))), (I((P13, 1)), I((P13, 1), (P37, 1))),
        (I((P13, 1), (P37, 1)), I((P13, 1), (P61, 1))), (I((P37, 1)), I((P61, 1))), (I((P37, 1)), I((P13, 1), (P61, 1))),
        (I((P13, 2), (P37, 1)), I((P13, 1), (P61, 1))), (I((P13, 3)), I((P13, 1))), (I((P13, 1), (P37, 1)), I((P13, 1), (P37, 1))),
        (I((P37, 1), (P61, 1)), I((P73, 1))),
    ]
    npass = nfail = 0; recip = []
    rng = np.random.default_rng(1)
    for u, v in cases:
        uv = Ideal(merge(u.exps, v.exps))
        js = [(0, 0, 0, 0), ONE, (2, 0, 0, 0), (0, 1, 0, 0), (1, 1, 1, 0)]
        js += [mul(P.pi, (1, 2, 0, 1)) for P, e in uv.exps]
        js += [mul(mul(P.pi, P.pi), (3, 0, 1, 0)) for P, e in uv.exps]
        js += [uv.gen, tuple(int(t) for t in rng.integers(-50, 50, 4))]
        for j in js:
            Fb = F_bf(u, v, j)
            Ff = F_fourier(u, v, j)
            cf, rho = closed_form(u, v, j)
            okA = abs(Fb-Ff) < 1e-6*max(1, abs(Fb))
            okC = abs(Fb-cf*rho) < 1e-6*max(1, abs(Fb))
            okC0 = abs(Fb-cf) < 1e-6*max(1, abs(Fb))
            if j == (0, 0, 0, 0):
                same = dict(u.exps) == dict(v.exps)
                phi = 1
                for P, e in u.exps: phi *= P.p**(e-1)*(P.p-1)
                okB = abs(Fb-(phi if same else 0)) < 1e-6
            else: okB = True
            if abs(rho-1) > 1e-9: recip.append((u.N, v.N, rho))
            ok = okA and okB and okC
            npass += ok; nfail += (not ok)
            print(('PASS' if ok else 'FAIL'), f'N(u)={u.N} N(v)={v.N} j={j}: F={Fb:.6g} fourier_ok={okA} j0_ok={okB} closed(with rho)_ok={okC} closed(rho=1)_ok={okC0} rho={rho:.4g}', flush=True)
    # (d) complete support: D,E with common prime 13, a at 37, b at 61
    for D, E, a, b in [(I((P13, 1)), I((P13, 1)), I((P37, 1)), I((P61, 1))), (I((P13, 1)), I((P13, 2)), I((P37, 1)), I((P61, 1))),
                       (I((P13, 2)), I((P13, 1)), I((P37, 1)), I((P73, 1)))]:
        Da = Ideal(D.exps+a.exps); Eb = Ideal(E.exps+b.exps)
        for j in [ONE, (0, 1, 0, 0), mul(P13.pi, (1, 1, 0, 0)), mul(P13.pi, (2, 0, 1, 0)), mul(mul(P13.pi, P13.pi), ONE), mul(P37.pi, P13.pi), (5, 3, 0, 2)]:
            lhs = F_bf(Da, Eb, j)
            base = F_bf(D, E, j)*a.chi(j)*np.conj(b.chi(neg(j)))
            unit = (E.chi(a.gen)/a.chi(E.gen))*(b.chi(D.gen)/D.chi(b.gen))*(b.chi(a.gen)/a.chi(b.gen))
            ok1 = abs(lhs-base) < 1e-6*max(1, abs(lhs)); ok2 = abs(lhs-base*unit) < 1e-6*max(1, abs(lhs)); ok3 = abs(lhs-base/unit) < 1e-6*max(1, abs(lhs))
            ok = ok1 or ok2 or ok3
            npass += ok; nfail += (not ok)
            print(('PASS' if ok else 'FAIL'), f'support N(D)={D.N} N(E)={E.N} N(a)={a.N} N(b)={b.N} j={j}: lhs={lhs:.6g} plain={ok1} times_unit={ok2} over_unit={ok3} unit={unit:.4g}', flush=True)
    # c = 6 (principal power) single prime: u = v = pi^6, solve y from x (y = x - k mod p^6)
    P = P13; p = 13; N6 = p**6
    r6 = None
    u6 = Ideal([(P, 6)])
    for k0 in [(1, 0, 0, 0), (2, 1, 0, 0), mul(P.pi, ONE)]:
        j = mul(pw(P.pi, 6), k0)
        # brute: x - y = j/pi^6 = k0 (mod pi^6); in Z/p^6 via root r6
        r = u6.mods[0][3]
        kk = ev(k0, r, N6)
        x = np.arange(N6, dtype=np.int64)
        val = np.sum(u6.chi_int(x)*np.conj(u6.chi_int((x-kk) % N6)))
        cf, rho = closed_form(u6, u6, j)
        ok = abs(val-cf) < 1e-6*max(1, abs(cf))
        npass += ok; nfail += (not ok)
        print(('PASS' if ok else 'FAIL'), f'c=6 u=v=pi13^6 k={k0}: F={val:.8g} closed={cf:.8g} (P^5(P-2)={p**5*(p-2)}, P^5(P-1)={p**5*(p-1)})', flush=True)
    print('reciprocity factors != 1 for our generators:', recip[:6])
    print('SUMMARY_CORR', npass, 'pass', nfail, 'fail')

if __name__ == '__main__':
    main()

