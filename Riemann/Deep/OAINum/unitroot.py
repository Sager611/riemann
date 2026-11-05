
# unitroot.py -- Lemma ff:unitroot (root extraction from a primary generator) over k = F_p
from ffpoly import *

def ddf(H, p):
    # distinct-degree factorization of monic square-free H over F_p: dict d -> H_d
    out = {}
    Hr = list(H)
    x = [0, 1]
    xp = x
    d = 0
    while deg(Hr) > 0:
        d += 1
        xp = ppowmod(xp, p, Hr, p)  # X^{p^d} mod Hr
        g = pgcd(psub(xp, x, p), Hr, p)
        if deg(g) > 0:
            out[d] = g
            Hr = pdiv_exact(Hr, g, p)
            xp = pmod(xp, Hr, p) if deg(Hr) > 0 else []
        if 2 * (d + 1) > deg(Hr) and deg(Hr) > 0:
            out[deg(Hr)] = Hr if deg(Hr) not in out else pmul(out[deg(Hr)], Hr, p)
            break
    return out

def crt2(r1, m1, r2, m2, p):
    # r = r1 mod m1, r2 mod m2 (coprime)
    inv = pinvmod(m1, m2, p)
    t = pmod(pmul(psub(r2, r1, p), inv, p), m2, p)
    return pmod(padd(r1, pmul(m1, t, p), p), pmul(m1, m2, p), p), pmul(m1, m2, p)

def frob_matrix(Hd, e, p):
    # matrix of v -> v^(p^e) mod Hd on basis 1..X^{n-1}; columns = images
    n = deg(Hd)
    xq = ppowmod([0, 1], p ** e, Hd, p)
    cols = []
    cur = [1]
    for j in range(n):
        cols.append(cur + [0] * (n - len(cur)))
        cur = pmod(pmul(cur, xq, p), Hd, p)
    # image of X^j is (X^{p^e})^j since Frobenius is a ring hom fixing F_p
    import numpy as np
    return np.array(cols, dtype=np.int64).T % p

def unitroot(a, H, q, zeta, omega, p, stats=None):
    """Return b with b^q = a mod H, a unit promised to be a q-th power. k = F_p, Q = p."""
    H = monic(H, p)
    h = deg(H)
    a = pmod(a, H, p)
    Q = p
    pieces_all = []  # list of (root mod m, m)
    for d, Hd in sorted(ddf(H, p).items()):
        n = Q ** d - 1
        s = 0
        t = n
        while t % q == 0:
            t //= q
            s += 1
        assert s >= 1
        if d % q != 0:
            u = [omega % p]
        else:
            Fm = frob_matrix(Hd, d // q, p)
            nn = deg(Hd)
            import numpy as np
            Amat = (Fm - zeta * np.eye(nn, dtype=np.int64)) % p
            ker = nullspace(Amat, nn, p)
            r = len(ker)
            assert r >= 1
            found = None
            for c in range(0, h * h + 1):
                v = [0] * nn
                cj = 1
                for j in range(r):
                    for i in range(nn):
                        v[i] = (v[i] + cj * int(ker[j][i])) % p
                    cj = (cj * c) % p
                v = norm(v)
                if v and pgcd(v, Hd, p) == [1]:
                    found = v
                    if stats is not None:
                        stats['search_c_max'] = max(stats.get('search_c_max', 0), c)
                    break
            assert found is not None, "unit search failed"
            u = ppowmod(found, t, Hd, p)
        ad = pmod(a, Hd, p)
        qs = q ** s
        # exponents for primary / prime-to-q parts
        e_q = t * pow(t, -1, qs) if qs > 1 else 0
        a_q = ppowmod(ad, e_q, Hd, p)
        if t > 1:
            e_t = qs * pow(qs, -1, t)
            a_t = ppowmod(ad, e_t, Hd, p)
            root_t = ppowmod(a_t, pow(q, -1, t), Hd, p)
        else:
            root_t = [1]
        # digit logarithm with partition
        uinv = pinvmod(u, Hd, p)
        pieces = [(Hd, 0)]
        for j in range(s):
            newp = []
            for (Hp, ell) in pieces:
                z = pmod(pmul(pmod(a_q, Hp, p), ppowmod(pmod(uinv, Hp, p), ell, Hp, p), p), Hp, p)
                z = ppowmod(z, q ** (s - j - 1), Hp, p)
                zeta1 = ppowmod(pmod(u, Hp, p), q ** (s - 1), Hp, p)
                rest = Hp
                zc = [1]
                for c in range(q):
                    if deg(rest) == 0:
                        break
                    g = pgcd(psub(pmod(z, rest, p), pmod(zc, rest, p), p), rest, p) if True else None
                    if deg(g) > 0:
                        newp.append((g, ell + c * q ** j))
                        rest = pdiv_exact(rest, g, p)
                    zc = pmod(pmul(zc, zeta1, p), Hp, p)
                assert deg(rest) == 0, "digit partition incomplete"
            pieces = newp
        for (Hp, ell) in pieces:
            assert ell % q == 0, "log not divisible by q: not a q-th power"
            rq = ppowmod(pmod(u, Hp, p), ell // q, Hp, p)
            b = pmod(pmul(rq, pmod(root_t, Hp, p), p), Hp, p)
            pieces_all.append((b, Hp))
    # CRT
    r, m = pieces_all[0]
    for (b, Hp) in pieces_all[1:]:
        r, m = crt2(r, m, b, Hp, p)
    return r

