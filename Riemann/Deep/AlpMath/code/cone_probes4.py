# Probes N1 (HasseNormSt, F = Q), K1 (ChebotarevDegreeOne instances), J1 (BrauerDegreeZero on small groups)
import sys, math, itertools, random
import numpy as np
from sympy import factorint, primerange, isprime
random.seed(13)
def rep(tag, ok, msg): print(f"{tag}: {'PASS' if ok else 'FAIL'} {msg}"); sys.stdout.flush()
# ---------- N1 Hasse norm theorem over Q for L = Q(sqrt delta)
def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
def val(a, p):
    v = 0
    while a % p == 0: a //= p; v += 1
    return v, a
def hilb(a, b, p):
    if p == 0: return -1 if (a < 0 and b < 0) else 1
    va, ua = val(a, p); vb, ub = val(b, p)
    if p != 2:
        return (-1) ** (va * vb * ((p - 1) // 2)) * leg(ua, p) ** vb * leg(ub, p) ** va
    e = lambda u: ((u - 1) // 2) % 2; w = lambda u: ((u * u - 1) // 8) % 2
    return (-1) ** (e(ua) * e(ub) + va * w(ub) + vb * w(ua))
issq = lambda n: n >= 0 and math.isqrt(n) ** 2 == n
cases = 0; bad = 0; locfail = 0
for delta in range(-30, 31):
    if delta == 0 or issq(delta): continue
    for x in range(-30, 31):
        if x == 0: continue
        ps = set(factorint(abs(x * delta))) | {2}
        if any(hilb(x, delta, p) != 1 for p in list(ps) + [0]): locfail += 1; continue
        cases += 1; found = False
        for Z in range(1, 121):
            for Y in range(0, 121):
                if issq(delta * Y * Y + x * Z * Z): found = True; break
            if found: break
        if not found: bad += 1
rep("N1 HasseNormSt (F = Q): local norms everywhere => global norm from Q(sqrt delta)", bad == 0, f"{cases} pairs (delta, x) in [-30,30]^2 locally norms everywhere: rational solution of X^2 - delta Y^2 = x Z^2 found for all (|Y|,|Z| <= 120); {locfail} pairs fail locally")
# ---------- K1 Chebotarev degree-one instances
ok = True
for m in range(3, 31):
    for a in range(1, m):
        if math.gcd(a, m) != 1: continue
        if not any(p % m == a for p in primerange(500, 20000)): ok = False
types = {}
for p in primerange(500, 5000):
    if p in (2, 3): continue
    r = sum(1 for t in range(p) if (t ** 3 - 2) % p == 0)
    types.setdefault({3: 'id', 1: 'transp', 0: '3-cycle'}[r], p)
ok &= set(types) == {'id', 'transp', '3-cycle'}
rep("K1 ChebotarevDegreeOne instances", ok, f"Q(zeta_m), m <= 30: every class a mod m has a prime p > 500; splitting field of x^3-2 (S3): first p > 500 per Frobenius class {types}")
# ---------- J1 Brauer induction in degree zero: chi - chi(1) in Z-span of Ind_H(lam) - Ind_H(1)
def mul(p, q): return tuple(p[i] for i in q)  # (p*q)(i) = p(q(i))
def inv(p):
    r = [0] * len(p)
    for i, j in enumerate(p): r[j] = i
    return tuple(r)
def gen_group(gens):
    e = tuple(range(len(gens[0]))); G = {e}; fr = [e]
    while fr:
        x = fr.pop()
        for g in gens:
            y = mul(g, x)
            if y not in G: G.add(y); fr.append(y)
    return sorted(G)
def classes(G):
    cl = {}; reps = []
    for g in G:
        if g in cl: continue
        orb = {mul(mul(h, g), inv(h)) for h in G}
        for o in orb: cl[o] = len(reps)
        reps.append(g)
    return cl, reps
def char_table(G):
    cl, reps = classes(G); k = len(reps); order = len(G)
    sizes = np.zeros(k)
    for g in G: sizes[cl[g]] += 1
    c = np.zeros((k, k, k))
    for kk, z in enumerate(reps):
        for x in G:
            c[cl[x], cl[mul(inv(x), z)], kk] += 1
    A = np.tensordot(np.random.randn(k), c, axes=(0, 0))
    ev, V = np.linalg.eig(A)
    idc = cl[tuple(range(len(G[0])))]
    chars = []
    for col in range(k):
        w = V[:, col] / V[idc, col]
        deg = np.sqrt((order / np.sum(np.abs(w) ** 2 / sizes)).real)
        chars.append(w * deg / sizes)
    return cl, reps, sizes, np.array(chars), idc
def hnf_solve(gens, t):
    # integer lattice membership: t in Z-span(gens) ? (small dims) via row-reduction (Hermite) of the generator matrix
    M = [list(v) for v in gens if any(v)]
    n = len(t); rows = []
    col = 0
    while M and col < n:
        M = [r for r in M if any(r)]
        piv = [r for r in M if r[col] != 0]
        if not piv: col += 1; continue
        while len([r for r in M if r[col] != 0]) > 1:
            nz = sorted([r for r in M if r[col] != 0], key=lambda r: abs(r[col]))
            p0 = nz[0]
            for r in nz[1:]:
                q = r[col] // p0[col]
                for j in range(n): r[j] -= q * p0[j]
        p0 = next(r for r in M if r[col] != 0)
        rows.append((col, p0[:])); M.remove(p0); col += 1
    t = list(t)
    for c0, r in rows:
        if t[c0] % r[c0]: return False
        q = t[c0] // r[c0]
        for j in range(n): t[j] -= q * r[j]
    return all(v == 0 for v in t)
def group_by_name(name):
    if name == 'S3': return gen_group([(1, 0, 2), (1, 2, 0)])
    if name == 'D4': return gen_group([(1, 2, 3, 0), (0, 3, 2, 1)])
    if name == 'A4': return gen_group([(1, 2, 0, 3), (1, 0, 3, 2)])
    if name == 'S4': return gen_group([(1, 2, 3, 0), (1, 0, 2, 3)])
    if name == 'Q8':  # regular representation of Q8 = {+-1, +-i, +-j, +-k}
        els = [(s, u) for u in '1ijk' for s in (1, -1)]
        tab = {('1', '1'): (1, '1'), ('1', 'i'): (1, 'i'), ('1', 'j'): (1, 'j'), ('1', 'k'): (1, 'k'), ('i', '1'): (1, 'i'), ('i', 'i'): (-1, '1'), ('i', 'j'): (1, 'k'), ('i', 'k'): (-1, 'j'),
               ('j', '1'): (1, 'j'), ('j', 'i'): (-1, 'k'), ('j', 'j'): (-1, '1'), ('j', 'k'): (1, 'i'), ('k', '1'): (1, 'k'), ('k', 'i'): (1, 'j'), ('k', 'j'): (-1, 'i'), ('k', 'k'): (-1, '1')}
        def m(a, b):
            s, u = tab[(a[1], b[1])]; return (a[0] * b[0] * s, u)
        idx = {e: i for i, e in enumerate(els)}
        perms = [tuple(idx[m(g, h)] for h in els) for g in [(1, 'i'), (1, 'j')]]
        return gen_group(perms)
allok = True; summary = []
for name in ['S3', 'D4', 'Q8', 'A4', 'S4']:
    G = group_by_name(name); Gs = set(G)
    cl, reps, sizes, chars, idc = char_table(G); k = len(reps); order = len(G)
    subs = set()
    for a in G:
        for b in G:
            subs.add(tuple(gen_group([a, b])))
    gensv = []
    for H in subs:
        H = list(H); Hs = set(H)
        clH, repsH, sizesH, charsH, idcH = char_table(H)
        lin = [ch for ch in charsH if abs(ch[idcH] - 1) < 1e-9]
        def ind(fH):
            vals = np.zeros(k, dtype=complex)
            for ci, g in enumerate(reps):
                s = 0
                for x in G:
                    y = mul(mul(x, g), inv(x))
                    if y in Hs: s += fH[clH[y]]
                vals[ci] = s / len(H)
            return vals
        one = ind(np.ones(len(repsH)))
        for lam in lin:
            v = ind(lam) - one
            coords = [(np.sum(sizes * v * np.conj(ch)) / order) for ch in chars]
            ci = [int(round(z.real)) for z in coords]
            assert max(abs(z - c) for z, c in zip(coords, ci)) < 1e-6
            gensv.append(ci)
    triv = next(i for i in range(k) if all(abs(chars[i] - 1) < 1e-9))
    okG = True
    for i in range(k):
        t = [0] * k; t[i] += 1; t[triv] -= int(round(chars[i][idc].real))
        if not hnf_solve([g[:] for g in gensv], t): okG = False
    allok &= okG; summary.append(f"{name}({order}):{'ok' if okG else 'FAIL'}")
rep("J1 BrauerDegreeZero: chi - chi(1) in Z-span{Ind_H lam - Ind_H 1}", allok, ", ".join(summary) + " (all 2-generated subgroups, all linear characters)")
