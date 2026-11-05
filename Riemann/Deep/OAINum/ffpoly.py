
# ffpoly.py -- polynomial arithmetic over F_p (coefficients low -> high), linear algebra mod p
import random
import numpy as np

class PF:
    def __init__(self, p):
        self.p = p
    def inv(self, a):
        return pow(a % self.p, self.p - 2, self.p)

def norm(f):
    while f and f[-1] == 0:
        f.pop()
    return f

def padd(f, g, p):
    n = max(len(f), len(g))
    r = [((f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)) % p for i in range(n)]
    return norm(r)

def psub(f, g, p):
    n = max(len(f), len(g))
    r = [((f[i] if i < len(f) else 0) - (g[i] if i < len(g) else 0)) % p for i in range(n)]
    return norm(r)

def pneg(f, p):
    return [(-c) % p for c in f]

def pscal(f, c, p):
    c %= p
    if c == 0:
        return []
    return norm([(x * c) % p for x in f])

def pmul(f, g, p):
    if not f or not g:
        return []
    if len(f) > 40 and len(g) > 40:
        # use numpy convolution in chunks to avoid overflow: split into 15-bit limbs
        return pmul_np(f, g, p)
    r = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a == 0:
            continue
        for j, b in enumerate(g):
            r[i + j] += a * b
    return norm([x % p for x in r])

def pmul_np(f, g, p):
    # exact via object-free method: split coefficients into 16-bit halves
    fa = np.array(f, dtype=np.int64); ga = np.array(g, dtype=np.int64)
    B = 1 << 16
    f0, f1 = fa % B, fa // B
    g0, g1 = ga % B, ga // B
    # each convolution entry < len * 2^32 < 2^63 for len < 2^31
    c00 = np.convolve(f0, g0) % p
    c01 = (np.convolve(f0, g1) % p + np.convolve(f1, g0) % p) % p
    c11 = np.convolve(f1, g1) % p
    Bp = B % p
    r = (c00 + (c01 * Bp) % p + (c11 * ((Bp * Bp) % p)) % p) % p
    return norm([int(x) for x in r])

def pdivmod(f, g, p):
    if not g:
        raise ZeroDivisionError
    f = list(f)
    dg = len(g) - 1
    inv = pow(g[-1], p - 2, p)
    if len(f) <= dg:
        return [], norm(f)
    q = [0] * (len(f) - dg)
    for i in range(len(f) - 1, dg - 1, -1):
        c = (f[i] * inv) % p
        if c:
            q[i - dg] = c
            for j in range(dg + 1):
                f[i - dg + j] = (f[i - dg + j] - c * g[j]) % p
    return norm(q), norm(f[:dg])

def pmod(f, g, p):
    return pdivmod(f, g, p)[1]

def pdiv_exact(f, g, p):
    q, r = pdivmod(f, g, p)
    assert not r, "non-exact division"
    return q

def monic(f, p):
    if not f:
        return []
    inv = pow(f[-1], p - 2, p)
    return [(c * inv) % p for c in f]

def pgcd(f, g, p):
    f, g = norm(list(f)), norm(list(g))
    while g:
        f, g = g, pmod(f, g, p)
    return monic(f, p)

def pxgcd(f, g, p):
    # returns (d, s, t) with s f + t g = d monic
    r0, r1 = norm(list(f)), norm(list(g))
    s0, s1 = [1], []
    t0, t1 = [], [1]
    while r1:
        q, r = pdivmod(r0, r1, p)
        r0, r1 = r1, r
        s0, s1 = s1, psub(s0, pmul(q, s1, p), p)
        t0, t1 = t1, psub(t0, pmul(q, t1, p), p)
    if not r0:
        return [], [], []
    inv = pow(r0[-1], p - 2, p)
    return pscal(r0, inv, p), pscal(s0, inv, p), pscal(t0, inv, p)

def pinvmod(a, m, p):
    d, s, t = pxgcd(a, m, p)
    if d != [1]:
        raise ValueError("not invertible")
    return pmod(s, m, p)

def pderiv(f, p):
    return norm([(i * f[i]) % p for i in range(1, len(f))])

def ppowmod(a, e, m, p):
    r = [1]
    a = pmod(a, m, p)
    while e > 0:
        if e & 1:
            r = pmod(pmul(r, a, p), m, p)
        e >>= 1
        if e:
            a = pmod(pmul(a, a, p), m, p)
    return r if m != [1] else []

def ppow(a, e, p):
    r = [1]
    while e > 0:
        if e & 1:
            r = pmul(r, a, p)
        e >>= 1
        if e:
            a = pmul(a, a, p)
    return r

def peval(f, x, p):
    r = 0
    for c in reversed(f):
        r = (r * x + c) % p
    return r

def deg(f):
    return len(f) - 1

def squarefree_decomp(f, p):
    # Yun's algorithm, char p > deg f assumed; returns list [(g_i, i)] with f = lc * prod g_i^i
    f = monic(f, p)
    out = []
    if deg(f) <= 0:
        return out
    fp = pderiv(f, p)
    a = pgcd(f, fp, p)
    b = pdiv_exact(f, a, p)
    c = pdiv_exact(fp, a, p) if a != [1] else fp
    # monic adjust: Yun with monic polys
    d = psub(c, pderiv(b, p), p)
    i = 1
    while deg(b) > 0:
        a = pgcd(b, d, p)
        if deg(a) > 0:
            out.append((a, i))
        b = pdiv_exact(b, a, p)
        c = pdiv_exact(d, a, p)
        d = psub(c, pderiv(b, p), p)
        i += 1
    return out

# ---------- linear algebra mod p (numpy int64, p < 2^31) ----------
def rref(Ain, p, inplace=False):
    # memory-lean Gauss-Jordan: in-place, chunked row updates (temporaries <= CH x cols)
    if inplace and isinstance(Ain, np.ndarray) and Ain.dtype == np.int64:
        A = Ain
        A %= p
    else:
        A = np.array(Ain, dtype=np.int64) % p
    rows, cols = A.shape
    piv = []
    r = 0
    CH = 256
    for c in range(cols):
        if r >= rows:
            break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        k = r + nz[0]
        if k != r:
            tmp = A[r].copy(); A[r] = A[k]; A[k] = tmp
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        pr = A[r]
        nzr = np.nonzero(A[:, c])[0]
        nzr = nzr[nzr != r]
        for s0 in range(0, len(nzr), CH):
            idx = nzr[s0:s0 + CH]
            blk = A[idx]
            blk -= (blk[:, c:c+1] * pr[None, :]) % p
            blk %= p
            A[idx] = blk
        piv.append(c)
        r += 1
    return A[:r], piv

def nullspace(Ain, ncols, p):
    if len(Ain) == 0:
        return [np.eye(ncols, dtype=np.int64)[i] for i in range(ncols)]
    R, piv = rref(Ain, p, inplace=True)
    free = [c for c in range(ncols) if c not in set(piv)]
    basis = []
    pivset = {c: i for i, c in enumerate(piv)}
    for f in free:
        v = np.zeros(ncols, dtype=np.int64)
        v[f] = 1
        for c, i in pivset.items():
            v[c] = (-R[i, f]) % p
        basis.append(v)
    return basis

def solve_lin(Ain, bvec, p):
    # find one solution x of A x = b (mod p) or None
    A = np.array(Ain, dtype=np.int64) % p
    b = np.array(bvec, dtype=np.int64).reshape(-1, 1) % p
    Ab = np.concatenate([A, b], axis=1)
    R, piv = rref(Ab, p)
    n = A.shape[1]
    if n in piv:
        return None
    x = np.zeros(n, dtype=np.int64)
    for i, c in enumerate(piv):
        x[c] = R[i, n]
    return x

