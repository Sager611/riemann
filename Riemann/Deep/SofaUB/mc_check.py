"""Independent float/rational helpers for verify.py: admissibility, true tube volume, root coverage."""
import math, subprocess
from fractions import Fraction as F
from exact_check import Dfam


def need(cond, msg):
    """explicit check (not stripped by python -O)"""
    if not cond:
        raise RuntimeError(msg)


def tube_geom(p):
    kx, kz, t, X1, X2 = p
    s = math.sqrt(1 + kx*kx + kz*kz)
    f = 1 / (1 + s)
    A = ((1 + f*kx*kx, f*kx*kz), (f*kx*kz, 1 + f*kz*kz))
    ct, st = (1 - t*t) / (1 + t*t), 2*t / (1 + t*t)
    g1 = (A[0][0]*ct + A[0][1]*st, A[1][0]*ct + A[1][1]*st)
    g2 = (-A[0][0]*st + A[0][1]*ct, -A[1][0]*st + A[1][1]*ct)
    return (kx, kz), (X1, X2), g1, g2


def slice_poly(p, y, ya):
    k, X, g1, g2 = tube_geom(p)
    cx, cz = X[0] + (y - ya)*k[0], X[1] + (y - ya)*k[1]
    return [(cx + a*g1[0]/2 + b*g2[0]/2, cz + a*g1[1]/2 + b*g2[1]/2) for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]


def sat_gap(P, Q):
    """max separation over edge normals (> 0: disjoint, < 0: overlap depth)"""
    best = -1e300
    for poly in (P, Q):
        n = len(poly)
        for i in range(n):
            ex, ez = poly[(i+1) % n][0] - poly[i][0], poly[(i+1) % n][1] - poly[i][1]
            nx, nz = -ez, ex
            L = math.hypot(nx, nz)
            if L == 0: continue
            nx, nz = nx/L, nz/L
            pa = [q[0]*nx + q[1]*nz for q in P]; pb = [q[0]*nx + q[1]*nz for q in Q]
            best = max(best, min(pb) - max(pa), min(pa) - max(pb))
    return best


SQ = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]


def admissible(p, u, ya, margin=1e-9):
    """True if the tube is admissible with margin (tilt and both windows)"""
    kx, kz = p[0], p[1]
    K2 = 8*u*u / (u*u - 2)**2
    if kx*kx + kz*kz > K2*(1 - 1e-12): return False
    c = 1 + u
    for y in (c, 1.0):
        if sat_gap(slice_poly(p, y, ya), SQ) > -margin: return False
    return True


def tube_volume(p, u, ya):
    k, X, g1, g2 = tube_geom(p)
    fam = Dfam(k, X, g1, g2, 0, 0, 0, 0, ya)
    return float(fam.integral(1 + u))


def check_roots(u, ya, grid, meta):
    exe = meta["_exe"]
    out = subprocess.run([exe, "roots", repr(u), repr(ya)] + [str(g) for g in grid], capture_output=True, text=True, check=True).stdout
    lines = out.split("\n")
    U = F(u); C = 1 + U; YA = F(ya)
    K2 = 8*U*U / (U*U - 2)**2
    KHI = F(float(lines[0].split()[1]))   # exact double
    need(KHI*KHI >= K2, "KHI below K")
    cells = [l.split() for l in lines if l.startswith("CELL")]
    nkx, nkz, nt, nx = grid
    need(len(cells) == nkx*nkz, "check failed: len(cells) == nkx*nkz")
    # k-cells tile [-KHI,KHI] x [0,KHI]
    if any(c[3] == "EMPTY" for c in cells):
        raise AssertionError("empty cell needs a proof")  # not expected for our grids
    q = lambda v: F(float(v))          # exact value of the printed double (%.17g round-trips)
    cell = {(int(c[1]), int(c[2])): c for c in cells}
    for jx in range(nkx):
        for jz in range(nkz):
            c = cell[(jx, jz)]
            need(q(c[3]) == (-KHI if jx == 0 else q(cell[(jx - 1, jz)][4])), "kx cells do not tile")
            need(jx < nkx - 1 or q(c[4]) == KHI, "kx cells do not reach KHI")
            need(q(c[5]) == (0 if jz == 0 else q(cell[(jx, jz - 1)][6])), "kz cells do not tile")
            need(jz < nkz - 1 or q(c[6]) == KHI, "kz cells do not reach KHI")
    for c in cells:
        klo, khi, zlo, zhi = map(q, c[3:7])
        Xlo, Xhi, Zlo, Zhi = map(q, c[8:12])
        ax = max(abs(klo), abs(khi)); az = max(abs(zlo), abs(zhi))
        _, s_hi = frac_sqrt(1 + ax*ax + az*az)
        _, r_hi = frac_sqrt(F(1, 2))
        e = s_hi * r_hi
        for lo_c, hi_c, k_lo, k_hi in ((Xlo, Xhi, klo, khi), (Zlo, Zhi, zlo, zhi)):
            # need [lo_c, hi_c] to contain the exact union bounds of both window constraints' intersection
            l1 = -e + min((YA - C)*k_lo, (YA - C)*k_hi); h1 = 1 + e + max((YA - C)*k_lo, (YA - C)*k_hi)
            l2 = -e + min((YA - 1)*k_lo, (YA - 1)*k_hi); h2 = 1 + e + max((YA - 1)*k_lo, (YA - 1)*k_hi)
            need(lo_c <= max(l1, l2) and hi_c >= min(h1, h2), "X range too small")
    tc = [l.split() for l in lines if l.startswith("TCELL")]
    need(q(tc[0][2]) == 0 and q(tc[-1][3]) == 1, "check failed: q(tc[0][2]) == 0 and q(tc[-1][3]) == 1")
    for a, b in zip(tc[:-1], tc[1:]):
        need(q(a[3]) == q(b[2]), "t cells do not tile")
    return True


def frac_sqrt(x, bits=96):
    x = F(x)
    s = math.isqrt(x.numerator * 4**bits // x.denominator)
    lo = F(s, 2**bits)
    while lo*lo > x: lo -= F(1, 2**bits)
    hi = lo + F(1, 2**bits)
    while hi*hi < x: hi += F(1, 2**bits)
    return lo, hi
