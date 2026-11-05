#!/usr/bin/env python3
"""Exact verifier for the upper bound  C_41b <= T  (ambidextrous moving sofa constant, constant 41b of
teorth/optimizationproblems), T = 19/10 in sofa_certificate.json.  Standard library only; all checks are in exact
rational/integer arithmetic.  The mathematics is in PROOF.md; in short:

  H = R x [0,1], L0 = (-inf,1]x[0,1] u [0,1]x(-inf,1] (unit-width hallway, Kallus-Romik notation),
  C(psi,u) = Rot_psi(L0 + u).  For a sofa S (connected, placed in H at the start) every rigid motion through the
  Z-shaped hallway (and also the "turn right or turn left from the same placement" version) forces, unless
  area(S) <= T, the containments  S c C(psi_j, u_j)  for all hallway angles psi_j in at least two full "quadrant
  windows" [j*pi/2 + eps, (j+1)*pi/2 - eps], sin(eps) = 1/T  (orientation lift + intermediate value theorem +
  Gerver's strip bounds 1/|cos|, 1/|sin|).  The window pairs that can occur are
      Z-hallway:        {0,2} {0,3} {1,2} {1,3}
      either direction: {0,2} {0,3} {2,3} {1,3}  and two angles in quadrant 3 alone.
  For each such case this file certifies  sup_u area(K cap INTERSECT_j C(psi_j,u_j)) <= T,  where K = [0,W]x[0,1]
  (S connected => its x-extent is <= (1+|cos|)/|sin| + (1+|sin|)/|cos| <= W after a horizontal shift) and every u_j
  may be clamped to [min_K p - 1, max_K p] x [min_K q - 1, max_K q] (p, q = rotated coordinates) without changing
  K cap C(psi_j,u_j) up to a null set.  The sup over this compact box is bounded by branch and bound: on a box B
  of translations, C(psi,u) c Pt(B) u Qt(B) (expanded arms), and area(K cap INTERSECT_j (Pt_j u Qt_j)) is computed
  exactly by inclusion-exclusion over {Pt, Qt, Pt&Qt}^k with exactly clipped convex polygons.
Usage:  python3 verify.py [sofa_certificate.json]      (about 4 minutes; exit status 0 iff certified)
"""
import sys, json, base64, hashlib, time, random
from fractions import Fraction as F

# ---------------- exact convex clipping: integer lines, homogeneous integer vertices ----------------
def meet(l1, l2):
    """Intersection of the boundary lines of a1 x + b1 y <= c1 and a2 x + b2 y <= c2 as (X, Y, Z), Z > 0."""
    a1, b1, c1 = l1; a2, b2, c2 = l2
    Z = a1 * b2 - a2 * b1
    if Z == 0: raise ValueError('parallel consecutive edges')
    X = c1 * b2 - c2 * b1; Y = a1 * c2 - a2 * c1
    return (X, Y, Z) if Z > 0 else (-X, -Y, -Z)

def clip(edges, verts, l):
    """Convex polygon given by its cyclic edge lines; verts[i] = edges[i] cap edges[i+1].
    Returns the polygon cut by the half-plane l, or None if the result has zero area."""
    a, b, c = l
    s = [a * X + b * Y - c * Z for (X, Y, Z) in verts]
    if all(v <= 0 for v in s): return edges, verts
    if not any(v < 0 for v in s): return None
    n = len(verts)
    p = next(i for i in range(n) if s[i] > 0 and s[i - 1] <= 0)
    q = p
    while s[(q + 1) % n] > 0: q = (q + 1) % n
    new_edges = []; i = (q + 1) % n
    while True:
        new_edges.append(edges[i])
        if i == p: break
        i = (i + 1) % n
    new_edges.append(l)
    m = len(new_edges)
    return new_edges, [meet(new_edges[k], new_edges[(k + 1) % m]) for k in range(m)]

def area(verts):
    tot = F(0); n = len(verts)
    for i in range(n):
        X1, Y1, Z1 = verts[i]; X2, Y2, Z2 = verts[(i + 1) % n]
        tot += F(X1 * Y2 - X2 * Y1, Z1 * Z2)
    return abs(tot) / 2

def K_polygon(W):
    edges = [(0, -1, 0), (W.denominator, 0, W.numerator), (0, 1, 1), (-1, 0, 0)]   # y>=0, x<=W, y<=1, x>=0
    return edges, [meet(edges[i], edges[(i + 1) % 4]) for i in range(4)]

# ---------------- corridors ----------------
def quadrant(j, c, s):
    """(cos, sin) of j*pi/2 + phi where cos phi = c, sin phi = s."""
    return [(c, s), (-s, c), (-c, -s), (s, -c)][j]

def arm_lines(cs, box, scale):
    """Expanded arms of C(psi,u), u in box=(a1,b1,a2,b2): Pt = {p<=b1+1, a2<=q<=b2+1}, Qt = {a1<=p<=b1+1, q<=b2},
    PQ = Pt cap Qt, with p = cos x + sin y, q = -sin x + cos y.  Integer lines after scaling by den*scale."""
    (c, s) = cs
    den = c.denominator * s.denominator
    cn, sn = int(c * den), int(s * den)
    def L(kp, kq, r):    # kp*p + kq*q <= r  ->  integer line
        a = (kp * cn - kq * sn) * scale; b = (kp * sn + kq * cn) * scale; rr = r * den * scale
        assert rr.denominator == 1
        return (a, b, int(rr))
    a1, b1, a2, b2 = box
    Pt = [L(1, 0, b1 + 1), L(0, -1, -a2), L(0, 1, b2 + 1)]
    Qt = [L(-1, 0, -a1), L(1, 0, b1 + 1), L(0, 1, b2)]
    PQ = [L(-1, 0, -a1), L(1, 0, b1 + 1), L(0, -1, -a2), L(0, 1, b2)]
    return Pt, Qt, PQ

def region_area(rots, boxes, K):
    """Exact area of K cap INTERSECT_j (Pt_j u Qt_j) by inclusion-exclusion."""
    scale = 1
    for bx in boxes:
        for v in bx: scale = max(scale, v.denominator)
    pcs = [arm_lines(cs, bx, scale) for cs, bx in zip(rots, boxes)]
    total = F(0)
    stack = [(0, K[0], K[1], 1)]
    while stack:
        j, ed, ve, sg = stack.pop()
        if j == len(pcs):
            total += sg * area(ve); continue
        for lines, s2 in ((pcs[j][0], 1), (pcs[j][1], 1), (pcs[j][2], -1)):
            e, v, ok = ed, ve, True
            for l in lines:
                r = clip(e, v, l)
                if r is None: ok = False; break
                e, v = r
            if ok: stack.append((j + 1, e, v, sg * s2))
    return total

# ---------------- self-test: exact area vs. grid sampling of the same region ----------------
def self_test(W):
    rng = random.Random(7)
    K = K_polygon(W)
    worst = 0.0
    for trial in range(6):
        js = [rng.randrange(4), rng.randrange(4)]
        rots = [quadrant(j, F(21, 29), F(20, 29)) for j in js]
        boxes = []
        for _ in rots:
            a1 = F(rng.randint(-64, 160), 64); a2 = F(rng.randint(-200, 30), 64)
            boxes.append((a1, a1 + F(rng.randint(1, 64), 64), a2, a2 + F(rng.randint(1, 64), 64)))
        exact = float(region_area(rots, boxes, K))
        nx, ny, cnt = 600, 150, 0
        for ix in range(nx):
            x = (ix + 0.5) * float(W) / nx
            for iy in range(ny):
                y = (iy + 0.5) / ny; inside = True
                for (c, s), (a1, b1, a2, b2) in zip(rots, boxes):
                    p = float(c) * x + float(s) * y; q = -float(s) * x + float(c) * y
                    inP = p <= b1 + 1 and a2 <= q <= b2 + 1
                    inQ = a1 <= p <= b1 + 1 and q <= b2
                    if not (inP or inQ): inside = False; break
                if inside: cnt += 1
        grid = cnt * float(W) / (nx * ny)
        worst = max(worst, abs(exact - grid))
    assert worst < 0.02, 'self-test failed: exact area disagrees with grid sampling (%g)' % worst
    sq = [(0, -1, 0), (1, 0, 1), (0, 1, 1), (-1, 0, 0)]
    e, v = clip(sq, [meet(sq[i], sq[(i + 1) % 4]) for i in range(4)], (1, 1, 1))   # unit square cap x+y<=1
    assert area(v) == F(1, 2)
    print('self-test ok (exact clipping vs grid sampling, max deviation %.4f)' % worst)

# ---------------- branch-and-bound tree check ----------------
def check_case(case, T, W, verbose=True):
    corr = [(int(j), F(c), F(s)) for j, c, s in case['corridors']]
    for j, c, s in corr:
        assert j in (0, 1, 2, 3) and c > 0 and s > 0 and c * c + s * s == 1, 'bad corridor angle'
        assert 1 / T <= s and 1 / T <= c, 'corridor base angle outside the window [eps, pi/2-eps]'
    # x-extent bound for a connected subset of H cap (one hallway of this case)
    assert W >= min((1 + c) / s + (1 + s) / c for j, c, s in corr), 'W too small'
    rots = [quadrant(j, c, s) for j, c, s in corr]
    dom = [[F(v) for v in d] for d in case['domain']]
    Kv = [(F(0), F(0)), (W, F(0)), (W, F(1)), (F(0), F(1))]
    for (c, s), (a1, b1, a2, b2) in zip(rots, dom):
        ps = [c * x + s * y for x, y in Kv]; qs = [-s * x + c * y for x, y in Kv]
        assert a1 <= min(ps) - 1 and max(ps) <= b1 and a2 <= min(qs) - 1 and max(qs) <= b2, 'domain too small'
        assert a1 < b1 and a2 < b2
    raw = base64.b64decode(case['tree_b64'])
    nbits = int.from_bytes(raw[:8], 'little'); bits = raw[8:]
    assert len(bits) == (nbits + 7) // 8, 'tree length mismatch'
    K = K_polygon(W)
    boxes = dom; pos = 0; leaves = 0; worst = F(0); t0 = time.time()
    stack = ['visit']
    while stack:
        op = stack.pop()
        if op != 'visit':
            i, k, val = op; boxes[i][k] = val; continue
        assert pos < nbits, 'tree truncated'
        bit = (bits[pos >> 3] >> (pos & 7)) & 1; pos += 1
        if bit == 0:
            ub = region_area(rots, [tuple(b) for b in boxes], K)
            assert ub <= T, 'leaf bound %.12f > T at %s' % (float(ub), boxes)
            worst = max(worst, ub); leaves += 1
            if verbose and leaves % 100000 == 0:
                print('   ... %d leaves, worst %.12f, %.0f s' % (leaves, float(worst), time.time() - t0), flush=True)
            continue
        best = None                       # widest coordinate, first index on ties, split at the midpoint
        for i, (a1, b1, a2, b2) in enumerate(boxes):
            for k, w in ((0, b1 - a1), (2, b2 - a2)):
                if best is None or w > best[0]: best = (w, i, k)
        _, i, k = best
        lo, hi = boxes[i][k], boxes[i][k + 1]; mid = (lo + hi) / 2
        stack += [(i, k, lo), 'visit', (i, k, mid), (i, k + 1, hi), 'visit']
        boxes[i][k + 1] = mid             # lower half first
    assert pos == nbits, 'tree not fully used'
    return leaves, worst, sorted(j for j, c, s in corr)

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'sofa_certificate.json'
    raw = open(path, 'rb').read()
    print('certificate', path, 'sha256', hashlib.sha256(raw).hexdigest())
    d = json.loads(raw)
    T = F(d['T']); W = F(d['W'])
    assert T * T > 2, 'need T > sqrt(2) (eps < pi/4)'
    self_test(W)
    done = []
    t0 = time.time()
    for case in d['cases']:
        leaves, worst, quads = check_case(case, T, W)
        done.append(quads)
        print('%-16s quadrants %s: %7d leaves, max leaf bound %.12f <= %s   (%.0f s)' % (
            case['name'], quads, leaves, float(worst), T, time.time() - t0), flush=True)
    have = set(tuple(q) for q in done)
    z_need = [(0, 2), (0, 3), (1, 2), (1, 3)]
    e_need = [(0, 2), (0, 3), (2, 3), (1, 3), (3, 3)]
    assert all(q in have for q in z_need), 'missing a case of the Z-hallway analysis'
    assert all(q in have for q in e_need), 'missing a case of the either-direction analysis'
    print('CERTIFIED: every ambidextrous sofa (Z-shaped hallway, or right and left turn from one placement) has area <= %s = %s'
          % (T, float(T)))

if __name__ == '__main__':
    main()
