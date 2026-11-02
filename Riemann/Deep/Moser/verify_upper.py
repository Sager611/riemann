#!/usr/bin/env python3
"""Independent exact checker for upper-bound certificates (support/length method).

Re-implements, from the statements in arXiv:2609.21968 Sec. 4 (Lemmas 4.2-4.4 and
the leaf conditions (4.7)), every finite check needed by the proof of Theorem 4.1,
for a convex polygon K = {x : n_i . x <= h_i}, i = 0..k-1 (k >= 4), n_i = n(sigma_i),
with h_0 = h_1 = 0 and every n_j (j >= 2) in the open negative cone of n_0, n_1.
For k = 4 this is exactly the paper's setting; for k > 4 the translation criterion
(Lemma 4.4) holds verbatim (one escape alternative per side j >= 2), and an
orientation node of the case tree has k - 2 children, child r for side 2 + r.

Leaf conditions are checked in homogeneous form: with B = sum of mu_k h_j over the
selected escape rows, require B > 0, sum_lambda f_lambda = 0 and |sum_{I} f|^2 <= B^2
for every I in F(R); dividing all multipliers by B gives the paper's (4.7) exactly.
Accepts the paper's file format (multipliers-*.jsonl.gz, records [i, den, [[key, num],..]])
and the compact format (multipliers-*.jsonl.xz, records [i, den, key_deltas, nums]).

Python standard library only; exact rational arithmetic (fractions.Fraction).
Usage: verify_upper.py CERT_DIR      (exit status 0 iff accepted)
"""
import gzip
import hashlib
import json
import lzma
import sys
from fractions import Fraction as Q
from pathlib import Path
from typing import NoReturn


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


def unit(t):
    return ((1 - t * t) / (1 + t * t), 2 * t / (1 + t * t))


def fail(msg) -> NoReturn:
    raise SystemExit("REJECT: " + msg)


def angle_key_less(u, v):
    """Strict counterclockwise angle comparison with angles in [0, 2pi)."""
    def sec(w):
        if w[1] == 0:
            return 0 if w[0] > 0 else 2
        return 1 if w[1] > 0 else 3
    a, b = sec(u), sec(v)
    return a < b if a != b else cross(u, v) > 0


def check_shape(T, H):
    k = len(T)
    if k < 4 or len(H) != k:
        fail("need k >= 4 sides")
    N = [unit(t) for t in T]
    for n in N:
        if n[0] * n[0] + n[1] * n[1] != 1:
            fail("normal not unit")
    if H[0] != 0 or H[1] != 0:
        fail("h0 = h1 = 0 required")
    if any(h <= 0 for h in H[2:]):
        fail("h_j > 0 required for j >= 2")
    # consecutive normals strictly counterclockwise, each turn in (0, pi)
    for i in range(k):
        if cross(N[i], N[(i + 1) % k]) <= 0:
            fail("normals not in strict convex counterclockwise order")
    # total turning exactly 2 pi, checked exactly (no floating point): every turn
    # N[i] -> N[i+1] lies in (0, pi) by the check above, so with angles taken in
    # [0, 2 pi) a step descends (angle_key_less fails) iff it wraps past angle 0,
    # and the total turning is 2 pi * (number of cyclic descents). Require one.
    descents = sum(1 for i in range(k) if not angle_key_less(N[i], N[(i + 1) % k]))
    if descents != 1:
        fail("normals wind more than once")
    det = cross(N[0], N[1])
    rays = []
    for j in range(2, k):
        v = (-N[j][0], -N[j][1])
        a, b = cross(v, N[1]) / det, cross(N[0], v) / det
        if a <= 0 or b <= 0:
            fail("side %d normal not in open negative cone of sides 0,1" % j)
        if (a * N[0][0] + b * N[1][0] + N[j][0], a * N[0][1] + b * N[1][1] + N[j][1]) != (0, 0):
            fail("ray identity")
        rays.append((j, a, b))
    # vertices V_i = line (i-1) cap line i, all sides active, convex polygon
    V = []
    for i in range(k):
        p, q = N[i - 1], N[i]
        hp, hq = H[i - 1], H[i]
        d = cross(p, q)
        V.append(((hp * q[1] - p[1] * hq) / d, (p[0] * hq - hp * q[0]) / d))
    for i, v in enumerate(V):
        for j in range(k):
            s = N[j][0] * v[0] + N[j][1] * v[1]
            if j in ((i - 1) % k, i):
                if s != H[j]:
                    fail("vertex not on its lines")
            elif s >= H[j]:
                fail("side %d not strictly active (vertex %d on/outside)" % (j, i))
    area2 = sum(cross(V[i], V[(i + 1) % k]) for i in range(k))
    if area2 <= 0:
        fail("nonpositive area")
    diam2 = max((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for a in V for b in V)
    if diam2 < 1:
        fail("K contains no unit segment")
    return N, rays, V, area2 / 2


def matrix(spec, N):
    kind = spec[0]
    if kind == "flush":
        side, flip, mirror = spec[1], spec[2], spec[3]
        if side not in range(len(N)) or flip not in (0, 1) or mirror not in (0, 1):
            fail("bad flush spec")
        x, y = N[side]
        m = ((y, -x), (x, y)) if flip else ((-y, x), (-x, -y))
    elif kind in ("rot", "vrot"):
        if len(spec) != 3 or spec[2] not in (0, 1):
            fail("bad rot spec")
        x, y = unit(Q(spec[1]))
        m = ((x, -y), (y, x))
        mirror = spec[2]
    else:
        fail("unknown spec kind")
    if mirror:
        m = ((-m[0][0], -m[0][1]), m[1])
    (a, b), (c, d) = m
    if (a * a + c * c, b * b + d * d, a * b + c * d) != (1, 1, 0):
        fail("orientation not orthogonal")
    return m


def apply(m, v):
    return (m[0][0] * v[0] + m[0][1] * v[1], m[1][0] * v[0] + m[1][1] * v[1])


def main(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, digest in manifest["sha256"].items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != digest:
            fail("sha256 mismatch " + name)
    header = json.loads((directory / "header.json").read_text())
    T = [Q(t) for t in header["t"]]
    H = [Q(h) for h in header["h"]]
    k = len(T)
    N, rays, V, area = check_shape(T, H)
    specs = header["specs"]
    mats = [matrix(s, N) for s in specs]
    # contact directions
    dirs = set()
    for s, m in zip(specs, mats):
        if s[0] != "vrot":
            for n in N:
                dirs.add(apply(m, n))
    for u in dirs:
        if u[0] == 0 and u not in ((0, 1), (0, -1)):
            fail("vertical direction not unit")
    right = sorted((u for u in dirs if u[0] > 0), key=lambda u: u[1])
    left = sorted((u for u in dirs if u[0] < 0), key=lambda u: -u[1])
    for lst, sgn in ((right, 1), (left, -1)):
        for a, b in zip(lst, lst[1:]):
            if not (sgn * (b[1] - a[1]) > 0):
                fail("contact directions not strictly ordered")
    kR, kL = len(right), len(left)
    if (kR, kL) != (header["kR"], header["kL"]):
        fail("kR/kL mismatch")
    normals = [(Q(0), Q(-1)), (Q(0), Q(-1)), (Q(0), Q(1))] + right + left
    n = len(normals)
    label = {u: i for i, u in enumerate(normals)}  # (0,-1) -> 1 (P+), (0,1) -> 2
    named = sorted(label, key=__import__("functools").cmp_to_key(
        lambda u, v: -1 if angle_key_less(u, v) else (1 if angle_key_less(v, u) else 0)))

    def expansion(u):
        if u in label:
            return [(u, Q(1))]
        i = next((i for i, v in enumerate(named) if angle_key_less(u, v)), 0)
        w1, w2 = named[i - 1], named[i]
        d = cross(w1, w2)
        a, b = cross(u, w2) / d, cross(w1, u) / d
        if a < 0 or b < 0 or (a * w1[0] + b * w2[0], a * w1[1] + b * w2[1]) != u:
            fail("invalid nonnegative expansion")
        return [(w1, a), (w2, b)]

    nalt = k - 2
    erows = {}

    def escape(mi, r):
        key = (mi, r)
        if key in erows:
            return erows[key]
        j, a0, a1 = rays[r]
        m = mats[mi]
        coef = {}
        for side, rho in ((0, a0), (1, a1), (j, Q(1))):
            d = apply(m, N[side])
            parts = expansion(d) if specs[mi][0] == "vrot" else [(d, Q(1))]
            for w, c in parts:
                if w not in label:
                    fail("escape direction is not a contact direction")
                lab = label[w]
                x, y = coef.get(lab, (0, 0))
                coef[lab] = (x + rho * c * w[0], y + rho * c * w[1])
        sx = sum(v[0] for v in coef.values())
        sy = sum(v[1] for v in coef.values())
        if (sx, sy) != (0, 0):
            fail("unbalanced escape row")
        erows[key] = (coef, H[j])
        return erows[key]

    # tree walk: collect leaf obligations
    with gzip.open(directory / "tree.json.gz", "rt") as f:
        tree = json.load(f)
    obligations = {}

    def walk(node, assigned, rect, depth):
        if depth > 10000:
            fail("tree too deep")
        kind = node[0]
        if kind == "l":
            if node[1] in obligations:
                fail("leaf reused")
            obligations[node[1]] = (assigned, rect)
        elif kind == "p":
            mi = node[1]
            if type(mi) is not int or not 0 <= mi < len(specs):
                fail("bad placement index")
            if len(node) != 2 + nalt:
                fail("placement node must have one child per escape alternative")
            for r in range(nalt):
                walk(node[2 + r], assigned | {(mi, r)}, rect, depth + 1)
        elif kind == "s":
            plo, phi, qlo, qhi = rect
            cells = node[1]
            covered = set()
            for cell, child in cells:
                a, b, c, d = cell
                if not (plo <= a <= b <= phi and qlo <= c <= d <= qhi):
                    fail("split cell outside parent")
                covered |= {(p, q) for p in range(a, b + 1) for q in range(c, d + 1)}
                walk(child, assigned, tuple(cell), depth + 1)
            if len(covered) != (phi - plo + 1) * (qhi - qlo + 1):
                fail("split does not cover parent rectangle")
        else:
            fail("unknown node")

    walk(tree, frozenset(), (0, kR, 0, kL), 0)
    floor = n * n
    checked = 0
    seen = set()

    def stream():
        for name in manifest["files"]:
            if name.endswith(".xz"):
                with lzma.open(directory / name, "rt") as f:
                    for line in f:
                        idx, den, dk, vals = json.loads(line)
                        if len(dk) != len(vals):
                            fail("compact leaf record malformed")
                        keys, acc = [], 0
                        for d in dk:
                            if type(d) is not int or (keys and d <= 0):
                                fail("compact keys must be strictly increasing")
                            acc += d
                            keys.append(acc)
                        yield [idx, den, [list(p) for p in zip(keys, vals)]]
            else:
                with gzip.open(directory / name, "rt") as f:
                    for line in f:
                        yield json.loads(line)

    for index, den, pairs in stream():
        if index in seen or index not in obligations:
            fail("leaf set mismatch")
        seen.add(index)
        assigned, (plo, phi, qlo, qhi) = obligations.pop(index)
        if type(den) is not int or den <= 0:
            fail("bad denominator")
        force = [[Q(0), Q(0)] for _ in range(n)]
        B = Q(0)
        for key, num in pairs:
            if type(key) is not int or type(num) is not int or num <= 0:
                fail("multiplier must be a positive integer numerator")
            mu = Q(num, den)  # only ratios matter (homogeneous check)
            if 0 <= key < floor:
                d_, p_ = divmod(key, n)
                if d_ == p_:
                    fail("diagonal support row")
                u = normals[d_]
                force[d_][0] += mu * u[0]; force[d_][1] += mu * u[1]
                force[p_][0] -= mu * u[0]; force[p_][1] -= mu * u[1]
            elif key == floor:
                force[1][0] += mu
                force[0][0] -= mu
            else:
                e = key - floor - 1
                mi, r = divmod(e, nalt)
                if not 0 <= mi < len(specs) or (mi, r) not in assigned:
                    fail("escape row not assumed on this branch")
                coef, hj = escape(mi, r)
                for lab, v in coef.items():
                    force[lab][0] += mu * v[0]; force[lab][1] += mu * v[1]
                B += mu * hj
        # Homogeneous form of the leaf conditions (4.7): B > 0 and |sum_I f|^2 <= B^2.
        # Dividing every multiplier by B gives exactly B = 1, |sum_I f| <= 1.
        if B <= 0:
            fail("leaf %d: B <= 0" % index)
        B2 = B * B
        if sum(v[0] for v in force) != 0 or sum(v[1] for v in force) != 0:
            fail("leaf %d: forces do not sum to zero" % index)

        def ok(v):
            return v[0] * v[0] + v[1] * v[1] <= B2

        Lc = [(Q(0), Q(0))]
        for v in force[3 + kR:]:
            Lc.append((Lc[-1][0] + v[0], Lc[-1][1] + v[1]))
        Rc = [(Q(0), Q(0))]
        for v in force[3:3 + kR]:
            Rc.append((Rc[-1][0] + v[0], Rc[-1][1] + v[1]))
        sets = []
        for a in range(qlo + 1, kL + 1):          # {L_a..L_b}
            for b in range(a, kL + 1):
                sets.append((Lc[b][0] - Lc[a - 1][0], Lc[b][1] - Lc[a - 1][1]))
        for j in range(1, qhi + 2):               # {P-, L_j..L_kL}
            jj = min(j, kL + 1)
            sets.append((Lc[kL][0] - Lc[jj - 1][0] + force[0][0], Lc[kL][1] - Lc[jj - 1][1] + force[0][1]))
        for s in range(plo, kR + 1):              # {P+, R_1..R_s}
            sets.append((force[1][0] + Rc[s][0], force[1][1] + Rc[s][1]))
        for a in range(0, phi + 1):               # {R_{a+1}..R_b}
            for b in range(a, phi + 1):
                sets.append((Rc[b][0] - Rc[a][0], Rc[b][1] - Rc[a][1]))
        for v in sets:
            if not ok(v):
                fail("leaf %d: prefix norm exceeds B" % index)
        checked += 1
    if obligations or checked != header.get("leaves", checked):
        fail("tree leaves without multipliers")
    print("ACCEPT: %d leaves, %d placements, kR=%d kL=%d, k=%d sides" % (checked, len(specs), kR, kL, k))
    print("area(K) = %s = %.15f" % (area, float(area)))
    return area


if __name__ == "__main__":
    main(sys.argv[1])

