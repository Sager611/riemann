"""Exact geometry for the support/length method for k-gons K = {n_i.x <= h_i}
(h0 = h1 = 0; sides j >= 2 in the open negative cone of sides 0, 1).

Generalizes scripts/upper/exact.py of ethan-keller/moser-worm-improved-bounds
(Apache-2.0): for k = 4 the rows are identical. Escape alternative r <-> side 2 + r;
escape row key = n^2 + 1 + (k-2)*spec + r.
"""
from fractions import Fraction as Q
from functools import cmp_to_key
from math import lcm

ZERO = (Q(0), Q(0))


def add(u, v):
    return u[0] + v[0], u[1] + v[1]


def scale(a, u):
    return a * u[0], a * u[1]


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


def unit(t):
    return (1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)


def rays_for(normals):
    n0, n1 = normals[0], normals[1]
    det = cross(n0, n1)
    rays = []
    for j in range(2, len(normals)):
        v = (-normals[j][0], -normals[j][1])
        a, b = cross(v, n1) / det, cross(n0, v) / det
        if a <= 0 or b <= 0:
            raise ValueError("side %d not in the negative cone of sides 0,1" % j)
        rays.append((a, b))
    return rays


def sector(u):
    if u[1] == 0:
        return 0 if u[0] > 0 else 2
    return 1 if u[1] > 0 else 3


def angle_lt(u, v):
    a, b = sector(u), sector(v)
    return a < b if a != b else cross(u, v) > 0


def vertices(T, H):
    k = len(T)
    normals = [unit(Q(t)) for t in T]
    H = [Q(h) for h in H]
    vs = []
    for i in range(k):
        a, b = normals[i - 1], normals[i]
        ha, hb = H[i - 1], H[i]
        det = cross(a, b)
        vs.append(((ha * b[1] - a[1] * hb) / det, (a[0] * hb - ha * b[0]) / det))
    return vs  # vs[i] = line (i-1) cap line i


def area(T, H):
    vs = vertices(T, H)
    k = len(vs)
    return sum(cross(vs[i], vs[(i + 1) % k]) for i in range(k)) / 2


class Geometry:
    def __init__(self, T, H, specs):
        self.T = [Q(t) for t in T]
        self.H = [Q(h) for h in H]
        self.k = len(self.T)
        self.nalt = self.k - 2
        if self.H[0] != 0 or self.H[1] != 0:
            raise ValueError("sides 0,1 must pass through the origin")
        self.KN = [unit(t) for t in self.T]
        self.RAYS = rays_for(self.KN)
        for spec in specs:
            self.check_spec(spec)
        self.specs = list(specs)
        dirs = {u for s in self.specs if s[0] != "vrot" for u in self.directions(s)}
        self.right = sorted((u for u in dirs if u[0] > 0), key=lambda u: u[1])
        self.left = sorted((u for u in dirs if u[0] < 0), key=lambda u: -u[1])
        self.kR, self.kL = len(self.right), len(self.left)
        self.normals = [(Q(0), Q(-1)), (Q(0), Q(-1)), (Q(0), Q(1))]
        self.normals += self.right + self.left
        self.n = len(self.normals)
        self.labels = {u: i for i, u in enumerate(self.normals)}
        self.named = sorted(
            self.labels,
            key=cmp_to_key(lambda u, v: -1 if angle_lt(u, v) else (1 if angle_lt(v, u) else 0)),
        )
        self.erows = [self.escape(s, r) for s in self.specs for r in range(self.nalt)]
        self._rows = {}

    def check_spec(self, spec):
        if spec[0] == "flush":
            if len(spec) != 4 or spec[1] not in range(self.k) or any(x not in (0, 1) for x in spec[2:]):
                raise ValueError("invalid flush placement")
        elif spec[0] in ("rot", "vrot"):
            if len(spec) != 3 or spec[2] not in (0, 1):
                raise ValueError("invalid rotation placement")
            Q(spec[1])
        else:
            raise ValueError("unknown placement type")

    def add_vrot(self, spec):
        """Append a vrot placement (labels unchanged); returns its index."""
        if spec[0] != "vrot":
            raise ValueError("only vrot placements can be appended")
        self.check_spec(spec)
        self.specs.append(spec)
        for r in range(self.nalt):
            self.erows.append(self.escape(spec, r))
        return len(self.specs) - 1

    def matrix(self, spec):
        kind, *args = spec
        if kind == "flush":
            side, flip, mirror = args
            x, y = self.KN[side]
            mat = ((y, -x), (x, y)) if flip else ((-y, x), (-x, -y))
        else:
            t, mirror = args
            x, y = unit(Q(t))
            mat = ((x, -y), (y, x))
        if mirror:
            mat = ((-mat[0][0], -mat[0][1]), mat[1])
        return mat

    def directions(self, spec):
        a, b = self.matrix(spec)
        return [(a[0] * x + a[1] * y, b[0] * x + b[1] * y) for x, y in self.KN]

    def expansion(self, u):
        if u in self.labels:
            return [(u, Q(1))]
        i = next((i for i, v in enumerate(self.named) if angle_lt(u, v)), 0)
        w1, w2 = self.named[i - 1], self.named[i]
        det = cross(w1, w2)
        a, b = cross(u, w2) / det, cross(w1, u) / det
        if a < 0 or b < 0 or add(scale(a, w1), scale(b, w2)) != u:
            raise ValueError("invalid direction expansion")
        return [(w1, a), (w2, b)]

    def escape(self, spec, r):
        coeffs = []
        ds = self.directions(spec)
        a0, a1 = self.RAYS[r]
        for i, rho in ((0, a0), (1, a1), (2 + r, Q(1))):
            parts = self.expansion(ds[i]) if spec[0] == "vrot" else [(ds[i], Q(1))]
            for u, a in parts:
                coeffs.append((self.labels[u], scale(rho * a, u)))
        if (sum(v[0] for _, v in coeffs), sum(v[1] for _, v in coeffs)) != ZERO:
            raise ValueError("unbalanced escape row")
        return coeffs, -self.H[2 + r]

    def row(self, key):
        if key in self._rows:
            return self._rows[key]
        floor = self.n ** 2
        if key < floor:
            d, p = divmod(key, self.n)
            if d == p:
                raise ValueError("diagonal support row")
            u = self.normals[d]
            row = [(d, u), (p, scale(-1, u))], Q(0)
        elif key == floor:
            row = [(1, (Q(1), Q(0))), (0, (Q(-1), Q(0)))], Q(0)
        else:
            row = self.erows[key - floor - 1]
        self._rows[key] = row
        return row

    def integer_rows(self, keys):
        rows = [self.row(k) for k in keys]
        den = 1
        for coeffs, b in rows:
            den = lcm(den, b.denominator)
            for _, v in coeffs:
                den = lcm(den, v[0].denominator, v[1].denominator)

        def integer(q):
            return q.numerator * (den // q.denominator)

        return [([(i, integer(v[0]), integer(v[1])) for i, v in coeffs], integer(-b)) for coeffs, b in rows]

    def valid(self, weights, rows, rect):
        force = [[0, 0] for _ in self.normals]
        bound = 0
        for weight, (coeffs, b) in zip(weights, rows):
            if weight < 0:
                return False
            if not weight:
                continue
            bound += weight * b
            for i, x, y in coeffs:
                force[i][0] += weight * x
                force[i][1] += weight * y
        if bound <= 0 or sum(v[0] for v in force) or sum(v[1] for v in force):
            return False
        bound2 = bound * bound

        def good(v):
            return v[0] * v[0] + v[1] * v[1] <= bound2

        a = [(0, 0)]
        b = [(0, 0)]
        for v in force[3 + self.kR:]:
            a.append(add(a[-1], v))
        for v in force[3: 3 + self.kR]:
            b.append(add(b[-1], v))
        lo, hi, qlo, qhi = rect

        def sub(u, v):
            return u[0] - v[0], u[1] - v[1]

        return (
            all(good(sub(a[j], a[i])) for i in range(qlo, self.kL) for j in range(i + 1, self.kL + 1))
            and all(good(add(sub(a[-1], a[j]), force[0])) for j in range(min(qhi, self.kL) + 1))
            and all(good(add(force[1], b[j])) for j in range(lo, self.kR + 1))
            and all(good(sub(b[j], b[i])) for i in range(hi + 1) for j in range(i, hi + 1))
        )

