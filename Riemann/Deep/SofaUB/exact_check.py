"""Exact rational re-computation of the box bound integral for sampled leaves (independent of the C core).

For a leaf, the C core reports the center slope kc, center Xc = X(ya), edge vectors g1, g2 (interval enclosures),
and push-outs PA0, PA1, PB0, PB1. The dilated slice is
    D(y) = {p : |cross(p - X(y), g2)| <= A(y), |cross(g1, p - X(y))| <= B(y)},
    X(y) = Xc + (y - ya) kc,  A(y) = det/2 + PA0 + PA1 |y - ya|,  B(y) = det/2 + PB0 + PB1 |y - ya|.
Here we take g1, g2 = interval midpoints as exact rationals (det = cross(g1, g2) exactly) and integrate
area(D(y) ^ R_y) over y in [0, c] EXACTLY: on each sub-range the vertices of D(y) move affinely, so the
area is piecewise quadratic; all candidate breakpoints are rational and Simpson's rule is exact on each piece.
"""
from fractions import Fraction as F


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def clip(poly, axis, sgn, h):
    """keep sgn*(coord - h) <= 0"""
    out = []
    n = len(poly)
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        da, db = sgn * (a[axis] - h), sgn * (b[axis] - h)
        if da <= 0:
            out.append(a)
        if (da < 0 < db) or (db < 0 < da):
            t = (h - a[axis]) / (b[axis] - a[axis])
            q = [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
            q[axis] = h
            out.append((q[0], q[1]))
    return out


def area_poly(p):
    s = F(0)
    for i in range(len(p)):
        j = (i + 1) % len(p)
        s += p[i][0] * p[j][1] - p[j][0] * p[i][1]
    return abs(s) / 2


def area_DR(V, reg):
    P = clip(V, 0, 1, F(1))
    if len(P) < 3: return F(0)
    P = clip(P, 1, -1, F(0))
    if len(P) < 3: return F(0)
    P = clip(P, 1, 1, F(1))
    if len(P) < 3: return F(0)
    if reg == 0:
        P = clip(P, 0, -1, F(0))
        if len(P) < 3: return F(0)
    return area_poly(P)


class Dfam:
    def __init__(self, kc, xc, g1, g2, PA0, PA1, PB0, PB1, ya):
        self.kc = [F(v) for v in kc]; self.xc = [F(v) for v in xc]
        self.g1 = [F(v) for v in g1]; self.g2 = [F(v) for v in g2]
        self.det = cross(self.g1, self.g2)
        if not self.det > 0:
            raise RuntimeError("edge vectors not counter-clockwise")
        self.PA0, self.PA1, self.PB0, self.PB1 = F(PA0), F(PA1), F(PB0), F(PB1)
        self.ya = F(ya)

    def verts(self, y):
        d = y - self.ya
        ad = abs(d)
        X = (self.xc[0] + d * self.kc[0], self.xc[1] + d * self.kc[1])
        al = F(1, 2) + (self.PA0 + self.PA1 * ad) / self.det
        be = F(1, 2) + (self.PB0 + self.PB1 * ad) / self.det
        g1, g2 = self.g1, self.g2
        out = []
        for s1, s2 in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            out.append((X[0] + s1 * al * g1[0] + s2 * be * g2[0], X[1] + s1 * al * g1[1] + s2 * be * g2[1]))
        return out

    def f(self, y, reg):
        return area_DR(self.verts(y), reg)

    def integ_sub(self, y0, y1, reg):
        # vertices affine on [y0,y1]: v(y) = v(y0) + (y-y0) v'
        V0, V1 = self.verts(y0), self.verts(y1)
        L = y1 - y0
        bps = {y0, y1}
        lines_x = [F(1), F(0)] if reg == 0 else [F(1)]
        for (a, b) in zip(V0, V1):
            vx, vz = (b[0] - a[0]) / L, (b[1] - a[1]) / L
            for xl in lines_x:
                if vx != 0:
                    t = y0 + (xl - a[0]) / vx
                    if y0 < t < y1: bps.add(t)
            for zl in (F(0), F(1)):
                if vz != 0:
                    t = y0 + (zl - a[1]) / vz
                    if y0 < t < y1: bps.add(t)
        corners = [(F(0), F(0)), (F(1), F(0)), (F(1), F(1)), (F(0), F(1))] if reg == 0 else [(F(1), F(0)), (F(1), F(1))]
        for e in range(4):
            a0, b0 = V0[e], V0[(e + 1) % 4]
            a1 = V1[e]
            n = (-(b0[1] - a0[1]), b0[0] - a0[0])        # edge normal (direction fixed in y)
            o0 = n[0] * a0[0] + n[1] * a0[1]
            o1 = n[0] * a1[0] + n[1] * a1[1]
            sl = (o1 - o0) / L
            if sl == 0: continue
            for p in corners:
                t = y0 + (n[0] * p[0] + n[1] * p[1] - o0) / sl
                if y0 < t < y1: bps.add(t)
        bl = sorted(bps)
        tot = F(0)
        for a, b in zip(bl[:-1], bl[1:]):
            m = (a + b) / 2
            tot += (b - a) * (self.f(a, reg) + 4 * self.f(m, reg) + self.f(b, reg)) / 6
        return tot

    def integral(self, c):
        pts = sorted({F(0), F(1), F(c), self.ya} if 0 < self.ya < F(c) else {F(0), F(1), F(c)})
        tot = F(0)
        for a, b in zip(pts[:-1], pts[1:]):
            reg = 1 if b <= 1 else 0
            tot += self.integ_sub(a, b, reg)
        return tot
