"""Search model (numerical proposals) for k-gon covers; exact acceptance via exactk.

Port of scripts/upper/search_model.py (ethan-keller/moser-worm-improved-bounds,
Apache-2.0) with: k-2 escape alternatives per placement, direct Clarabel calls
(no cvxpy), optional appended vrot placements, and adaptive orientation proposals.
"""
from fractions import Fraction
from functools import lru_cache
from math import cos, gcd, isfinite, radians, atan2, tan

import clarabel  # pyright: ignore[reportMissingImports]
import numpy as np
from scipy import sparse  # pyright: ignore[reportMissingImports]

from exactk import Geometry

SUPPORT_ANGLE_DEGREES = 90.0
ANGLE_TOLERANCE_DEGREES = 1e-9
PLACEMENT_TOLERANCE = 1e-9
BISECTION_IMPROVEMENT = 1e-4
CANDIDATE_TOLERANCE = 1e-6
COEFFICIENT_DROP = 1e-10
ROUNDING_DIGITS = 22
SOLVE_SECONDS = 30.0
MAX_ITERATIONS = 300
TIGHT_MAX_ITERATIONS = 50000


def partition(rect, right, left):
    def intervals(lo, hi, cuts):
        endpoints = [lo, *(i + 1 for i in sorted(set(cuts)) if lo <= i < hi), hi + 1]
        return [(a, b - 1) for a, b in zip(endpoints, endpoints[1:])]

    return [(a, b, c, d) for a, b in intervals(rect[0], rect[1], right) for c, d in intervals(rect[2], rect[3], left)]


class Model:
    def __init__(self, T, H, specs, adaptive=False):
        self.geometry = g = Geometry(T, H, specs)
        self.nalt = g.nalt
        self.adaptive = adaptive
        self.normals = np.array(g.normals, dtype=float)
        self.static_keys = []
        threshold = cos(radians(SUPPORT_ANGLE_DEGREES + ANGLE_TOLERANCE_DEGREES))
        for d in range(g.n):
            for p in range(g.n):
                if d != p and self.normals[d] @ self.normals[p] >= threshold:
                    self.static_keys.append(d * g.n + p)
        self.static_keys.append(g.n ** 2)
        self.H = np.array([float(h) for h in g.H])
        self._build_rows()
        self.cumulative = []
        for count, offset in ((g.kL, 3 + g.kR), (g.kR, 3)):
            rows, cols = [], []
            for j in range(count + 1):
                for i in range(j):
                    rows.append(j)
                    cols.append(offset + i)
            self.cumulative.append(sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(count + 1, g.n)))
        self.E = sparse.vstack(self.cumulative).tocsr()

    def _build_rows(self):
        g = self.geometry
        self.keys = self.static_keys + [g.n ** 2 + 1 + i for i in range(len(g.erows))]
        rows = [g.row(k) for k in self.keys]
        self.row_labels = [set(i for i, _ in terms) for terms, _ in rows]
        self.mx, self.my = self.coefficient_matrices(rows)
        self.constants = np.array([float(b) for _, b in rows])
        fit_rows, fit_cols, fit_values = [], [], []
        for r, (terms, bound) in enumerate(g.erows):
            for label, v in terms:
                weight = sum(a * b for a, b in zip(v, g.normals[label]))
                if weight < 0:
                    raise ValueError("negative support expansion")
                fit_rows.append(r)
                fit_cols.append(label)
                fit_values.append(float(weight / -bound))
        self.fit = sparse.csr_matrix((fit_values, (fit_rows, fit_cols)), shape=(len(g.erows), g.n))

    def sync(self, extra):
        """Append vrot placements beyond the base list (shared indices)."""
        g = self.geometry
        base = getattr(self, "_base", None)
        if base is None:
            self._base = base = len(g.specs)
        have = len(g.specs) - base
        if len(extra) > have:
            for spec in extra[have:]:
                g.add_vrot(spec)
            self._build_rows()

    def coefficient_matrices(self, rows):
        ri, ci, xs, ys = [], [], [], []
        for j, (terms, _) in enumerate(rows):
            for i, (x, y) in terms:
                ri.append(i)
                ci.append(j)
                xs.append(float(x))
                ys.append(float(y))
        shape = (self.geometry.n, len(rows))
        return sparse.csc_matrix((xs, (ri, ci)), shape=shape), sparse.csc_matrix((ys, (ri, ci)), shape=shape)

    @lru_cache(maxsize=64)
    def family_matrix(self, rect):
        g = self.geometry
        a, b = g.n, g.n + g.kL + 1
        lo, hi, qlo, qhi = rect
        rows, cols, values = [], [], []
        count = 0

        def add(entries):
            nonlocal count
            for col, value in entries:
                rows.append(count)
                cols.append(col)
                values.append(value)
            count += 1

        for i in range(qlo, g.kL):
            for j in range(i + 1, g.kL + 1):
                add([(a + j, 1), (a + i, -1)])
        for j in range(qhi + 1):
            add([(a + g.kL, 1), (a + j, -1), (0, 1)])
        for j in range(lo, g.kR + 1):
            add([(a + g.kL, 1), (0, 1), (2, 1), (b + g.kR, 1), (b + j, -1)])
        for i in range(hi + 1):
            for j in range(i, hi + 1):
                add([(a + g.kL, 1), (0, 1), (2, 1), (1, 1), (b + i, 1), (b + g.kR, 1), (b + j, -1)])
        return sparse.csr_matrix((values, (rows, cols)), shape=(count, g.n + g.kL + g.kR + 2))

    def settings(self, tight, variant):
        s = clarabel.DefaultSettings()
        s.verbose = False
        tol = 1e-11 if tight else 1e-8
        s.max_iter = TIGHT_MAX_ITERATIONS if tight else MAX_ITERATIONS
        s.time_limit = SOLVE_SECONDS
        s.max_threads = 1
        s.direct_solve_method = "qdldl"
        s.tol_gap_abs = tol
        s.tol_gap_rel = tol
        s.tol_feas = tol
        if tight:
            s.tol_ktratio = 1e-10
        if variant == 1:
            s.static_regularization_constant = 1e-7
        elif variant == 2:
            s.static_regularization_constant = 1e-6
            s.equilibrate_enable = False
        return s

    def solve(self, assigned, rect, *, tight=False, all_rows=False, subset=None):
        g = self.geometry
        nalt = self.nalt
        escapes = [len(self.static_keys) + nalt * m + ray for m, ray in assigned]
        involved = set().union(*(self.row_labels[i] for i in escapes)) if escapes else set()
        if subset is not None:
            selected = list(subset)
        else:
            selected = [i for i, key in enumerate(self.static_keys)
                        if all_rows or key == g.n ** 2 or self.row_labels[i] & involved] + escapes
        m = len(selected)
        n = g.n
        N = n + g.kL + g.kR + 2
        Ax = self.mx[:, selected]
        Ay = self.my[:, selected]
        fam = self.family_matrix(tuple(rect))
        R = fam.shape[0]
        I_n = sparse.identity(n, format="csc")
        Z = lambda r, c: sparse.csc_matrix((r, c))
        # variable order: mu (m), wx (N), wy (N)
        Ew = sparse.hstack([-self.E, sparse.identity(N - n, format="csc")])  # rows N-n, cols N
        blocks = [
            [-Ax, sparse.hstack([I_n, Z(n, N - n)]), Z(n, N)],
            [-Ay, Z(n, N), sparse.hstack([I_n, Z(n, N - n)])],
            [Z(N - n, m), Ew, Z(N - n, N)],
            [Z(N - n, m), Z(N - n, N), Ew],
            [-sparse.identity(m, format="csc"), Z(m, N), Z(m, N)],
        ]
        A_top = sparse.bmat(blocks, format="csc")
        # SOC rows interleaved: (t=1, famx, famy) per family row
        famc = fam.tocoo()
        rr = np.concatenate([3 * famc.row + 1, 3 * famc.row + 2])
        cc = np.concatenate([m + famc.col, m + N + famc.col])
        vv = np.concatenate([-famc.data, -famc.data])
        A_soc = sparse.csc_matrix((vv, (rr, cc)), shape=(3 * R, m + 2 * N))
        A = sparse.vstack([A_top, A_soc], format="csc")
        nz = 2 * n + 2 * (N - n)
        b = np.zeros(A.shape[0])
        b[nz + m::3] = 1.0
        cones = [clarabel.ZeroConeT(nz), clarabel.NonnegativeConeT(m)] + [clarabel.SecondOrderConeT(3)] * R
        q = np.concatenate([self.constants[selected], np.zeros(2 * N)])
        P = sparse.csc_matrix((m + 2 * N, m + 2 * N))
        for variant in range(3):
            try:
                solver = clarabel.DefaultSolver(P, q, A, b, cones, self.settings(tight, variant))
                sol = solver.solve()
            except Exception:
                continue
            st = str(sol.status)
            if "Solved" not in st:
                continue
            x = np.array(sol.x)
            zd = np.array(sol.z)
            mu = x[:m]
            value = float(-self.constants[selected] @ mu)
            points = -np.column_stack([zd[:n], zd[n:2 * n]])
            if not isfinite(value) or not np.isfinite(points).all() or not np.isfinite(mu).all():
                continue
            return value, mu, points, [self.keys[i] for i in selected]
        return None

    def round_leaf(self, keys, values, rect):
        """Small scale c with round(c*mu) passing the exact check: x4 search, then sqrt2 refine."""
        rows = self.geometry.integer_rows(keys)

        def ok(c):
            weights = [int(round(v * c)) for v in values]
            return weights if any(weights) and self.geometry.valid(weights, rows, rect) else None

        c = 4.0
        found = None
        while c < 1e23:
            found = ok(c)
            if found is not None:
                break
            c *= 4.0
        if found is None:
            return None, rows
        for f in (2.0, 1.4142135623730951):  # try c/2 and c/sqrt2 (search between c/4 and c)
            w = ok(c / f)
            if w is not None:
                return w, rows
        return found, rows

    def drop_small(self, keys, weights, rows, rect):
        mag = [w * sum(abs(x) + abs(y) for _, x, y in coeffs) for w, (coeffs, b) in zip(weights, rows)]
        order = sorted((i for i in range(len(weights)) if weights[i]), key=lambda i: mag[i])
        lo, hi, best = 0, len(order), 0
        while lo < hi:
            m = (lo + hi + 1) // 2
            drop = set(order[:m])
            ww = [0 if i in drop else w for i, w in enumerate(weights)]
            if self.geometry.valid(ww, rows, rect):
                best, lo = m, m
            else:
                hi = m - 1
        drop = set(order[:best])
        return [0 if i in drop else w for i, w in enumerate(weights)]

    def finish_leaf(self, prop, rect):
        _, values, _, kk = prop
        items = [(k, float(v)) for k, v in zip(kk, values) if v > COEFFICIENT_DROP]
        if not items:
            return None
        ks, vs = zip(*items)
        weights, rows = self.round_leaf(ks, vs, rect)
        if weights is None:
            return None
        weights = self.drop_small(ks, weights, rows, rect)
        common = 0
        for x in weights:
            common = gcd(common, x)
        # compact homogeneous form (den = 1): checker verifies B > 0, |sum_I f|^2 <= B^2
        return [1, [[k, x // common] for k, x in zip(ks, weights) if x]]

    def exact_leaf(self, proposal, rect, assigned=None):
        if proposal is None or proposal[0] < 1 - CANDIDATE_TOLERANCE:
            return None
        value, mu, _, keys = proposal
        if assigned is not None and len(mu) > 0:
            top = float(np.max(mu))
            pos = {k: i for i, k in enumerate(self.keys)}
            last = len(mu)
            for tau in (0.1, 0.03, 0.01):  # sparse re-solves, smallest row set first
                sub = [pos[k] for k, v in zip(keys, mu) if v >= tau * top]
                if len(sub) >= last:
                    continue
                last = len(sub)
                p2 = self.solve(assigned, rect, subset=sub)
                if p2 is not None and p2[0] >= 1 - CANDIDATE_TOLERANCE:
                    leaf = self.finish_leaf(p2, rect)
                    if leaf is not None:
                        return leaf
        return self.finish_leaf(proposal, rect)

    def placement_values(self, points):
        support = (points @ self.normals.T).max(axis=0)
        return np.asarray(self.fit @ support).reshape(-1, self.nalt).max(axis=1)

    def best_placement(self, assigned, points):
        values = self.placement_values(points)
        for m, _ in assigned:
            values[m] = np.inf
        index = int(values.argmin())
        return index if values[index] <= 1 + PLACEMENT_TOLERANCE else None

    def new_labels(self, assigned, placement, ray):
        offset = len(self.static_keys)
        nalt = self.nalt
        involved = set().union(*(self.row_labels[offset + nalt * m + r] for m, r in assigned)) if assigned else set()
        new = self.row_labels[offset + nalt * placement + ray] - involved
        right = sorted(i - 3 for i in new if 3 <= i < 3 + self.geometry.kR)
        left = sorted(i - 3 - self.geometry.kR for i in new if i >= 3 + self.geometry.kR)
        return right, left

    # ---- adaptive orientation: continuous search with vrot expansions ----
    def vrot_ratio(self, points, thetas, mirror):
        """Escape LHS/h ratio (max over alternatives) of vrot placements at angles thetas."""
        g = self.geometry
        P = points
        sup = (P @ self.normals.T).max(axis=0)  # support per label direction
        ang = np.arctan2(self.normals[:, 1], self.normals[:, 0]) % (2 * np.pi)
        order = np.argsort(ang, kind="stable")
        # drop duplicate directions (labels 0 and 1 are both downward): keep label 1
        uniq = []
        for idx in order:
            if uniq and abs(ang[idx] - ang[uniq[-1]]) < 1e-15:
                if idx == 1:
                    uniq[-1] = idx
                continue
            uniq.append(idx)
        uang = ang[uniq]
        KN = np.array([[float(x), float(y)] for x, y in g.KN])
        out = np.full(len(thetas), -np.inf)
        a = np.array([[float(x), float(y)] for x, y in g.RAYS])
        Hj = np.array([float(h) for h in g.H[2:]])
        c, s = np.cos(thetas), np.sin(thetas)
        hv = []
        for i in range(g.k):
            dx = c * KN[i, 0] - s * KN[i, 1]
            dy = s * KN[i, 0] + c * KN[i, 1]
            if mirror:
                dx = -dx
            phi = np.arctan2(dy, dx) % (2 * np.pi)
            pos = np.searchsorted(uang, phi, side="right")
            i2 = pos % len(uniq)
            i1 = (pos - 1) % len(uniq)
            w1 = self.normals[np.array(uniq)[i1]]
            w2 = self.normals[np.array(uniq)[i2]]
            det = w1[:, 0] * w2[:, 1] - w1[:, 1] * w2[:, 0]
            aa = (dx * w2[:, 1] - dy * w2[:, 0]) / det
            bb = (w1[:, 0] * dy - w1[:, 1] * dx) / det
            hv.append(aa * sup[np.array(uniq)[i1]] + bb * sup[np.array(uniq)[i2]])
        hv = np.array(hv)  # (k, G)
        ratios = []
        for r in range(g.nalt):
            ratios.append((a[r, 0] * hv[0] + a[r, 1] * hv[1] + hv[2 + r]) / Hj[r])
        return np.max(np.array(ratios), axis=0)

    def propose_orientation(self, points, existing):
        """Find a vrot placement (rational t, mirror) excluding the given points, or None."""
        best = None
        for mirror in (0, 1):
            th = np.linspace(-np.pi, np.pi, 7200, endpoint=False)
            r = self.vrot_ratio(points, th, mirror)
            i = int(r.argmin())
            t0, v0 = th[i], r[i]
            step = th[1] - th[0]
            for _ in range(3):
                tt = t0 + np.linspace(-step, step, 81)
                rr = self.vrot_ratio(points, tt, mirror)
                j = int(rr.argmin())
                if rr[j] < v0:
                    t0, v0 = tt[j], rr[j]
                step /= 40
            if best is None or v0 < best[0]:
                best = (v0, t0, mirror)
        v0, t0, mirror = best
        if v0 > 1 - 1e-7:
            return None
        # rational t = tan(theta/2) with small denominator, re-checked with exact rows
        g = self.geometry
        sup = (points @ self.normals.T).max(axis=0)
        tval = tan(t0 / 2)
        for den in (50, 200, 1000, 2000, 10000, 100000, 10 ** 6):
            t = Fraction(tval).limit_denominator(den)
            spec = ["vrot", str(t), mirror]
            if spec in existing:
                continue
            worst = 0.0
            for r in range(g.nalt):
                terms, bound = g.escape(spec, r)
                val = 0.0
                for label, v in terms:
                    wgt = float(sum(a * b for a, b in zip(v, g.normals[label])))
                    val += wgt * sup[label]
                worst = max(worst, val / float(-bound))
            if worst <= 1 - 1e-8:
                return spec
        return None

    def expand(self, state, existing_specs=None):
        assigned, rect, right, left = state
        proposal = self.solve(assigned, rect)
        leaf = self.exact_leaf(proposal, rect, assigned)
        if leaf is None and proposal is not None and proposal[0] >= 1 - CANDIDATE_TOLERANCE:
            for kw in ({"tight": True}, {"all_rows": True}):
                leaf = self.exact_leaf(self.solve(assigned, rect, **kw), rect, assigned)
                if leaf is not None:
                    break
        if leaf is not None:
            return ["l", *leaf]
        cells = partition(rect, right, [])
        if len(cells) > 1:
            return ["s", [[cell, [assigned, cell, [], left]] for cell in cells]]
        cells = partition(rect, [], left)
        if len(cells) > 1:
            return ["s", [[cell, [assigned, cell, [], []]] for cell in cells]]
        if proposal is None:
            return ["unresolved", "solver"]
        placement = self.best_placement(assigned, proposal[2])
        if placement is None:
            for axis in sorted((0, 2), key=lambda i: rect[i] - rect[i + 1]):
                if rect[axis] == rect[axis + 1]:
                    continue
                a, b = list(rect), list(rect)
                mid = (rect[axis] + rect[axis + 1]) // 2
                a[axis + 1], b[axis] = mid, mid + 1
                trials = [self.solve(assigned, cell) for cell in (a, b)]
                if any(p is not None and p[0] > proposal[0] + BISECTION_IMPROVEMENT for p in trials):
                    return ["s", [[cell, [assigned, cell, [], []]] for cell in (a, b)]]
            if self.adaptive:
                spec = self.propose_orientation(proposal[2], existing_specs or self.geometry.specs)
                if spec is not None:
                    return ["newspec", spec, proposal[0]]
            return ["unresolved", "no excluding placement", proposal[0], proposal[2].tolist()]
        children = []
        for ray in range(self.nalt):
            rnew, lnew = self.new_labels(assigned, placement, ray)
            children.append([assigned + [[placement, ray]], rect, rnew, lnew])
        return ["p", placement, *children]

