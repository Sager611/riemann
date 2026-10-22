#!/usr/bin/env python3
"""Exact margin checks for the kappa-dependent steps of OpenAI 7/8, Lemma 18.1.

(1) Paper (TeX l. 12947-12984): centered range A > 5M/6 with A + (6k-1)z <= M.
    sup z/M = 1/(6(6k-1)); A_comp <= 3M/2 - A + 2z; A_comp + (6k-1)z <= 5M/2 - 2A + 2z.
(2) Lean port, robust gate (Energy/LiveClippingDefect.lean, robust_reflected_low and
    robust_reflected_low_admissible): exact LP by vertex enumeration over the rationals,
    with clipping defect e and reflection slack xi (the port uses e = xi = M0/56).
Pure stdlib, runs in < 1 s, about 16 MB.
"""
import itertools
from fractions import Fraction as F

Row = list[F]


def solve(rows: list[Row], rhs: Row) -> Row | None:
    n = len(rows)
    m: list[Row] = [r[:] + [b] for r, b in zip(rows, rhs)]
    for c in range(n):
        p = next((r for r in range(c, n) if m[r][c] != 0), None)
        if p is None:
            return None
        m[c], m[p] = m[p], m[c]
        for r in range(n):
            if r != c and m[r][c] != 0:
                f = m[r][c] / m[c][c]
                m[r] = [a - f * b for a, b in zip(m[r], m[c])]
    return [m[i][n] / m[i][i] for i in range(n)]


def dot(u: Row, v: Row) -> F:
    return sum((a * b for a, b in zip(u, v)), F(0))


def lp_max(obj: Row, a: list[Row], b: Row) -> F:
    best: F | None = None
    for s in itertools.combinations(range(len(a)), len(obj)):
        x = solve([a[i] for i in s], [b[i] for i in s])
        if x is None:
            continue
        if all(dot(a[i], x) <= b[i] for i in range(len(a))):
            val = dot(obj, x)
            if best is None or val > best:
                best = val
    if best is None:
        raise ValueError("infeasible or unbounded vertex set")
    return best


def paper(k: F) -> tuple[F, F, F, F]:
    q = 6 * k - 1
    zs = F(1, 6) / q
    a1 = F(2, 3) + 2 * zs
    a2 = F(5, 6) + F(1, 3) / q
    return q, zs, a1, a2


def lean_gate(k: F, e: F, xi: F) -> tuple[F, F]:
    # variables: short, along, z, ell, total, refl   (width M = 1)
    k6 = 6 * k
    rows: list[Row] = []
    rhs: Row = []

    def le(c: list[F | int], v: F | int) -> None:
        rows.append([F(x) for x in c])
        rhs.append(F(v))

    le([1, 0, 0, 0, 0, 0], F(1, 4))          # short <= M/4
    le([-1, 0, 0, 0, 0, 0], 0)               # short >= 0
    le([0, 0, 1, -1, 0, 0], 0)               # z <= ell
    le([0, 0, -1, 0, 0, 0], 0)               # z >= 0
    le([0, 0, 0, -1, 0, 0], 0)               # ell >= 0
    le([0, 0, 0, -1, -1, 0], -(F(5, 6) - e))  # 5M/6 - e <= total + ell
    le([0, 0, 0, k6, 1, 0], 1)               # total + ell + (6k-1) ell <= M
    le([-1, -1, -1, 0, 0, 0], -F(5, 6))      # 5M/6 <= short + along + z (closure of <)
    le([0, 1, 0, 0, 0, 1], 1 + xi)           # refl <= M - along + xi
    le([0, 1, 0, 0, -1, 0], -F(1, 4))        # along <= total - M/4
    le([0, -1, 0, 0, 0, 0], 0)               # along >= 0
    one = [F(1), F(0), F(1), F(0), F(0), F(1)]
    six = [F(1), F(0), k6, F(0), F(0), F(1)]
    return lp_max(one, rows, rhs), lp_max(six, rows, rhs)


def main() -> None:
    kflat = F(2 * 874886647040 - 10**12, 10**12)
    cases: list[tuple[str, F]] = [("13/18", F(13, 18)), ("Xu-Sang k1", F(29991, 40000)),
                                  ("Liu k_flat", kflat), ("3/4", F(3, 4))]
    for name, k in cases:
        q, zs, a1, a2 = paper(k)
        print(f"{name:11s} k={float(k):.6f} 6k-1={float(q):.5f} z<M/{float(1 / zs):.3f} "
              f"Acomp<={a1}={float(a1):.5f} (5/6 margin {float(F(5, 6) - a1):.5f}) "
              f"Acomp+(6k-1)z<={float(a2):.5f} (margin {float(1 - a2):.5f})")
    for k in [F(3, 4), F(13, 18)]:
        for e in [F(0), F(1, 56)]:
            g1, g2 = lean_gate(k, e, e)
            print(f"Lean robust gate k={k} e=xi={e}: max refl+short+z={g1} "
                  f"({float(g1):.5f} <= 5/6: {g1 <= F(5, 6)}), "
                  f"max refl+short+6kz={g2} ({float(g2):.5f} <= 1: {g2 <= 1})")
    print("general comparison margin (q-2)/(6q), q=6k-1: positive iff k>1/2")


if __name__ == "__main__":
    main()
