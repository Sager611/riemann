#!/usr/bin/env python3
"""Exact LP for the reflected-energy exponent of OpenAI 7/8, Lemma 14.3, and the 13/15 floor.

Lemma 14.3 (as used in Liu 2610.12234, Lemmas 16.1/22.1): for rows of length m and remaining
marked prime length l (theta length Z^(1+l)), local exponents
  O,H,dH,A0,N0,S0,B0,za,v,lb,el >= 0,  O+H+dH = m,  2A0 <= O,  N0+B0+S0 <= A0 (or N0 <= A0 only),
  za <= l,  Td = 2H+2A0+2za-(1+l)-N0-3B0,  v+3lb+el <= Td,
  Eref = O/2 + max(H, v+lb) - S0 - B0 + [za-u] - lb - 2el/3 - (Td-v-3lb-el)/2,
  za-u = max(za-v, 0, (2za-v)/3)   (cubic large sieve factor N + L + (NL)^(2/3), Lemma 14.2).
k_LP(m,l) = max(Eref) - m.  Liu's closed form: k = max(0, (1+3l-2m)/4, m+l-1, (3l-1)/2).

Variant "trilinear" (sketch; not in the papers): the active slot tuple coefficient carries the
mark-mark phase (q/p)_3^2 (OpenAI (14.5) with j_p = j_q = 0).  Split the active slots P = P1*P2;
by reciprocity the P2-phase is conj((P1*n/P2)_3).  Cauchy over P1, then Heath-Brown's cubic large
sieve (OpenAI (14.6)) in (m = P1*n, P2) gives N*L1 + L2 + (N*L)^(2/3); the split L1 = (L/N)^(1/2)
makes this max(N, (NL)^(2/3)): the L term of Lemma 14.2 drops, i.e. za-u -> max(0, (2za-v)/3).

Run: python3 reflected_lp.py   (numpy only, ~10 s, <60 MB)
"""
import itertools
import numpy as np

V = dict(O=0, H=1, dH=2, A0=3, N0=4, S0=5, B0=6, za=7, v=8, lb=9, el=10, t=11)
NV = len(V)


def simplex_max(c, Aub, bub, Aeq, beq, tol=1e-11):
    """max c.x, Aub x <= bub, Aeq x = beq, x >= 0 (two-phase tableau, Bland's rule)."""
    rows = [(np.array(a, float), float(b), 1) for a, b in zip(Aub, bub)]
    rows += [(np.array(a, float), float(b), 0) for a, b in zip(Aeq, beq)]
    n, m = len(c), len(rows)
    ns = sum(r[2] for r in rows)
    N = n + ns + m
    T = np.zeros((m, N + 1))
    basis, si = [], 0
    for i, (a, b, ub) in enumerate(rows):
        T[i, :n] = a
        if ub:
            T[i, n + si] = 1.0
            si += 1
        T[i, -1] = b
        if b < 0:
            T[i, :] *= -1
        T[i, n + ns + i] = 1.0
        basis.append(n + ns + i)

    def pivot(r, col):
        T[r, :] /= T[r, col]
        for i in range(m):
            if i != r and abs(T[i, col]) > 0:
                T[i, :] -= T[i, col] * T[r, :]
        basis[r] = col

    def run(obj, allowed):
        for _ in range(5000):
            red = obj - obj[basis] @ T[:, :-1]
            cand = [j for j in range(N) if allowed[j] and red[j] > tol]
            if not cand:
                return True
            col = cand[0]
            rat = [(T[i, -1] / T[i, col], basis[i], i) for i in range(m) if T[i, col] > tol]
            if not rat:
                return False                      # unbounded
            mn = min(x[0] for x in rat)
            r = min((x for x in rat if x[0] <= mn + 1e-12), key=lambda x: x[1])[2]
            pivot(r, col)
        raise RuntimeError("cycling")

    obj1 = np.zeros(N)
    obj1[n + ns:] = -1.0
    run(obj1, [True] * N)
    if T[:, -1] @ (obj1[basis]) < -1e-9:
        return None                               # infeasible
    for i in range(m):                            # drive zero artificials out
        if basis[i] >= n + ns:
            nz = [j for j in range(n + ns) if abs(T[i, j]) > 1e-9]
            if nz:
                pivot(i, nz[0])
    obj2 = np.zeros(N)
    obj2[:n] = c
    allowed = [j < n + ns for j in range(N)]
    if not run(obj2, allowed):
        return np.inf
    x = np.zeros(N)
    for i, b in enumerate(basis):
        x[b] = T[i, -1]
    return float(obj2 @ x)


def k_lp(m, l, bilinear=False, b0_le_a0=True):
    best = -np.inf
    for col_branch, term in itertools.product((0, 1), range(3)):
        c = np.zeros(NV)
        c[V['O']] = 0.5
        c[V['H']] = -1.0 + (0 if col_branch else 1)
        c[V['v']] = 0.5 + (1 if col_branch else 0)
        c[V['lb']] = 0.5 + (1 if col_branch else 0)
        c[V['S0']] = -1.0
        c[V['B0']] = 0.5
        c[V['el']] = -1/6
        c[V['A0']] = -1.0
        c[V['za']] = -1.0
        c[V['N0']] = 0.5
        c[V['t']] = 1.0
        const = (1 + l) / 2 - m
        Aub, bub = [], []

        def row(**kw):
            a = np.zeros(NV)
            for k_, val in kw.items():
                a[V[k_]] = val
            return a
        Aub.append(row(A0=2, O=-1)); bub.append(0)
        Aub.append(row(N0=1, B0=1, S0=1, A0=-1) if b0_le_a0 else row(N0=1, A0=-1)); bub.append(0)
        Aub.append(row(za=1)); bub.append(l)
        Aub.append(row(v=1, lb=3, el=1, H=-2, A0=-2, za=-2, N0=1, B0=3)); bub.append(-(1 + l))
        if col_branch:                           # max(H, v+lb) = v+lb  (consistent branch)
            Aub.append(row(H=1, v=-1, lb=-1)); bub.append(0)
        else:
            Aub.append(row(v=1, lb=1, H=-1)); bub.append(0)
        if term == 0:
            if bilinear:                         # trilinear sieve: no L term, branch absent
                continue
            Aub.append(row(t=1, za=-1, v=1)); bub.append(0)
        elif term == 1:
            Aub.append(row(t=1)); bub.append(0)
        else:
            Aub.append(row(t=1, za=-2/3, v=1/3)); bub.append(0)
        val = simplex_max(c, Aub, bub, [row(O=1, H=1, dH=1)], [m])
        if val is not None:
            best = max(best, val + const)
    return max(best, 0.0)


def k_liu(m, l):
    return np.maximum.reduce([0*m + 0*l, (1 + 3*l - 2*m)/4, m + l - 1, (3*l - 1)/2])


def k_bil(m, l):                                 # closed form of the trilinear LP (checked below)
    return np.maximum.reduce([0*m + 0*l, (3 + 5*l - 6*m)/12, m + l - 1, (3*l - 1)/2])


def b_low(lx, ly, l, kfun):
    """Liu Prop. 22.2 direct boundary with all subsets d_I in [0,l] (grid)."""
    M, b = lx + ly, ly - lx
    d = np.linspace(0, l, 41)
    gram = np.maximum(np.maximum(0, b/12), b - ly/2 + d/2)
    return float((lx/2 - d + gram + kfun(M - 2*d, l - d)/2).max()) + 5/6 - lx/3 - l/6


def refine(f, x0, step=0.01, ok=lambda x: True):
    """Coordinate pattern search minimising f from x0."""
    x, fx = list(x0), f(x0)
    while step > 1e-6:
        moved = False
        for i in range(len(x)):
            for s in (step, -step):
                y = list(x)
                y[i] += s
                if ok(y):
                    fy = f(y)
                    if fy < fx - 1e-13:
                        x, fx, moved = y, fy, True
        if not moved:
            step /= 2
    return fx, tuple(round(t, 5) for t in x)


def floor(kfun):
    """min over (l, M) of the direct boundary (b = 0; b does not lower it), lx = ly = M/2 > l."""
    best = (9.0, None)
    for lv in np.arange(0.12, 0.38, 0.0025):
        for M in np.arange(0.5, 1.0, 0.0025):
            if M/2 > lv:
                B = b_low(M/2, M/2, lv, kfun)
                if B < best[0]:
                    best = (B, (round(lv, 4), round(M, 4)))
    return best


_dg = np.linspace(1/8, 5/6, 241)
_xg = np.linspace(0, 0.5, 26)
_D, _X = np.meshgrid(_dg, _xg, indexing="ij")


def B_needed_dh(lx, ly, l, kfun):
    """Smallest B with B >= B_low and max E <= 0 under a DH count R = 1-delta (Liu (14.7), d = h)."""
    a, q, R = (1 + _D)/2, _D*_X, 1 - _D
    h = 1 - lx + l
    B = b_low(lx, ly, l, kfun)
    for _ in range(4):                          # kappa = 2B-1 fixed point (enters theta only)
        c = 1/(3*(2*B - 1))
        Eav = a - ly/2 - h/6 - l/2 + l*q + h*(3*R + 1)/4 + (q*c/2)*(h/2 - ly)
        Ept = a*(1 - ly) - h/6 - l/2 + l*q + h*(R + _D/2)
        B = max(b_low(lx, ly, l, kfun), float(np.minimum(Eav, Ept).max()))
    return B


def system_floor(kfun):
    """min over free (lx, ly, l) of max(B_low, generic-row endpoint).  The endpoint is Liu (14.7) at
    d = h, delta -> 0 (a = 1/2, g = theta = 0, R = 1, K = 1): B >= 4/3 - ly/2 - 5lx/6 + l/3.  It uses no
    count or moment estimate (it holds even under GRH for every row L-function)."""
    def f(y):
        lv, lx, ly = y
        return max(b_low(lx, ly, lv, kfun), 4/3 - ly/2 - 5*lx/6 + lv/3)
    ok = lambda y: y[1] > y[0] > 0 and y[2] >= y[1]
    starts = [(lv, (1 - lv - b)/2 + e, (1 - lv + b)/2 + e) for lv in np.linspace(0.14, 0.32, 10)
              for b in (0.0, 0.05) for e in (-0.03, 0.0, 0.03)]
    best = min((f(s), s) for s in starts if ok(s))
    return refine(f, best[1], ok=ok)


def b_dh(kfun):
    best = (9.0, None)
    for lv in np.arange(0.15, 0.30, 0.01):
        for lx in np.arange(0.30, 0.46, 0.01):
            for ly in np.arange(lx, 0.52, 0.01):
                if lx > lv and ly > lv:
                    B = B_needed_dh(lx, ly, lv, kfun)
                    if B < best[0]:
                        best = (B, (round(lv, 3), round(lx, 3), round(ly, 3)))
    return best


if __name__ == "__main__":
    grid = [(m, l) for m in np.linspace(0.55, 1.1, 12) for l in np.linspace(0.0, 0.3, 7)]
    e1 = max(abs(k_lp(m, l) - k_liu(m, l)) for m, l in grid)
    e1b = max(abs(k_lp(m, l, b0_le_a0=False) - k_liu(m, l)) for m, l in grid)
    e2 = max(abs(k_lp(m, l, bilinear=True) - k_bil(m, l)) for m, l in grid)
    print(f"1) Lemma 14.3 LP vs Liu k(m,l) on {len(grid)} points: max |diff| = {e1:.1e} "
          f"(without B0<=A0: {e1b:.1e})  -> the (1+3l-2m)/4 excess is attained, not a bookkeeping loss")
    print(f"2) trilinear mark-mark LP vs max(0,(3+5l-6m)/12,m+l-1,(3l-1)/2): max |diff| = {e2:.1e}")
    f1, f2 = floor(k_liu), floor(k_bil)
    print(f"3) direct-side floor  min B_low: published {f1[0]:.6f} (13/15={13/15:.6f}) at (l,M)={f1[1]}")
    print(f"   with trilinear mark-mark step: {f2[0]:.6f} at (l,M)={f2[1]}; exact value at (3/11, 8/11): "
          f"{b_low(4/11, 4/11, 3/11, k_bil):.6f} = 28/33 = {28/33:.6f}")
    d1, d2 = b_dh(k_liu), b_dh(k_bil)
    print(f"4) hypothetical DH count R=1-delta: B = {d1[0]:.6f} (published floor), {d2[0]:.6f} (trilinear floor)"
          f"  at (l,lx,ly)={d2[1]}")
    s1, s2 = system_floor(k_liu), system_floor(k_bil)
    print(f"5) floor of the whole system with the count-free generic-row endpoint (delta->0 of (14.7)):"
          f" {s1[0]:.6f} (published k), {s2[0]:.6f} (trilinear k) at (l,lx,ly)={s2[1]}")
    print("   balanced lengths: endpoint = 2/3 + l + b/6 = C0 + B(l), so C0<0 forces l < 1/5 - 2b/15 and"
          " B(l) = 11/12 - l/4 > 13/15: the 13/15 floor is two-sided.")
