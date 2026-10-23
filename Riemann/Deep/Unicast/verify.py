#!/usr/bin/env python3
"""Audit of arXiv 2610.10108: claimed counterexample to the undirected multiple-unicast conjecture.

Exact arithmetic over F9.  Reads instance.json (R, pi, C, J copied verbatim from the e-print).
  python3 verify.py          full check under the paper's own conventions (Section 3)
  python3 verify.py --alt    degree check under alternative F9 encodings / orderings
  python3 verify.py --control  negative control: wrong line rule, the decoder checks must fail
"""
import json, os, sys
from collections import deque
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
INST = json.load(open(os.path.join(HERE, "instance.json")))


class Field:
    """F9 = F3[t]/(t^2 + c1 t + c0); element a + b t has index a + 3b."""

    def __init__(self, c1, c0):
        add = np.zeros((9, 9), np.uint8)
        mul = np.zeros((9, 9), np.uint8)
        for x in range(9):
            a1, b1 = x % 3, x // 3
            for y in range(9):
                a2, b2 = y % 3, y // 3
                add[x, y] = (a1 + a2) % 3 + 3 * ((b1 + b2) % 3)
                bb = b1 * b2  # t^2 = -c1 t - c0
                a = (a1 * a2 - c0 * bb) % 3
                b = (a1 * b2 + a2 * b1 - c1 * bb) % 3
                mul[x, y] = a + 3 * b
        self.add, self.mul = add, mul
        self.neg = np.array([(-(x % 3)) % 3 + 3 * ((-(x // 3)) % 3) for x in range(9)], np.uint8)
        inv = np.zeros(9, np.uint8)
        for x in range(1, 9):
            ys = [y for y in range(1, 9) if mul[x, y] == 1]
            assert len(ys) == 1, "modulus not irreducible"
            inv[x] = ys[0]
        self.inv = inv

    def pw(self, x, n):
        r = 1
        for _ in range(n):
            r = int(self.mul[r, x])
        return r

    def dot(self, x, y):
        s = 0
        for a, b in zip(x, y):
            s = int(self.add[s, self.mul[a, b]])
        return s

    def rref(self, M):
        M = M.copy()
        rows, cols = M.shape
        piv, r = [], 0
        for c in range(cols):
            if r == rows:
                break
            nz = np.nonzero(M[r:, c])[0]
            if len(nz) == 0:
                continue
            p = r + nz[0]
            if p != r:
                M[[r, p]] = M[[p, r]]
            M[r] = self.mul[self.inv[M[r, c]], M[r]]
            f = M[:, c].copy()
            f[r] = 0
            idx = np.nonzero(f)[0]
            if len(idx):
                M[idx] = self.add[M[idx], self.neg[self.mul[f[idx][:, None], M[r][None, :]]]]
            piv.append(c)
            r += 1
        return M, piv

    def nullspace(self, M):
        Rm, piv = self.rref(M)
        cols = M.shape[1]
        free = [c for c in range(cols) if c not in piv]
        basis = []
        for fc in free:
            v = np.zeros(cols, np.uint8)
            v[fc] = 1
            for i, pc in enumerate(piv):
                v[pc] = self.neg[Rm[i, fc]]
            basis.append(v)
        return np.array(basis, np.uint8).reshape(len(basis), cols)

    def matmul(self, A, B):
        out = np.zeros((A.shape[0], B.shape[1]), np.uint8)
        for k in range(A.shape[1]):
            out = self.add[out, self.mul[A[:, k][:, None], B[k][None, :]]]
        return out

    def inverse(self, A):
        n = A.shape[0]
        Rm, piv = self.rref(np.hstack([A, np.eye(n, dtype=np.uint8)]))
        assert piv == list(range(n)), "singular interpolation matrix"
        return Rm[:, n:]


def paper_eta(r):  # eta_r = (r mod 3) + floor(r/3) theta  ->  index r
    return (r % 3) + 3 * (r // 3)


def build(F, eta):
    T = [(0, 0, 1)] + [(0, 1, eta(r)) for r in range(9)] + \
        [(1, eta(r), eta(s)) for r in range(9) for s in range(9)]
    assert len(set(T)) == 91
    nb = [sorted(j for j in range(91) if F.dot(T[i], T[j]) == 0) for i in range(91)]
    assert all(len(x) == 10 for x in nb)
    edges = [(i, 91 + nb[i][r]) for i in range(91) for r in range(10)]
    inc = [[] for _ in range(182)]
    for q, (a, b) in enumerate(edges):
        inc[a].append(q)
        inc[b].append(q)
    return T, edges, inc


def local_generators(F, T, edges, inc, cube=True):
    """5 x 10 generator of the local code at every vertex (columns ordered as inc[v])."""
    Gs = []
    for v in range(182):
        z = T[v] if v < 91 else T[v - 91]
        zp = [t for t in T if F.dot(z, t) == 0]
        u, w = zp[0], zp[1]
        G = np.zeros((5, 10), np.uint8)
        for k, q in enumerate(inc[v]):
            i, l = edges[q]
            y = T[l - 91] if v < 91 else T[i]
            sol = [(la, mu) for la in range(9) for mu in range(9)
                   if tuple(int(F.add[F.mul[la, u[c]], F.mul[mu, w[c]]]) for c in range(3)) == y]
            assert len(sol) == 1
            la, mu = sol[0]
            if v >= 91 and cube:  # line vertex: Q(lambda^3, mu^3)
                la, mu = F.pw(la, 3), F.pw(mu, 3)
            for r in range(5):
                G[r, k] = F.mul[F.pw(la, 4 - r), F.pw(mu, r)]
        Gs.append(G)
    return Gs


def orientation_check(edges, inc, R, pi, C, J):
    pos = {v: k for k, v in enumerate(pi)}
    Rs, Cs, Js = set(R), set(C), set(J)
    bad, cj = 0, 0
    for v in range(182):
        din = dout = cout = jin = 0
        for q in inc[v]:
            a, b = edges[q]
            other = b if a == v else a
            later = pos[v] > pos[other]
            incoming = later if q not in Rs else not later
            if incoming:
                din += 1
                jin += q in Js
            else:
                dout += 1
                cout += q in Cs
        bad += not (din == 5 and dout == 5)
        cj += cout != jin
    return bad, cj


def main_paper(control=False):
    F = Field(0, 1)
    T, edges, inc = build(F, paper_eta)
    R, pi, C, J = INST["R"], INST["pi"], INST["C"], INST["J"]
    Rs, Cs, Js = set(R), set(C), set(J)
    assert len(Rs) == 175 and Cs <= Rs and Js <= Rs and not (Cs & Js) and sorted(pi) == list(range(182))
    B = [q for q in R if q not in Cs and q not in Js]
    phys = [q for q in range(910) if q not in Rs]
    print(f"[graph] |V|=182 |E(H)|=910 |R|={len(Rs)} |E(G)|={len(phys)} |B|={len(B)}")
    bad, cj = orientation_check(edges, inc, R, pi, C, J)
    print(f"[orient] vertices violating d-=d+=5: {bad}; violating dC+=dJ-: {cj}")

    # distances in G
    adj = [[] for _ in range(182)]
    for q in phys:
        a, b = edges[q]
        adj[a].append(b)
        adj[b].append(a)
    pos = {v: k for k, v in enumerate(pi)}
    st = {}
    for q in B:
        a, b = edges[q]
        st[q] = (a, b) if pos[a] < pos[b] else (b, a)
    dists = []
    for q in B:
        s, t = st[q]
        d = {s: 0}
        dq = deque([s])
        while dq:
            x = dq.popleft()
            for y in adj[x]:
                if y not in d:
                    d[y] = d[x] + 1
                    dq.append(y)
        dists.append(d.get(t, -1))
    from collections import Counter
    print(f"[dist] session distance histogram {dict(Counter(dists))}; sum={sum(dists)}; "
          f"flow bound 735/{sum(dists)} = {735/sum(dists):.6f}")

    # local codes and W
    Gs = local_generators(F, T, edges, inc, cube=not control)
    if control:
        print('[control] NEGATIVE CONTROL: line vertices use Q(lambda, mu) instead of Q(lambda^3, mu^3)')
    assert all(len(F.rref(G)[1]) == 5 for G in Gs)
    rows = []
    for v in range(182):
        for c in F.nullspace(Gs[v]):
            row = np.zeros(910, np.uint8)
            row[inc[v]] = c
            rows.append(row)
    Hg = np.array(rows, np.uint8)
    _, piv = F.rref(Hg)
    dimW = 910 - len(piv)
    print(f"[W] parity checks {Hg.shape}, rank {len(piv)}, dim W = {dimW}")

    # appendix: the 166 explicit assignments h(p^(3), l)
    mons = []
    for a0 in range(5):
        for a1 in range(5 - a0):
            for b0 in range(5):
                for b1 in range(5 - b0):
                    a = (a0, a1, 4 - a0 - a1)
                    b = (b0, b1, 4 - b0 - b1)
                    if (a[0] >= 1 and b[0] >= 3) or (a[1] >= 3 and b[1] >= 1):
                        continue
                    mons.append((a, b))
    Wm = np.zeros((len(mons), 910), np.uint8)
    for m, (a, b) in enumerate(mons):
        for q, (i, l) in enumerate(edges):
            p, L = T[i], T[l - 91]
            val = 1
            for c in range(3):
                val = int(F.mul[val, F.pw(F.pw(p[c], 3), a[c])])
                val = int(F.mul[val, F.pw(L[c], b[c])])
            Wm[m, q] = val
    inW = not F.matmul(Hg, Wm.T).any()
    print(f"[appendix] {len(mons)} monomials; all in W: {inW}; rank {len(F.rref(Wm)[1])}")

    # causal simulation with symbolic inputs: value = vector in F9^157
    nB = len(B)
    bidx = {q: k for k, q in enumerate(B)}
    unit = np.eye(nB, dtype=np.uint8)
    zero = np.zeros(nB, np.uint8)
    sent = {}
    computed = {}  # (q, v) -> value computed at v
    worst = 0
    known_at = {}
    for v in pi:
        known = []
        for q in inc[v]:
            a, b = edges[q]
            other = b if a == v else a
            earlier_other = pos[other] < pos[v]
            if q not in Rs and earlier_other:
                known.append((q, sent[q]))
            elif q in bidx and not earlier_other:  # v = s_q
                known.append((q, unit[bidx[q]]))
            elif q in Cs:
                known.append((q, zero))
        worst = max(worst, abs(len(known) - 5))
        assert len(known) == 5, (v, len(known))
        known_at[v] = known
        cols = [inc[v].index(q) for q, _ in known]
        A = Gs[v][:, cols].T
        coef = F.matmul(F.inverse(A), np.array([x for _, x in known], np.uint8))
        allv = F.matmul(Gs[v].T, coef)
        for k, q in enumerate(inc[v]):
            computed[(q, v)] = allv[k]
            a, b = edges[q]
            other = b if a == v else a
            if q not in Rs and pos[other] > pos[v]:
                sent[q] = allv[k]
    n_bad = sum(not np.array_equal(computed[(q, st[q][1])], unit[bidx[q]]) for q in B)
    dec_ok = n_bad == 0
    j_ok = all(np.array_equal(computed[(q, edges[q][0])], computed[(q, edges[q][1])]) for q in J)
    c_ok = all(not computed[(q, v)].any() for q in C for v in edges[q])
    w = np.zeros((910, nB), np.uint8)
    mism = 0
    for q in range(910):
        a, b = edges[q]
        mism += not np.array_equal(computed[(q, a)], computed[(q, b)])
        w[q] = computed[(q, a)]
    in_W = not F.matmul(Hg, w).any()
    rk = len(F.rref(w.T.copy())[1])
    print(f"[code] every vertex has exactly 5 known values: {worst == 0}; all 157 sinks decode their own message: {dec_ok} (failures {n_bad}); edges whose two ends disagree: {mism}")
    print(f"[code] J-edge values agree at both ends: {j_ok}; C-edges stay 0: {c_ok}; "
          f"global assignment in W: {in_W}; rank of input->assignment map: {rk}")
    # model-level decodability: e_q lies in the span of what t_q holds (incoming physical symbols + own messages)
    span_ok = True
    for q in B:
        A = np.array([x for _, x in known_at[st[q][1]]], np.uint8)
        r1 = len(F.rref(A)[1])
        r2 = len(F.rref(np.vstack([A, unit[bidx[q]][None, :]]))[1])
        span_ok &= r1 == r2
    print(f"[code] each x_q is in the span of the symbols available at t_q: {span_ok}")
    n_src_pts = sum(st[q][0] < 91 for q in B)
    print(f"[code] sources on point side {n_src_pts}, on line side {len(B) - n_src_pts}; "
          f"max sessions per source vertex {max(Counter(st[q][0] for q in B).values())}")
    print(f"[flow] common rate <= 735/785 = 147/157; total throughput <= 735/5 = 147 < 157; "
          f"with distinct pendant terminals (unit edges) <= {735 + 2 * nB}/{nB * 7} = {(735 + 2 * nB) / (7 * nB):.4f}")
    print(f"[code] symbols per physical edge per execution: 1 (one direction, earlier->later); edges used {len(sent)}")


def main_alt():
    R, pi, C, J = INST["R"], INST["pi"], INST["C"], INST["J"]
    polys = {"t^2+1": (0, 1), "t^2+t+2": (1, 2), "t^2+2t+2": (2, 2)}
    etas = {"paper eta_r=(r%3)+(r//3)t": paper_eta,
            "swapped eta_r=(r//3)+(r%3)t": lambda r: (r // 3) + 3 * (r % 3)}
    for pn, (c1, c0) in polys.items():
        F = Field(c1, c0)
        for en, eta in etas.items():
            T, edges, inc = build(F, eta)
            bad, cj = orientation_check(edges, inc, R, pi, C, J)
            print(f"[alt] {pn:9s} {en:30s} vertices violating d-=d+=5: {bad:3d}; dC+!=dJ-: {cj:3d}")


if __name__ == "__main__":
    main_alt() if "--alt" in sys.argv else main_paper(control="--control" in sys.argv)

