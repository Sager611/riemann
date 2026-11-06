#!/usr/bin/env python3
"""lincomb3 numerics step 2 (adversarial test of the dominance mechanism), q = 7, three odd primitive characters.
Input: Z_q7_n{3,5,6}.npy from zgrid.py (Arb Hardy Z on t = 5 + 0.03 k, t <= 2400).
G(t) = sum_j b_j Z_j(t) is e^{i vartheta(t)} F(1/2+it) for F = sum_j b_j' e^{i theta_j} L(s,chi_j) (b_j' = +-b_j).
Prints: on-line zeros of F (sign changes of G) against the main term N(T) = (T/2pi) log(qT/(2 pi e)) + 7/8-type O(1),
block dominance statistics, how often a zero of F sits next to a zero of the dominant L_j, and an empirical test of the
missing lemma (are non-isolated zeros of L_j over-represented in the blocks where L_j dominates?).
Run: uv run --no-project --with numpy python3 analyse.py
"""
import numpy as np, math, sys
q, T0, H = 7, 5.0, 0.03
ns = (3, 5, 6)
Z = [np.load('Z_q7_n%d.npy' % n) for n in ns]
K = min(len(z) for z in Z); Z = np.array([z[:K] for z in Z])
t = T0 + H * np.arange(K)

def sign_changes(v):
    s = np.sign(v)
    idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    # linear interpolation of the crossing
    return t[idx] - v[idx] * (t[idx + 1] - t[idx]) / (v[idx + 1] - v[idx])

def Nmain(T):  # zeros of L(s,chi) (or F) with 0 < gamma <= T: (T/2pi) log(qT/(2 pi e)) (+ O(log T))
    return T / (2 * math.pi) * math.log(q * T / (2 * math.pi * math.e))

zj = [sign_changes(Z[j]) for j in range(3)]
print('grid: t in [%.1f, %.1f], step %.3f, %d points' % (t[0], t[-1], H, K))
for j, n in enumerate(ns):
    for (a, b) in ((T0, 600), (600, 1200), (1200, 2400)):
        c = np.sum((zj[j] > a) & (zj[j] <= b))
        print('L(s,chi_%d): sign changes in (%5.0f,%5.0f]: %5d   main term %8.1f   ratio %.4f' % (n, a, b, c, Nmain(b) - Nmain(a), c / (Nmain(b) - Nmain(a))))

def isolated_frac(z, a, b, kappa):
    """fraction of sign changes in (a,b] whose two neighbouring sign changes are > kappa mean spacings away (normalised)."""
    sp = 2 * math.pi / np.log(q * z / (2 * math.pi))      # local mean spacing
    g = np.diff(z)
    left = np.r_[np.inf, g]; right = np.r_[g, np.inf]
    iso = (left > kappa * sp) & (right > kappa * sp)
    m = (z > a) & (z <= b)
    return iso[m].mean(), iso

for kappa in (0.05, 0.2, 0.35):
    print('kappa = %.2f: fraction of on-line zeros of L_j isolated from other on-line zeros, t in (1200,2400]: %s'
          % (kappa, ', '.join('%.4f' % isolated_frac(zj[j], 1200, 2400, kappa)[0] for j in range(3))))

M = 10.0   # block length M / log t (in t), u-length M
for b in ((1.0, 1.0, 1.0), (1.0, -2.0, 0.5), (3.0, 1.0, -1.0)):
    bv = np.array(b)
    G = (bv[:, None] * Z).sum(0)
    zg = sign_changes(G)
    print('\n== b = %s' % (b,))
    for (a, bb) in ((T0, 600), (600, 1200), (1200, 2400), (T0, 2400)):
        c = np.sum((zg > a) & (zg <= bb))
        print('F: sign changes of G in (%5.0f,%5.0f]: %5d   main term %8.1f   proportion on line ~ %.4f' % (a, bb, c, Nmain(bb) - Nmain(a), c / (Nmain(bb) - Nmain(a))))
    # blocks
    logA = np.log(np.abs(bv[:, None] * Z) + 1e-300)
    edges = [T0]
    while edges[-1] < t[-1]:
        edges.append(edges[-1] + M / math.log(edges[-1]))
    edges = np.array(edges[:-1])
    bi = np.searchsorted(edges, t, side='right') - 1
    nb = len(edges) - 1
    med = np.full((3, nb), -np.inf)
    for k in range(nb):
        sl = slice(np.searchsorted(t, edges[k]), np.searchsorted(t, edges[k + 1]))
        if sl.stop - sl.start < 5: continue
        med[:, k] = np.median(logA[:, sl], axis=1)
    order = np.argsort(-med, axis=0)
    dom = order[0]; gap = med[order[0], np.arange(nb)] - med[order[1], np.arange(nb)]
    hi = edges[1:] > 1200
    print('blocks (u-length %.0f) in (1200,2400]: %d; dominant index shares %s; median top-two gap %.2f; gap > 1 in %.2f of blocks'
          % (M, hi.sum(), [round(float(np.mean(dom[hi] == j)), 3) for j in range(3)], float(np.median(gap[hi])), float(np.mean(gap[hi] > 1))))
    # tracking: zeros of F next to a zero of the dominant L_j
    sp = 2 * math.pi / np.log(q * zg / (2 * math.pi))
    kb = np.clip(np.searchsorted(edges, zg, side='right') - 1, 0, nb - 1)
    near = np.zeros(len(zg), bool)
    for j in range(3):
        m = dom[kb] == j
        if m.sum() == 0: continue
        pos = np.searchsorted(zj[j], zg[m])
        d1 = np.abs(zj[j][np.clip(pos, 0, len(zj[j]) - 1)] - zg[m]); d0 = np.abs(zj[j][np.clip(pos - 1, 0, len(zj[j]) - 1)] - zg[m])
        near[m] = np.minimum(d0, d1) < 0.15 * sp[m]
    for thr in (0.0, 1.0, 2.0):
        m = (zg > 1200) & (gap[kb] > thr)
        print('zeros of F in (1200,2400] in blocks with top-two gap > %.0f: %5d, of which within 0.15 mean spacings of a zero of the dominant L_j: %.3f'
              % (thr, m.sum(), near[m].mean() if m.sum() else float('nan')))
    # per-block mechanism: #zeros of F in block >= #(isolated zeros of dominant L_j in block) - 2 ?
    fails = 0; tot = 0
    for k in np.nonzero(hi & (gap > 2.0))[0]:
        a, bb = edges[k], edges[k + 1]
        j = dom[k]
        _, iso = isolated_frac(zj[j], 0, 1e9, 0.35)
        cj = np.sum(iso & (zj[j] > a) & (zj[j] <= bb))
        cg = np.sum((zg > a) & (zg <= bb))
        tot += 1; fails += cg < cj - 2
    print('blocks in (1200,2400] with top-two gap > 2: %d; per-block count < (0.35-isolated zeros of dominant L_j) - 2 in %d' % (tot, fails))
    # empirical test of the missing lemma: share of non-isolated zeros of L_j lying in L_j's dominance blocks
    for kappa in (0.35, 0.5):
        out = []
        for j in range(3):
            _, iso = isolated_frac(zj[j], 0, 1e9, kappa)
            kz = np.clip(np.searchsorted(edges, zj[j], side='right') - 1, 0, nb - 1)
            m = zj[j] > 1200
            bad = m & ~iso; good = m & iso
            th = np.mean(dom[kz[bad]] == j) if bad.sum() else float('nan')
            sh = np.mean(dom[kz[m]] == j)
            out.append('j=%d: |H_j|=%d, share in own dominance blocks %.3f (share of all its zeros %.3f)' % (ns[j], bad.sum(), th, sh))
        print('missing-lemma test, kappa = %.2f: ' % kappa + '; '.join(out))
