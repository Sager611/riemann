#!/usr/bin/env python3
"""verify.py -- rigorous verifier for the certified lower bound V*(K_4) >= <value> (3D snaky-corridor sofa).

Usage: python3 verify.py <path.json> <certificate.npz>   (numpy and a C compiler 'cc'; verify_core.c next to this file)

Certified object: a body B = union of convex polytopes ("pieces") and a motion given by the pose table of the
certificate (translations T_k, quaternions Q_k); between consecutive poses: linear translation and uniform rotation
about the fixed axis of R(Q_k+1) R(Q_k)^T, where R(Q) is the exact rotation matrix of Q/|Q| (rational in Q).
Every inequality is rigorous (NOTE.md, section "Rounding"):
 - exact rational arithmetic (fractions.Fraction): rotation matrices and their float roundings, rotation angles,
   sag bounds, the core's error bound, the final volume;
 - outward-rounded interval arithmetic (numpy nextafter after every operation): piece corners and vertices, |q|
   bounds, start/end containment, piece validity and disjointness, connectivity;
 - a written worst-case error bound for the floating-point core verify_core.c: its inputs are bounded here and
   its total rounding error E (computed here exactly) must stay below ETA/1000, ETA = its safety margin.
"""
import json, math, os, subprocess, sys, tempfile, time
from fractions import Fraction as Fr
import numpy as np

ETA = 0.5e-9; BS = 128; FLIM = 4.0; C12 = 1 + 1e-12
U = Fr(1, 2 ** 53); TAU = Fr(1, 2 ** 1074)
SIGNS = [(-1, -1), (-1, 1), (1, -1), (1, 1)]

def dn(x): return np.nextafter(x, -np.inf)
def up(x): return np.nextafter(x, np.inf)
def fup(r):
    x = float(r)
    return x if Fr(x) >= r else float(np.nextafter(x, np.inf))

def smooth(path, tt, tr):
    """the pipeline's smooth_path_to_tolerances (float64)."""
    out = [path[0].tolist()]
    for i in range(len(path) - 1):
        p1 = path[i]; p2 = path[i + 1]
        td = np.linalg.norm(p2[:3] - p1[:3]); rd = np.max(np.abs(p2[3:] - p1[3:]))
        nt = int(np.ceil(td / tt)) if td > tt else 0
        nr = int(np.ceil(rd / tr)) if rd > tr else 0
        n = max(1, nt, nr)
        for j in range(1, n):
            out.append((p1 + (j / n) * (p2 - p1)).tolist())
        out.append(list(p2))
    return np.array(out)

def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2)

def euler_quat(yaw, pitch, roll):
    """quaternion of Rx(roll) Ry(pitch) Rz(yaw) (degrees; scipy 'zyx' convention)."""
    h = lambda d: math.radians(d) / 2
    return qmul((math.cos(h(roll)), math.sin(h(roll)), 0.0, 0.0),
                qmul((math.cos(h(pitch)), 0.0, math.sin(h(pitch)), 0.0), (math.cos(h(yaw)), 0.0, 0.0, math.sin(h(yaw)))))

def rot_exact(q):
    """exact rotation matrix (row-major, 9 Fractions) of the nonzero float quaternion q = (w,x,y,z)."""
    w, x, y, z = (Fr(float(v)) for v in q)
    ww, xx, yy, zz = w * w, x * x, y * y, z * z
    n = ww + xx + yy + zz
    xy, xz, yz, wx, wy, wz = x * y, x * z, y * z, w * x, w * y, w * z
    return [(ww + xx - yy - zz) / n, 2 * (xy - wz) / n, 2 * (xz + wy) / n,
            2 * (xy + wz) / n, (ww - xx + yy - zz) / n, 2 * (yz - wx) / n,
            2 * (xz - wy) / n, 2 * (yz + wx) / n, (ww - xx - yy + zz) / n]

def iscale(c, lo, hi):
    p = c * lo; q = c * hi
    return dn(np.minimum(p, q)), up(np.maximum(p, q))

def imul(alo, ahi, blo, bhi):
    p = [alo * blo, alo * bhi, ahi * blo, ahi * bhi]
    return dn(np.minimum(np.minimum(p[0], p[1]), np.minimum(p[2], p[3]))), up(np.maximum(np.maximum(p[0], p[1]), np.maximum(p[2], p[3])))

def ilin(c0, ca, cb, ulo, uhi, vlo, vhi):
    """interval of c0 + ca*u + cb*v (c0, ca, cb floats; u, v intervals)."""
    t1 = iscale(ca, ulo, uhi); t2 = iscale(cb, vlo, vhi)
    return dn(dn(c0 + t1[0]) + t2[0]), up(up(c0 + t1[1]) + t2[1])

def main():
    t0 = time.time()
    d = json.load(open(sys.argv[1])); poses = np.array(d['poses'] if isinstance(d, dict) else d, float)
    cert = np.load(sys.argv[2])
    planes = cert['planes'].astype(float); cells = cert['cells'].astype(np.int64)
    T = cert['T'].astype(float); Q = cert['Q'].astype(float)
    nc0, L, off, tt, tr = cert['meta']; nc0 = int(nc0); L = int(L); off = float(off)
    n = len(planes); m = len(T) - 1
    # ---- 1. motion: pose table = pipeline smoothing of the path (translations bit-exact, quaternions ~ Euler angles)
    assert poses[0].tolist() == [-3, 0, 0, 0, 0, 0] and poses[-1].tolist() == [0, 6, 8, 90, 90, 90]
    sm = smooth(poses, tt, tr)
    assert sm.shape[0] == m + 1 and np.array_equal(sm[:, :3], T), 'translations differ from the smoothed path'
    Qe = np.array([euler_quat(*p[3:]) for p in sm])
    sgn = np.sign((Qe * Q).sum(1)); assert (sgn != 0).all()
    qdev = float(np.abs(Q - sgn[:, None] * Qe).max()); assert qdev < 1e-12, qdev
    assert np.isfinite(T).all() and np.isfinite(Q).all()
    # ---- 2. exact rotations, their float roundings and the rotation angle of every step (Fractions)
    Rf = np.zeros((m + 1, 9)); dR = Fr(0); ratio = []; Qfr = [[Fr(float(v)) for v in q] for q in Q]
    for k in range(m + 1):
        Rx = rot_exact(Q[k])
        for i, e in enumerate(Rx):
            Rf[k, i] = float(e); dR = max(dR, abs(Fr(Rf[k, i]) - e))
        if k < m:
            a, b = Qfr[k], Qfr[k + 1]
            C = sum(a[i] * b[i] for i in range(4)); assert C > 0, 'quaternion signs not aligned'
            S = sum((a[i] * b[j] - a[j] * b[i]) ** 2 for i in range(4) for j in range(i + 1, 4))
            ratio.append(S / (C * C))                 # tan^2(theta/2) of step k
    assert max(ratio) < Fr(1, 10 ** 4)
    dRf = fup(dR); t1 = time.time()
    # ---- 3. pieces: corners, vertices, |q| bounds by outward-rounded intervals
    lev, J, K = cells[:, 0], cells[:, 1], cells[:, 2]
    assert (lev >= 0).all() and (lev <= L).all() and (J >= 0).all() and (K >= 0).all()
    w = (np.int64(1) << (L - lev)); wf = w.astype(float); N = nc0 << L; Nf = float(N)
    assert ((J + 1) * w <= N + 2).all() and ((K + 1) * w <= N + 2).all()
    hf = wf / (2 * Nf); hlo, hhi = dn(hf), up(hf)
    def pos(num):
        s = num + off
        return s / Nf, dn(dn(s) / Nf), up(up(s) / Nf)
    A0, aA, bA, B0, aB, bB = (planes[:, i] for i in range(6))
    # geometry: distinct dyadic cells (none inside another); pieces of one cell, ordered by A0, pairwise disjoint
    key = lambda l, j, k: (l << 50) | (j << 25) | k
    ukey, inv = np.unique(key(lev, J, K), return_inverse=True); uset = set(ukey.tolist())
    ul = ukey >> 50; uj = (ukey >> 25) & ((1 << 25) - 1); uk = ukey & ((1 << 25) - 1)
    for a_ in range(1, L + 1):
        s_ = ul >= a_
        assert not any(x in uset for x in key(ul[s_] - a_, uj[s_] >> a_, uk[s_] >> a_).tolist()), 'nested cells'
    del uset, ukey, ul, uj, uk
    order = np.lexsort((A0, inv)); same = inv[order][1:] == inv[order][:-1]
    i1, i2 = order[:-1][same], order[1:][same]; del order, same, inv
    sep = np.ones(len(i1), bool)
    # start pose (entry arm: x <= 0) and end pose (exit arm: z >= 1) at all 8 vertices, exact rotations as intervals
    Rs = [(dn(Rf[0, 0:3] - dRf), up(Rf[0, 0:3] + dRf)), (dn(Rf[m, 6:9] - dRf), up(Rf[m, 6:9] + dRf))]
    def world(Rl, Rh, Tc, x, y, z):
        p0 = imul(Rl[0], Rh[0], *x); p1 = imul(Rl[1], Rh[1], *y); p2 = imul(Rl[2], Rh[2], *z)
        return dn(dn(dn(p0[0] + p1[0]) + p2[0]) + Tc), up(up(up(p0[1] + p1[1]) + p2[1]) + Tc)
    start_ok = np.ones(n, bool); end_ok = np.ones(n, bool); valid = np.ones(n, bool)
    ea = ed = ey = ez = 0.0; X = D = Y = Z = QR2 = 0.0
    tdir = tempfile.mkdtemp(); fpc = os.path.join(tdir, 'pieces.bin'); pf = open(fpc, 'wb')
    for c, (sy, sz) in enumerate(SIGNS):
        u = (-hhi, -hlo) if sy < 0 else (hlo, hhi); v = (-hhi, -hlo) if sz < 0 else (hlo, hhi)
        Alo, Ahi = ilin(A0, aA, bA, u[0], u[1], v[0], v[1]); Blo, Bhi = ilin(B0, aB, bB, u[0], u[1], v[0], v[1])
        af = A0 + aA * (sy * hf) + bA * (sz * hf); bf = B0 + aB * (sy * hf) + bB * (sz * hf); df = bf - af; del bf
        dlo, dhi = dn(Blo - Ahi), up(Bhi - Alo)
        yf, ylo, yhi = pos(((J + (sy > 0)) * w).astype(float)); zf, zlo, zhi = pos(((K + (sz > 0)) * w).astype(float))
        valid &= dlo > 0
        sep &= Bhi[i1] <= Alo[i2]
        ea = max(ea, float(np.maximum(up(af - Alo), up(Ahi - af)).max())); ed = max(ed, float(np.maximum(up(df - dlo), up(dhi - df)).max()))
        ey = max(ey, float(np.maximum(up(yf - ylo), up(yhi - yf)).max())); ez = max(ez, float(np.maximum(up(zf - zlo), up(zhi - zf)).max()))
        X = max(X, float(np.abs(Alo).max()), float(np.abs(Ahi).max()), float(np.abs(Blo).max()), float(np.abs(Bhi).max()))
        D = max(D, float(np.abs(dlo).max()), float(np.abs(dhi).max()))
        Y = max(Y, float(np.abs(ylo).max()), float(np.abs(yhi).max())); Z = max(Z, float(np.abs(zlo).max()), float(np.abs(zhi).max()))
        ya = np.maximum(np.abs(ylo), np.abs(yhi)); za = np.maximum(np.abs(zlo), np.abs(zhi))
        for lo_, hi_ in ((Alo, Ahi), (Blo, Bhi)):
            xa = np.maximum(np.abs(lo_), np.abs(hi_))
            QR2 = max(QR2, float(up(up(up(xa * xa) + up(ya * ya)) + up(za * za)).max()))
            start_ok &= world(*Rs[0], T[0, 0], (lo_, hi_), (ylo, yhi), (zlo, zhi))[1] <= 0
            end_ok &= world(*Rs[1], T[m, 2], (lo_, hi_), (ylo, yhi), (zlo, zhi))[0] >= 1
        for arr in (af, df, yf, zf): pf.write(np.ascontiguousarray(arr, float).tobytes())
        del Alo, Ahi, Blo, Bhi, af, df, dlo, dhi, yf, ylo, yhi, zf, zlo, zhi, ya, za
    pf.close()
    assert sep.all(), 'pieces in one cell not provably disjoint: %d' % int((~sep).sum())
    QR = float(up(np.sqrt(QR2)))
    # ---- 4. sag bound per step: eps_k >= QR theta_k^2/8, using theta <= 2 tan(theta/2)  (Fractions, rounded up)
    eps = np.array([fup(Fr(QR) * r / 2) for r in ratio]); epsmax = float(eps.max())
    # ---- 5. worst-case rounding error of the core (NOTE.md "Rounding"), exact rational evaluation
    Rm = float(np.abs(Rf).max()); Tm = float(np.abs(T).max())
    F_, X_, D_, Y_, Z_, R_, d_, T_ = (Fr(v) for v in (FLIM, X, D, Y, Z, Rm, dRf, Tm))
    ea_, ed_, ey_, ez_ = (Fr(v) for v in (ea, ed, ey, ez)); u, t = U, TAU; g4 = 4 * u / (1 - 4 * u)
    Ea = d_ * (X_ + Y_ + Z_ + ea_ + ey_ + ez_) + (R_ + d_) * (ea_ + ey_ + ez_) + g4 * (R_ * (X_ + Y_ + Z_ + ea_ + ey_ + ez_) + T_) + 3 * t
    Eb = d_ * (D_ + ed_) + (R_ + d_) * ed_ + u * R_ * (D_ + ed_) + t
    AL = (R_ * (X_ + Y_ + Z_ + ea_ + ey_ + ez_) + T_) * (1 + g4) + 3 * t
    BE = R_ * (D_ + ed_) * (1 + u) + t
    Ediv = ((6 + F_) * (1 + u) + AL) * (2 * u + u * u) + BE * t
    EH = u * (6 + F_)
    E_core = u * (Fr(epsmax) + Fr(ETA)) + EH + Ediv + Ea + Eb
    E_block = EH + Ediv + Ea + Eb + 2 * d_ * (X_ + Y_ + Z_) + 5 * t + u * Fr(ETA)
    assert Fr(C12) - 1 >= Fr(1, 10 ** 13) and Fr(epsmax) + Fr(ETA) <= F_
    assert E_core <= Fr(ETA) / 1000 and E_block <= Fr(ETA) / 1000, (float(E_core), float(E_block))
    # ---- 6. run the core
    C = np.zeros((m + 1, 26))
    mp = lambda Rr, Tt: np.concatenate([Rr[:, [0, 3, 6]], Rr[:, [1, 4, 7]], Rr[:, [2, 5, 8]], Tt], 1)
    C[0, 1:13] = mp(Rf[:1], T[:1])[0]
    C[1:, 0] = 1; C[1:, 1:13] = mp(Rf[:-1], T[:-1]); C[1:, 13:25] = mp(Rf[1:], T[1:]); C[1:, 25] = eps
    td = tdir; exe = os.path.join(td, 'vcore'); fi = os.path.join(td, 'in.bin'); fo = os.path.join(td, 'out.bin')
    subprocess.run(['cc', '-O3', '-ffp-contract=off', '-o', exe, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'verify_core.c')], check=True)
    with open(fi, 'wb') as f:
        f.write(np.array([m + 1, n, ETA, X, Y, Z, BS, FLIM, C12], float).tobytes()); f.write(C.tobytes())
    tc = time.time(); r = subprocess.run([exe, fi, fpc, fo], capture_output=True, text=True); tc = time.time() - tc
    assert r.returncode == 0, r.stderr
    full = np.frombuffer(open(fo, 'rb').read(), np.uint8).astype(bool); assert len(full) == n
    for fn in (fi, fo, exe, fpc):
        if os.path.exists(fn): os.remove(fn)
    os.rmdir(td)
    good = full & valid & start_ok & end_ok
    # ---- 9. connectivity of the good pieces (face-adjacent pieces sharing a point), intervals
    gi = np.nonzero(good)[0]
    Y0 = J * w; Y1 = Y0 + w; Z0 = K * w; Z1 = Z0 + w
    yc2 = (2 * J + 1) * w; zc2 = (2 * K + 1) * w          # twice the cell centres, lattice units
    parent = np.arange(n)
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    nedge = 0
    for axis in (0, 1):
        A_0, A_1, B_0, B_1 = (Y0, Y1, Z0, Z1) if axis == 0 else (Z0, Z1, Y0, Y1)
        ol = gi[np.lexsort((B_0[gi], A_1[gi]))]; orr = gi[np.lexsort((B_0[gi], A_0[gi]))]
        vl = A_1[ol]; vr = A_0[orr]
        for v_ in np.intersect1d(vl, vr).tolist():
            Ls = ol[np.searchsorted(vl, v_, 'left'):np.searchsorted(vl, v_, 'right')]
            Rs = orr[np.searchsorted(vr, v_, 'left'):np.searchsorted(vr, v_, 'right')]
            rb0 = B_0[Rs]; rb1 = B_1[Rs]
            lo_i = np.searchsorted(rb1, B_0[Ls], 'right'); hi_i = np.searchsorted(rb0, B_1[Ls], 'left')
            cntp = np.maximum(hi_i - lo_i, 0)
            if cntp.sum() == 0: continue
            pi_ = np.repeat(Ls, cntp)
            st_ = np.repeat(lo_i - np.r_[0, np.cumsum(cntp)[:-1]], cntp) + np.arange(int(cntp.sum()))
            pq_ = Rs[st_]
            plo = np.maximum(B_0[pi_], B_0[pq_]); phi = np.minimum(B_1[pi_], B_1[pq_])
            ok = phi > plo; pi_, pq_, plo, phi = pi_[ok], pq_[ok], plo[ok], phi[ok]
            conn = np.zeros(len(pi_), bool)
            for frac in (0.0, 0.5, 1.0):
                t2 = 2 * plo + frac * 2 * (phi - plo)      # twice the transverse lattice coordinate (exact)
                ab = []
                for p in (pi_, pq_):
                    ny = (2 * v_ - yc2[p]) if axis == 0 else (t2 - yc2[p]); nz = (t2 - zc2[p]) if axis == 0 else (2 * v_ - zc2[p])
                    dy = ny.astype(float) / (2 * Nf); dz = nz.astype(float) / (2 * Nf)
                    dyl, dyh, dzl, dzh = dn(dy), up(dy), dn(dz), up(dz)
                    ab.append((ilin(A0[p], aA[p], bA[p], dyl, dyh, dzl, dzh), ilin(B0[p], aB[p], bB[p], dyl, dyh, dzl, dzh)))
                (Ai, Bi), (Aq, Bq) = ab
                conn |= np.maximum(Ai[1], Aq[1]) <= np.minimum(Bi[0], Bq[0])
            nedge += int(conn.sum())
            for i, q in zip(pi_[conn].tolist(), pq_[conn].tolist()):
                ra, rb = find(i), find(q)
                if ra != rb: parent[ra] = rb
    roots = np.array([find(i) for i in gi.tolist()])
    # ---- 10. volume: piece = cell area * (B0 - A0); exact rational sum of lower bounds
    contrib = (wf * wf) * dn(B0 - A0)                      # = N^2 * (lower bound of the piece volume), exact scaling
    comps, cinv = np.unique(roots, return_inverse=True)
    approx = np.bincount(cinv, weights=contrib[gi]); big = int(np.argmax(approx))
    vol_big = sum(Fr(x) for x in contrib[gi[cinv == big]].tolist()) / (N * N)
    vol_all = sum(Fr(x) for x in contrib[gi].tolist()) / (N * N)
    claim = Fr(math.floor(vol_big * 10 ** 7), 10 ** 7)
    res = dict(pieces=n, core_full=int(full.sum()), valid=int(valid.sum()), start_ok=int(start_ok.sum()), end_ok=int(end_ok.sum()),
               good=int(good.sum()), components=len(comps), adjacency_edges=nedge, refined_poses=m + 1, quaternion_vs_euler=qdev,
               rot_rounding=dRf, max_abs_q=QR, max_eps=epsmax, input_errors=[ea, ed, ey, ez], E_core=float(E_core), E_block=float(E_block),
               margin_over_E=round(ETA / float(max(E_core, E_block))), volume_all_lower=float(vol_all), largest_component_lower=float(vol_big),
               seconds=round(time.time() - t0, 1), fraction_seconds=round(t1 - t0, 1), core_seconds=round(tc, 1))
    print(json.dumps(res, indent=1))
    print('CERTIFIED: V*(K_4) >= %d.%07d' % (claim.numerator // claim.denominator, (claim - claim.numerator // claim.denominator) * 10 ** 7))

if __name__ == '__main__':
    main()
