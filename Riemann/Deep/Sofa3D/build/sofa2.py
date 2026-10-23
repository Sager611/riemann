"""sofa2.py -- motions (checks) + pieces for carve2 (rigorous s-carve of sheared prisms).
Body coordinates q = (x,y,z) in the block [-4,0]x[0,1]x[0,1]; a pose (x,y,z,yaw,pitch,roll)
places q at R q + T, R = Rx(roll) Ry(pitch) Rz(yaw) (scipy 'zyx' extrinsic, degrees)."""
import json, math, os, struct, subprocess, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CARVE2 = os.path.join(HERE, "carve2")
QR = math.sqrt(18.0)            # max |q| over the block
ETA = 1e-9                      # absolute safety margin (world units)

def load_path(fn):
    d = json.load(open(fn))
    return np.array(d["poses"] if isinstance(d, dict) else d, dtype=np.float64)

def smooth(path, tt, tr):
    """Exact replica of the pipeline's smooth_path_to_tolerances (float64)."""
    out = [path[0].tolist()]
    for i in range(len(path) - 1):
        p1 = path[i]; p2 = path[i + 1]
        td = np.linalg.norm(p2[:3] - p1[:3]); rd = np.max(np.abs(p2[3:] - p1[3:]))
        nt = int(np.ceil(td / tt)) if td > tt else 0
        nr = int(np.ceil(rd / tr)) if rd > tr else 0
        n = max(1, nt, nr)
        if n > 1:
            for j in range(1, n):
                out.append((p1 + (j / n) * (p2 - p1)).tolist())
        out.append(list(p2))
    return np.array(out)

def euler_R(ang_deg):
    a = np.radians(np.asarray(ang_deg, dtype=np.float64))
    y, p, r = a[..., 0], a[..., 1], a[..., 2]
    cy, sy, cp, sp, cr, sr = np.cos(y), np.sin(y), np.cos(p), np.sin(p), np.cos(r), np.sin(r)
    R = np.empty(a.shape[:-1] + (3, 3))
    R[..., 0, 0] = cp * cy;  R[..., 0, 1] = -cp * sy; R[..., 0, 2] = sp
    R[..., 1, 0] = sr * sp * cy + cr * sy; R[..., 1, 1] = -sr * sp * sy + cr * cy; R[..., 1, 2] = -sr * cp
    R[..., 2, 0] = -cr * sp * cy + sr * sy; R[..., 2, 1] = cr * sp * sy + sr * cy; R[..., 2, 2] = cr * cp
    return R

def rel_rot(R1, R2):
    """unit axis and angle in [0,pi] of R2 R1^T (batched)."""
    Rr = R2 @ np.swapaxes(R1, -1, -2)
    v = np.stack([Rr[..., 2, 1] - Rr[..., 1, 2], Rr[..., 0, 2] - Rr[..., 2, 0], Rr[..., 1, 0] - Rr[..., 0, 1]], -1)
    nv = np.linalg.norm(v, axis=-1)
    th = np.arctan2(nv / 2, (np.trace(Rr, axis1=-2, axis2=-1) - 1) / 2)
    ax = np.where(nv[..., None] > 0, v / np.where(nv > 0, nv, 1)[..., None], np.array([1.0, 0, 0]))
    return ax, th

def rodrigues(ax, phi):
    K = np.zeros(ax.shape[:-1] + (3, 3))
    K[..., 0, 1] = -ax[..., 2]; K[..., 0, 2] = ax[..., 1]; K[..., 1, 0] = ax[..., 2]
    K[..., 1, 2] = -ax[..., 0]; K[..., 2, 0] = -ax[..., 1]; K[..., 2, 1] = ax[..., 0]
    s = np.sin(phi)[..., None, None]; c = (1 - np.cos(phi))[..., None, None]
    return np.eye(3) + s * K + c * (K @ K)

def eps_of(theta):
    """chord deviation bound for uniform rotation by theta (+ linear translation):
    |p(t) - chord(t)| <= r (theta^2/8 + theta^3/24), r <= sqrt(18)."""
    return QR * (theta ** 2 / 8 + theta ** 3 / 24) * (1 + 1e-6) + 1e-12

def maps(R, T):
    R = np.asarray(R); T = np.asarray(T)
    return np.concatenate([R[..., :, 0], R[..., :, 1], R[..., :, 2], T], -1)

def pack(typ, m1, m2, eps):
    n = len(typ)
    return np.concatenate([np.asarray(typ, float)[:, None], m1, m2, np.asarray(eps, float)[:, None]], 1).reshape(n, 26)

def checks_stair(sm):
    """the pipeline's staircase: per step x-, y-, z-translation at R_k, then rotation at T_{k+1}."""
    R = euler_R(sm[:, 3:]); T = sm[:, :3]; n = len(sm) - 1
    ax, th = rel_rot(R[:-1], R[1:])
    tX = np.stack([T[1:, 0], T[:-1, 1], T[:-1, 2]], 1); tY = np.stack([T[1:, 0], T[1:, 1], T[:-1, 2]], 1)
    z = np.zeros((n, 12))
    cX = pack(np.zeros(n), maps(R[:-1], tX), z, np.zeros(n))
    cY = pack(np.zeros(n), maps(R[:-1], tY), z, np.zeros(n))
    cR = pack(np.ones(n), maps(R[:-1], T[1:]), maps(R[1:], T[1:]), eps_of(th))
    C = np.stack([cX, cY, cR], 1).reshape(3 * n, 26)
    c0 = pack([0], maps(R[:1], T[:1]), np.zeros((1, 12)), [0])
    return np.vstack([c0, C]), th

def checks_linear(sm):
    """simultaneous motion: linear translation + uniform rotation about the fixed axis per step."""
    R = euler_R(sm[:, 3:]); T = sm[:, :3]; n = len(sm) - 1
    ax, th = rel_rot(R[:-1], R[1:])
    C = pack(np.ones(n), maps(R[:-1], T[:-1]), maps(R[1:], T[1:]), eps_of(th))
    c0 = pack([0], maps(R[:1], T[:1]), np.zeros((1, 12)), [0])
    return np.vstack([c0, C]), th

def checks_hybrid(sm, order=(0, 1, 2)):
    """per step: three substeps, each translating along one axis (in 'order') while
    rotating by theta/3 about the fixed axis."""
    R = euler_R(sm[:, 3:]); T = sm[:, :3]; n = len(sm) - 1
    ax, th = rel_rot(R[:-1], R[1:])
    Rs = [R[:-1], rodrigues(ax, th / 3) @ R[:-1], rodrigues(ax, 2 * th / 3) @ R[:-1], R[1:]]
    Ts = [T[:-1].copy()]
    cur = T[:-1].copy()
    for a in order:
        cur = cur.copy(); cur[:, a] = T[1:, a]; Ts.append(cur)
    e = eps_of(th / 3)
    cs = [pack(np.ones(n), maps(Rs[j], Ts[j]), maps(Rs[j + 1], Ts[j + 1]), e) for j in range(3)]
    C = np.stack(cs, 1).reshape(3 * n, 26)
    c0 = pack([0], maps(R[:1], T[:1]), np.zeros((1, 12)), [0])
    return np.vstack([c0, C]), th

def run(checks, pieces, eta=ETA, qx=4.0, qy=1.0, qz=1.0, ts=()):
    pieces = np.ascontiguousarray(pieces, dtype=np.float64)
    td = tempfile.mkdtemp(dir=os.path.join(HERE, "tmp"))
    fi = os.path.join(td, "in.bin"); fo = os.path.join(td, "out.bin")
    with open(fi, "wb") as f:
        f.write(np.array([len(checks), len(pieces), eta, qx, qy, qz, len(ts)] + list(ts), np.float64).tobytes())
        f.write(np.ascontiguousarray(checks, np.float64).tobytes()); f.write(pieces.tobytes())
    t0 = time.time()
    p = subprocess.run([CARVE2, fi, fo], capture_output=True, text=True)
    dt = time.time() - t0
    os.remove(fi)
    if p.returncode != 0:
        raise RuntimeError("carve2 rc=%d %s" % (p.returncode, p.stderr))
    raw = open(fo, "rb").read(); os.remove(fo); os.rmdir(td)
    NP, tot, _, ovf = struct.unpack("4i", raw[:16])
    cnt = np.frombuffer(raw, dtype=np.int32, count=NP, offset=16).copy()
    iv = np.frombuffer(raw, dtype=np.float64, offset=16 + 4 * NP).reshape(-1, 2).copy()
    assert len(iv) == tot == int(cnt.sum())
    return cnt, iv, dt, p.stderr.strip()

def piece_len(cnt, iv):
    off = np.zeros(len(cnt) + 1, np.int64); np.cumsum(cnt, out=off[1:])
    L = iv[:, 1] - iv[:, 0]
    return np.add.reduceat(np.concatenate([L, [0.0]]), off[:-1]) * (cnt > 0), off

def prism_pieces(nc):
    h = 1.0 / nc; c = (np.arange(nc) + 0.5) * h
    Y, Z = np.meshgrid(c, c, indexing="ij")
    P = np.zeros((nc * nc, 10)); P[:, 0] = Y.ravel(); P[:, 1] = Z.ravel(); P[:, 2] = P[:, 3] = h / 2
    P[:, 4] = -4.0; P[:, 7] = 0.0
    return P

def thin_pieces(ys, zs):
    Y, Z = np.meshgrid(ys, zs, indexing="ij")
    P = np.zeros((Y.size, 10)); P[:, 0] = Y.ravel(); P[:, 1] = Z.ravel(); P[:, 4] = -4.0; P[:, 7] = 0.0
    return P

def volume(P, cnt, iv):
    s, _ = piece_len(cnt, iv)
    return math.fsum((s * (P[:, 7] - P[:, 4]) * (2 * P[:, 2]) * (2 * P[:, 3])).tolist())

def merge_gaps(L, g):
    if len(L) < 2 or g <= 0: return L
    out = [list(L[0])]
    for a, b in L[1:]:
        if a - out[-1][1] < g: out[-1][1] = max(out[-1][1], b)
        else: out.append([a, b])
    return np.array(out)

def corner_lattice(nc):
    g = np.arange(nc + 1) / nc
    return thin_pieces(g, g)

def sheared_pieces(nc, cnt, iv, delta=0.0, fallback=True, gmerge=0.0):
    """Build sheared prisms per cell from thin-line intervals on the (nc+1)^2 corner
    lattice (line (j,k) at y=j/nc, z=k/nc, index j*(nc+1)+k; x = -4 + 4 s).
    Returns pieces (n,10) and cell index per piece.  Correctness is NOT assumed:
    the s-carve re-verifies every piece."""
    h = 1.0 / nc; hh = h / 2; n1 = nc + 1
    off = np.zeros(len(cnt) + 1, np.int64); np.cumsum(cnt, out=off[1:])
    X = -4.0 + 4.0 * iv
    out = []; cells = []
    J, K = np.meshgrid(np.arange(nc), np.arange(nc), indexing="ij"); J = J.ravel(); K = K.ravel()
    cidx = [J * n1 + K, (J + 1) * n1 + K, J * n1 + K + 1, (J + 1) * n1 + K + 1]   # (--),(+-),(-+),(++)
    single = np.ones(len(J), bool)
    for c in cidx: single &= (cnt[c] == 1)
    def mk(a, b, yc, zc):
        # a, b: (m,4) corner values in order (--),(+-),(-+),(++)
        ma = a.mean(1); ua = (a[:, 1] + a[:, 3] - a[:, 0] - a[:, 2]) / (4 * hh); va = (a[:, 2] + a[:, 3] - a[:, 0] - a[:, 1]) / (4 * hh)
        ta = (a[:, 0] - a[:, 1] - a[:, 2] + a[:, 3]) / 4
        mb = b.mean(1); ub = (b[:, 1] + b[:, 3] - b[:, 0] - b[:, 2]) / (4 * hh); vb = (b[:, 2] + b[:, 3] - b[:, 0] - b[:, 1]) / (4 * hh)
        tb = (b[:, 0] - b[:, 1] - b[:, 2] + b[:, 3]) / 4
        A0 = ma + np.abs(ta) + delta; B0 = mb - np.abs(tb) - delta
        A0 = np.maximum(A0, -4.0 + np.abs(ua) * hh + np.abs(va) * hh)
        B0 = np.minimum(B0, 0.0 - np.abs(ub) * hh - np.abs(vb) * hh)
        P = np.stack([yc, zc, np.full(len(yc), hh), np.full(len(yc), hh), A0, ua, va, B0, ub, vb], 1)
        # B - A > 0 at all corners
        ok = np.ones(len(yc), bool)
        for sy in (-1, 1):
            for sz in (-1, 1):
                ok &= (B0 - A0) + (ub - ua) * sy * hh + (vb - va) * sz * hh > 1e-9
        return P, ok
    yc = (J + 0.5) * h; zc = (K + 0.5) * h
    # fast path: every corner has exactly one interval
    sel = np.nonzero(single)[0]
    a = np.stack([X[off[c[sel]], 0] for c in cidx], 1); b = np.stack([X[off[c[sel]], 1] for c in cidx], 1)
    ov = np.minimum.reduce([b[:, i] for i in range(4)]) - np.maximum.reduce([a[:, i] for i in range(4)])
    P, ok = mk(a, b, yc[sel], zc[sel]); ok &= ov > 0
    out.append(P[ok]); cells.append(sel[ok])
    done = np.zeros(len(J), bool); done[sel[ok]] = True
    # slow path
    for q in np.nonzero(~single)[0]:
        lists = [merge_gaps(X[off[c[q]]:off[c[q] + 1]], gmerge) for c in cidx]
        if any(len(L) == 0 for L in lists):
            continue
        tups = set()
        for ci in range(4):
            for I in range(len(lists[ci])):
                lo, hi = lists[ci][I]
                tup = []
                for cj in range(4):
                    if cj == ci: tup.append(I); continue
                    L = lists[cj]
                    o = np.minimum(L[:, 1], hi) - np.maximum(L[:, 0], lo)
                    bi = int(np.argmax(o))
                    if o[bi] <= 0: tup = None; break
                    tup.append(bi)
                if tup is not None: tups.add(tuple(tup))
        if not tups: continue
        tups = sorted(tups)
        a = np.array([[lists[c][t[c], 0] for c in range(4)] for t in tups]); b = np.array([[lists[c][t[c], 1] for c in range(4)] for t in tups])
        P, ok = mk(a, b, np.full(len(tups), yc[q]), np.full(len(tups), zc[q]))
        P = P[ok]
        if len(P) == 0: continue
        # resolve overlaps: keep largest first, add others only if x-ranges over the cell are disjoint
        xlo = P[:, 4] - np.abs(P[:, 5]) * hh - np.abs(P[:, 6]) * hh; xhi = P[:, 7] + np.abs(P[:, 8]) * hh + np.abs(P[:, 9]) * hh
        order = np.argsort(-(P[:, 7] - P[:, 4]))
        keep = []
        for i in order:
            if all(xhi[i] < xlo[j] or xhi[j] < xlo[i] for j in keep): keep.append(i)
        out.append(P[keep]); cells.append(np.full(len(keep), q)); done[q] = True
    PS = np.vstack(out); CS = np.concatenate(cells)
    if fallback:
        rest = np.nonzero(~done)[0]
        if len(rest):
            F = np.zeros((len(rest), 10)); F[:, 0] = yc[rest]; F[:, 1] = zc[rest]; F[:, 2] = F[:, 3] = hh; F[:, 4] = -4.0
            PS = np.vstack([PS, F]); CS = np.concatenate([CS, rest])
    return PS, CS

