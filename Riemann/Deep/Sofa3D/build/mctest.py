"""mctest.py -- Monte Carlo end-to-end test: sample points of certified sub-pieces and follow the
TRUE continuous motion (linear translation + uniform rotation about the fixed axis per smoothed step),
sampled at many times per step; every sample must lie in K."""
import sys, math, time, numpy as np
sys.path.insert(0, '.')
import sofa2 as S
sub = np.load(sys.argv[1]); path = sys.argv[2]; npts = int(sys.argv[3]); nt = int(sys.argv[4]); seed = int(sys.argv[5])
tt = float(sys.argv[6]) if len(sys.argv) > 6 else 0.005; tr = float(sys.argv[7]) if len(sys.argv) > 7 else 0.05
rng = np.random.default_rng(seed)
P = S.load_path('../gist/%s_path.json' % path)
sm = S.smooth(P, tt, tr)
R = S.euler_R(sm[:, 3:]); T = sm[:, :3]
ax, th = S.rel_rot(R[:-1], R[1:])
vol = 4 * sub[:, 2] * sub[:, 3] * (sub[:, 7] - sub[:, 4])
idx = rng.choice(len(sub), size=npts, p=vol / vol.sum())
Pc = sub[idx]
# half of the points on piece boundaries (faces), half interior
s = rng.random(npts); u = rng.uniform(-1, 1, npts); v = rng.uniform(-1, 1, npts)
face = rng.integers(0, 7, npts)   # 0..5 faces, 6 interior
s = np.where(face == 0, 0.0, np.where(face == 1, 1.0, s))
u = np.where(face == 2, -1.0, np.where(face == 3, 1.0, u))
v = np.where(face == 4, -1.0, np.where(face == 5, 1.0, v))
uu = u * Pc[:, 2]; vv = v * Pc[:, 3]
A = Pc[:, 4] + Pc[:, 5] * uu + Pc[:, 6] * vv; B = Pc[:, 7] + Pc[:, 8] * uu + Pc[:, 9] * vv
q = np.stack([A + s * (B - A), Pc[:, 0] + uu, Pc[:, 1] + vv], 1)

def inK(w, tol=1e-12):
    x, y, z = w[..., 0], w[..., 1], w[..., 2]
    zin = (z >= -tol) & (z <= 1 + tol); yin = (y >= -tol) & (y <= 1 + tol)
    a = (x <= 1 + tol) & yin & zin
    b = (x >= -tol) & (x <= 1 + tol) & (y >= -tol) & (y <= 6 + tol) & zin
    c = (x >= -tol) & (x <= 1 + tol) & (y >= 5 - tol) & (y <= 6 + tol) & (z >= -tol)
    return a | b | c

ts_ = np.linspace(0, 1, nt + 1)
bad = np.zeros(npts, bool); worst = np.zeros(npts)
t0 = time.time()
w0 = q @ R[0].T + T[0]
startSX = (w0[:, 0] <= 0) & (w0[:, 1] >= 0) & (w0[:, 1] <= 1) & (w0[:, 2] >= 0) & (w0[:, 2] <= 1)
for k in range(len(sm) - 1):
    Rt = S.rodrigues(np.broadcast_to(ax[k], (nt + 1, 3)), ts_ * th[k]) @ R[k]     # (nt+1,3,3)
    Tt = (1 - ts_)[:, None] * T[k] + ts_[:, None] * T[k + 1]
    w = np.einsum('tij,nj->tni', Rt, q) + Tt[:, None, :]
    ok = inK(w).all(0)
    bad |= ~ok
wE = q @ R[-1].T + T[-1]
endSZ = (wE[:, 0] >= 0) & (wE[:, 0] <= 1) & (wE[:, 1] >= 5) & (wE[:, 1] <= 6) & (wE[:, 2] >= 0)
endR = np.abs(Rt[-1] - R[-1]).max()
print('points', npts, 'samples/step', nt + 1, 'steps', len(sm) - 1, 'violations', int(bad.sum()), 'faces of bad', np.bincount(face[bad], minlength=7).tolist(),
      'start in SX', bool(startSX.all()), 'end in SZ', bool(endSZ.all()), 'min end z', float(wE[:, 2].min()), 'rot err', endR, 'time %.0fs' % (time.time() - t0))
