"""mkcert.py -- write a certificate: pieces (planes, cells, meta) + the pose table of the motion (T, Q).
Usage: python3 mkcert.py <run.npz from run_adapt4.py | existing certificate.npz> <path.json> <out.npz>
Pose table: the pipeline's smoothing of the path at the certificate's tolerances; Q_k = quaternion of
Rx(roll) Ry(pitch) Rz(yaw), signs chosen so that consecutive quaternions have positive dot product."""
import json, math, sys
import numpy as np

def smooth(path, tt, tr):
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
    h = lambda d: math.radians(d) / 2
    return qmul((math.cos(h(roll)), math.sin(h(roll)), 0.0, 0.0),
                qmul((math.cos(h(pitch)), 0.0, math.sin(h(pitch)), 0.0), (math.cos(h(yaw)), 0.0, 0.0, math.sin(h(yaw)))))

src, pathfile, out = sys.argv[1:4]
z = np.load(src)
if 'sub' in z.files:
    planes = z['sub'][:, 4:10]; cells = z['cells'].astype(np.int32); meta = z['meta']
else:
    planes = z['planes']; cells = z['cells'].astype(np.int32); meta = z['meta']
d = json.load(open(pathfile)); P = np.array(d['poses'] if isinstance(d, dict) else d, float)
sm = smooth(P, float(meta[3]), float(meta[4]))
Q = np.array([euler_quat(*p[3:]) for p in sm])
for k in range(1, len(Q)):
    if np.dot(Q[k - 1], Q[k]) < 0: Q[k] = -Q[k]
np.savez_compressed(out, planes=planes, cells=cells, meta=meta, T=sm[:, :3].copy(), Q=Q)
print('certificate', out, 'pieces', len(planes), 'poses', len(sm))
