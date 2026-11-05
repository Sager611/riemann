# bb.py: verify the torus envelope min_theta Phi(theta) >= -beta_t for the S-unit extension found by torus_lp.py, by
# branch-and-bound on T^D with a second-order Taylor bound (float64; margin >> rounding). Phi = -2 sum w cos(V.th) + 2 sum e cos(F.th).
import sys, math, numpy as np
d = np.load(sys.argv[1]); slack = float(sys.argv[2])
F, e, beta, V, W = d['F'], d['e'], float(d['beta']), d['Vfix'], d['Wfix']
A = np.vstack([V, F]) if len(F) else V          # frequency vectors
C = np.concatenate([-2 * W, 2 * e])             # Phi = sum C_j cos(A_j . th)
D = A.shape[1]
bt = beta + slack
def check(cent, h):
    ph = cent @ A.T                              # (n, J)
    val = np.cos(ph) @ C
    grad = -(np.sin(ph) * C) @ A                # (n, D)
    lin = np.abs(grad) @ np.full(D, 1.0) * h     # |grad| . h (h same in all dims)
    quad = 0.5 * np.sum(np.abs(C) * (np.abs(A).sum(axis=1) * h) ** 2)
    lower = val - lin - quad - 1e-9
    return lower >= -bt, val
# symmetry: Phi(-th) = Phi(th) -> half torus suffices, but keep it simple: full torus [-pi, pi)^D
n0 = 24; h = math.pi / n0
g = (np.arange(n0) + 0.5) * (2 * math.pi / n0) - math.pi
boxes = np.array(np.meshgrid(*([g] * D), indexing='ij')).reshape(D, -1).T
total = 0; minval = 1e9; level = 0
while len(boxes):
    bad = []
    for i in range(0, len(boxes), 100000):
        ok, val = check(boxes[i:i + 100000], h)
        minval = min(minval, float(val.min()))
        bad.append(boxes[i:i + 100000][~ok])
    total += len(boxes)
    bad = np.vstack(bad)
    if np.any(bad.shape[0] and (np.cos(bad @ A.T) @ C).min() < -bt):
        print('ENVELOPE VIOLATED at a box centre: min', (np.cos(bad @ A.T) @ C).min(), 'vs', -bt); break
    print(f'level {level}: boxes {len(boxes)}, unresolved {len(bad)}, h = {h:.2e}, min centre value {minval:.4f}', flush=True)
    if len(bad) == 0: break
    if len(bad) > 3_000_000 or level > 12: print('too many boxes; stop'); break
    h /= 2; level += 1
    offs = np.array(np.meshgrid(*([[-h, h]] * D), indexing='ij')).reshape(D, -1).T
    boxes = (bad[:, None, :] + offs[None, :, :]).reshape(-1, D)
else:
    pass
print(f'beta(LP) = {beta:.4f}, verified envelope: Phi >= -{bt:.4f} on the whole torus' if len(bad) == 0 else 'not verified', f'(boxes checked {total})')
