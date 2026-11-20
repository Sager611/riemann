# explore61 s4: from s1 npz files: relative in-W1 constant c_e = L Q/R along the hybrid optimizer, the pencil sup_e L e'(S-Om)e / e'Om e,
# the Lemma-B cosine R/sqrt((D^2-R) Q), and N-stability of the scaled head profile L*e.  usage: python3 s4_pencil.py
import glob, math, numpy as np
from scipy.linalg import eigh
fs = sorted(glob.glob('logs/s1_N*_Q*_K*.npz'), key=lambda f: [int(x) for x in f.replace('.npz','').split('_')[1:] for x in [x[1:]]])
prof = {}
for f in fs:
    z = np.load(f); N = int(f.split('_N')[1].split('_')[0]); Qg = int(f.split('_Q')[1].split('_')[0]); K = int(f.split('_K')[1].split('.')[0])
    L = math.log(N); S, Om, e = z['S'], z['Om'], z['e']; P, R = float(z['P']), float(z['R']); Qd = P - R; D2 = float(z['d2']) - P
    w = eigh(S - Om, Om, eigvals_only=True)
    cosB = R / math.sqrt((D2 - R) * Qd)
    print(f'N={N:5d} Q={Qg:2d} K0={K:2d}: L^2P={L*L*P:.4f} L^2R={L*L*R:.4f} L^3Q={L**3*Qd:.4f} c_e=LQ/R={L*Qd/R:.3f} '
          f'pencil L*lam_max={L*w.max():.2f} L*lam_min={L*w.min():.4f} cosB={cosB:.4f} L|e|_1={L*np.abs(e).sum():.3f} L|e|_2={L*np.linalg.norm(e):.3f}')
    prof[(N, Qg, K)] = (z['th'], L * e)
for (a, b) in [((640, 8, 5), (2000, 8, 5)), ((640, 16, 5), (2000, 16, 5)), ((160, 8, 5), (2000, 8, 5)), ((40, 8, 5), (2000, 8, 5))]:
    if a in prof and b in prof:
        x, y = prof[a][1], prof[b][1]
        print(f'profile L*e: N={a[0]} vs N={b[0]} (Q={a[1]}): corr {np.corrcoef(x, y)[0,1]:.4f}  |diff|/|e| {np.linalg.norm(x-y)/np.linalg.norm(y):.3f}')
th, le = prof[(2000, 16, 5)]
o = np.argsort(th)
print('L*e at N=2000, Q=16 (theta: value):', ' '.join(f'{th[i]:.3f}:{le[i]:+.3f}' for i in o[::3]))
