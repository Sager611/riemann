# validate weilleg against WeilSimple's Arb Fourier-basis matrix (eigs_quick values pasted by the shell script) and test convergence
import numpy as np, sys
sys.path.insert(0, '.')
from weilleg import WeilLeg
for c in [2, 3, 5]:
    a = np.log(c)/2
    out = []
    for K in [8, 12, 16, 24, 32]:
        W = WeilLeg(a, K).build()
        ev = np.linalg.eigvalsh(W.Q)
        out.append((K, ev[0], ev[1], np.abs(W.G - np.eye(K)).max()))
    print(f"c={c} a={a:.6f}:", ' | '.join(f"K={K}: l0={e0:.10e} l1={e1:.6e} (G err {g:.1e})" for K, e0, e1, g in out))
