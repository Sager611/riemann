# explore18 T3 diagnostic: where are the off-line zeros of the symmetrized models? (near the line or near Re s = 0, 1)
import sys
from t3_bulk import *
mp.mp.dps = 15
for kind, par, s_hi in (('euler', 3, 0.97), ('euler', 29, 0.97), ('zpart', 10, 2.0), ('dhpart', 10, 2.0), ('dheuler', 29, 0.97)):
    G, B, F = make(kind, par)
    a = zeros_in_rect(F, 0.505, 0.8, 1, 250, 0.1)
    b = zeros_in_rect(F, 0.8, s_hi, 1, 250, 0.1)
    print(f"{kind:8s} {par:3d}: zeros of F with 0.505 < Re s < 0.8, 1 < t < 250: {mp.nstr(a, 5):>6s};  0.8 < Re s < {s_hi}: {mp.nstr(b, 5):>6s}", flush=True)
