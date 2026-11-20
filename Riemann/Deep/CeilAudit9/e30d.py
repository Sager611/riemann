# ceilaudit9 E30d: anchored LP (block ends at S = S_1 = 0) at t -> infinity (S/D <= 0.594, S_1/D <= 0.371), w = 0.03a, 16 phases.
import sys, math
sys.path.insert(0, 'Riemann/Deep/CeilAudit9')
from e30 import single_anchored
for a in (1.2, 1.3, 1.4, 1.5):
    print(f"t->inf anchored a={a} (c={math.exp(2*a):.2f}): W_norm {single_anchored(a, 0.03, 1.0, 0.0946*2*math.pi, 0.059*2*math.pi, npsi=16):+.4f}", flush=True)
