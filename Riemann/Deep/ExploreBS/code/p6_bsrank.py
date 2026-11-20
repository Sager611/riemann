# p6_bsrank.py (explore57): Birman-Schwinger order |K-| + |J-| of Chekmachev's even-sector certificates (his Table A.2,
# Zenodo 23848597, p. 34) against the Bohr-Jessen Weyl count W(c) = 2 ln(c) E[exp X_c] (p5_lattice.py).  No fit: c = 71 and
# every c are compared with the same formula.
import math, numpy as np, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from p5_lattice import EexpX
KM = {2: (2, 1), 3: (5, 1), 4: (6, 1), 5: (8, 1), 6: (10, 4), 7: (16, 1), 8: (21, 4), 9: (26, 4), 10: (21, 12), 11: (51, 4),
      12: (52, 18), 13: (78, 8), 14: (66, 35), 15: (56, 40), 16: (105, 28), 17: (166, 21), 18: (148, 64), 19: (271, 28),
      20: (214, 97), 21: (183, 121), 22: (202, 138), 23: (409, 48), 24: (345, 180), 25: (423, 135), 26: (340, 255),
      27: (497, 176), 28: (448, 262), 29: (955, 129), 30: (849, 416), 31: (1465, 198), 45: (2605, 1819), 47: (6301, 869),
      71: (33702, 4675)}
rows = []
for c, (km, jm) in KM.items():
    W = 2*math.log(c)*EexpX(c); o = km + jm; rows.append((c, o, W, o/W))
    print(f'c={c:3d} order={o:6d} W_BJ={W:10.1f} ratio={o/W:6.3f}', flush=True)
r = np.array([x[3] for x in rows if x[0] >= 10])
print(f'c >= 10: ratio min {r.min():.3f} max {r.max():.3f} geometric mean {math.exp(np.log(r).mean()):.3f}; '
      f'log-log slope of order vs W over c = 10..71: {np.polyfit(np.log([x[2] for x in rows if x[0] >= 10]), np.log([x[1] for x in rows if x[0] >= 10]), 1)[0]:.3f}')
for c in [100, 150, 200]:
    print(f'prediction c={c}: W_BJ = {2*math.log(c)*EexpX(c):.3g}')
