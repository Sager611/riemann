# T6c: band-limited Wirtinger gain vs run length m of maximal gaps, symmetric bands (real F = Z |M|^2: [-(1+2v)pi, (1+2v)pi], v = theta)
import numpy as np
from t5_bandwirt import R_config
PI = np.pi
for kap in [3.324]:
    for B in [1.0, 1.5, 2.0]:
        out = []
        for n, m in [(20, 1), (40, 1), (40, 2), (40, 3), (40, 5), (40, 8), (40, 11)]:
            rest = n - m*kap; k = n - m
            x = np.concatenate([kap*np.arange(m + 1), m*kap + rest*np.arange(1, k)/k])
            r = R_config(x, n, (-B*PI, B*PI), nq=121)
            out.append(f'm={m}(n={n}): {r/(PI/kap)**2:.4f}')
        print(f'kappa={kap} band [-{B}pi,{B}pi]: ' + '  '.join(out), flush=True)
