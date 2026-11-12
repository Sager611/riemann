"""Independent float cross-checks: (1) closed-form H(c) vs quadrature R(psi); (3) quadratic-form slack vs direct slack
at random windows. (The real-space check of the Poisson energy is crosscheck2.py.)
Run as: uv run ... python3 -c "import sys; sys.argv=['cover.py','227','153']; exec(open('crosscheck.py').read())" """
import json
import numpy as np
from tenwin import H_of, Khat
from winopt import Hc
from cover import Qmat, slack, S, P, NU0

rng = np.random.default_rng(1)
d = json.load(open(f"cover_{S}_{P}.json"))
for t in range(3):
    c = np.concatenate([[1.0], rng.normal(0, 0.02, 16)])
    print("H closed", Hc(c), " H quad", H_of(c))
x = np.array(d["x"][0])
Q = Qmat(x)
for t in range(2):
    cc = np.concatenate([[1.0], rng.normal(0, 0.01, 16)])
    ms = np.arange(1, S)
    kh = Khat(ms / S, cc); k0 = Khat(0.0, cc)[0]
    mu = np.exp(2j * np.pi * np.outer(ms, x) / S).sum(axis=1)
    e = (k0 * P * P + 2 * np.sum(kh * np.abs(mu) ** 2)) / S ** 2 - P / S
    print("slack direct", NU0 - H_of(cc) - e, " quadratic form", slack(cc, Q))
