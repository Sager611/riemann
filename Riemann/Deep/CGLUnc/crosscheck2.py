"""Real-space check of the Poisson energy formula: e(X) = (1/S) sum_{i in period} sum_{j != i, all copies} k(x_j - x_i)^2,
k(x) = psihat(2 pi x)/int psi in closed form; |d| <= R summed exactly, |d| > R by the far-field average."""
import json, sys
import numpy as np
from tenwin import TEN_C, OM, int_psi, psi

S, P = 227, 153
d = json.load(open(f"cover_{S}_{P}.json"))
x = np.array(d["x"][0]); c = TEN_C
I1 = int_psi(c)


def k(xx):
    z = 2 * np.pi * xx
    out = np.zeros_like(xx)
    for cj, w in zip(c, OM):
        a = z - w; b = z + w
        ta = np.where(np.abs(a) < 1e-12, 0.5, np.sin(a / 2) / np.where(np.abs(a) < 1e-12, 1, a))
        out += cj * (ta + np.sin(b / 2) / b)
    return out / I1


R = 4000.0
nmax = int(R / S) + 2
E = 0.0
for xi in x:
    dd = (x[None, :] - xi + S * np.arange(-nmax, nmax + 1)[:, None]).ravel()
    dd = dd[(np.abs(dd) > 1e-9) & (np.abs(dd) <= R)]
    E += np.sum(k(dd) ** 2)
psi_half = psi(np.array([0.5]), c)[0]
tail = P * (P / S) * psi_half ** 2 / (np.pi ** 2 * R * I1 ** 2)
from cover import Qmat, NU0
e_p = c @ Qmat(x) @ c - NU0
print("real-space e =", (E + tail) / S, "(tail part", tail / S, ")  Poisson e =", e_p)
