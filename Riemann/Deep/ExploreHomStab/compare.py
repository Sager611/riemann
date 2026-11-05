# explore9: compare sum* L(1/2,chi_8d)^k (d odd squarefree <= X) with the CFKRS recipe polynomial sum.
# Usage: python3 compare.py L1e6.npz Q1.json Q2.json Q3.json
import sys, json
import numpy as np, mpmath as mp
z = np.load(sys.argv[1]); d = z['d'].astype(float); Lv = z['L']
Qs = [json.load(open(f)) for f in sys.argv[2:]]
ell = np.log(8 * d / np.pi)
def W(u):  # smooth bump supported on (1/2, 1)
    s = 4 * u - 3; out = np.zeros_like(u); m = np.abs(s) < 1
    out[m] = np.exp(1 - 1 / (1 - s[m] ** 2)); return out
What = float(mp.quad(lambda u: mp.e ** (1 - 1 / (1 - (4 * u - 3) ** 2)) * u ** (-0.25), [0.5, 0.75, 1]))
# Diaconu-Whitehead residue of Z_0 at s = 3/4 (cubic moment, family chi_{2d} = chi_{8d})
def Pfac(x): return (1 - x) ** 5 * (1 + x) * (1 + 4 * x + 11 * x ** 2 + 10 * x ** 3 - 11 * x ** 4 + 11 * x ** 6 - 4 * x ** 7 - x ** 8)
s = np.ones(10 ** 7 + 1, bool); s[:2] = False
for i in range(2, int(10 ** 3.5) + 1):
    if s[i]: s[i * i::i] = False
pp = np.nonzero(s)[0]; pp = pp[pp > 2].astype(float)
logprod = np.sum(np.log(Pfac(pp ** -0.5)))  # tail beyond 1e7 ~ -14 sum p^{-3/2} ~ 1e-3 relative
res = float(mp.mpf(9) / (256 * mp.pi) * mp.mpf(2) ** 0.25 * (-181 + 128 * mp.sqrt(2)) * mp.gamma(0.25) ** 4 * mp.zeta(0.5) ** 7) * np.exp(logprod)
print(f"DW residue at 3/4: {res:.6f}; W-hat(3/4) = {What:.6f}")
Xs = [10 ** 3, 10 ** 4, 10 ** 5, 3 * 10 ** 5, 10 ** 6]
for Q in Qs:
    k = Q['k']; c = np.array(Q['coef_re'])
    pred = np.polyval(c[::-1], ell); lead = c[-1] * ell ** len(c[:-1]) if False else c[-1] * ell ** (len(c) - 1)
    data = Lv ** k
    print(f"\nk={k}: Q_k coefficients (powers L^0..): " + ", ".join(f"{x:.10g}" for x in c))
    print("  sharp:  X, data, recipe, data/recipe-1, diff/X^0.5, leading-only/recipe-1")
    for X in Xs:
        m = d <= X; D = data[m].sum(); P = pred[m].sum(); Lo = lead[m].sum()
        print(f"   {X:>8d} {D:.8e} {P:.8e} {D/P-1:+.3e} {(D-P)/X**0.5:+.3e} {Lo/P-1:+.3e}")
    print("  smooth W(d/X):  X, data, recipe, rel diff, diff/X^0.5" + (", DW 3/4-term" if k == 3 else ""))
    for X in [10 ** 4, 3 * 10 ** 4, 10 ** 5, 3 * 10 ** 5, 10 ** 6]:
        w = W(d / X); D = (data * w).sum(); P = (pred * w).sum()
        extra = f" {res*What*X**0.75:+.3e}" if k == 3 else ""
        print(f"   {X:>8d} {D:.8e} {P:.8e} {D/P-1:+.3e} {(D-P)/X**0.5:+.3e}{extra}")
