
# sanity.py: compare the assembled lower bound with the exact c_n for given m (must never exceed it)
import sys, math
import flint
from flint import arb, acb, fmpq, arb_series
import assembly as A
from fh2lib import Psie, S, PI
from math import comb, factorial
m = int(sys.argv[1]); nmax = int(sys.argv[2]); step = int(sys.argv[3]) if len(sys.argv) > 3 else 1
nu0 = arb('203.5'); A.load(nu0, 'tab203.json', 'lg203.pkl')
flint.ctx.prec = 64 + 4 * m + 2 * nmax
flint.ctx.cap = nmax + 10
nu = arb(m) + arb(1) / 2
b = [fmpq((-1)**j * comb(m, j) * factorial(2*m - j), factorial(2*m)) for j in range(m + 1)]
c = (1 / arb_series([arb(x) for x in b], prec=nmax + 1)).coeffs()
th0 = arb(factorial(2*m)) / (factorial(m) * 2**m)
Lam = th0 * nu.sqrt() * nu**(-nu)
logLam = Lam.log()
worst = None; viol = 0
for n in range(8, nmax + 1, step):
    tau = arb(n) / nu
    # best good l at this exact tau (as in check_box with a degenerate box)
    best = None
    for l in range(0, n // (2*m+1) + 1):
        a = 2*l + 1; t = tau - 2*l
        if not (t > 0) and l > 0: continue
        r = A.saddle(a, t)
        rf = float(r.mid())
        gf = A.good_factor(l, rf * 0.999999, rf * 1.000001 + (1e-9 if rf == 0 else 0), nu0)
        if gf is None: continue
        ph = A.Phi(a, t)
        val = nu * ph + gf      # log of kappa^l * J_l lower bound, at actual nu (gf at nu0 is smaller: ok)
        if best is None or val > best[0]: best = (val, l)
    if best is None:
        print('n', n, 'no good term'); continue
    # exact normalized value: (2 nu)^n c_n / Lam
    lhs = (2 * nu).log() * n + c[n].log() - logLam
    gap = lhs - best[0]
    if not (gap > 0): viol += 1; print('VIOLATION n=%d l=%d gap=%s' % (n, best[1], gap.str(5)))
    gf_ = float(gap.mid())
    if worst is None or gf_ < worst[0]: worst = (gf_, n, best[1])
print('m=%d n<=%d: violations %d ; min log(true/LB) = %.4f at n=%d (l=%d)' % (m, nmax, viol, worst[0], worst[1], worst[2]))
