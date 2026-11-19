# Transitions (iv): best cell margins just below and above loglog q = 3, 3.13, 4, with lwcheck.py's own F (B mode).
import sys
sys.argv = ['lwcheck.py', 'B']
exec(open('lwcheck.py').read().split('# ------------------------------------------------------------------ main')[0])
def marg(c, i, j, den=1000):
    ell = iv.exp(I(j, den)); r, P = optimise(float(ell.a)); v = F(ell, Pq(r))
    return float((v - iv.log(I(*c)/I(i, den))).a), P
for lab, c, cells in (('0.11 /max(3,.)', (11, 100), ((2990, 3000), (3000, 3001))),
                      ('0.115', (115, 1000), ((3120, 3121), (3125, 3126), (3129, 3130), (3130, 3131))),
                      ('0.116', (116, 1000), ((3990, 4000), (4000, 4001), (3800, 3801)))):
    for i, j in cells:
        m, P = marg(c, i, j); print(f'{lab:15s} cell [{i/1000}, {j/1000}]: best margin {m:+.6f}  P={P}')
ell3 = iv.exp(I(3)); r, P = optimise(float(ell3.a)); v = F(ell3, Pq(r))
print(f'Lam <= 3 regime: F(e^3) - log(0.11/3) = {float((v - iv.log(I(11, 300))).a):+.6f}')
