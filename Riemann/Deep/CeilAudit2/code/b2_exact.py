# B2 exactness check: BFlat's fixed point at M2 cost 2 vs B_2 = (777 - sqrt(1081))/852 at 30 digits (mpmath findroot in the sympy pipeline).
import os, io, contextlib
src = open(os.path.join(os.path.dirname(__file__), '..', '..', 'BFlat', 'bflat_sens.py')).read()
g = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src[:src.index("BF = mm.mpf")], 'bflat_sens_defs', 'exec'), g)
mm, Q, pub, lam, k = g['mm'], g['Q'], g['pub'], g['lam'], g['k']
mm.mp.dps = 40
for nm in ('E', 'TS', 'RS', 'coef_th'): g[nm] = g[nm].subs(lam, Q(2) / k)
dd, kk, bb, BB = g['solve_fp'](pub)
B2 = (777 - mm.sqrt(1081)) / 852; D2 = (61 - mm.sqrt(1081)) / 66; b2 = -mm.mpf(11) / 71 + 595 * mm.sqrt(1081) / 76751
print("fixed point B = %s\nB_2         = %s\ndiff %.3e ; delta diff %.3e ; b diff %.3e" % (mm.nstr(BB, 32), mm.nstr(B2, 32), float(BB - B2), float(dd - D2), float(bb - b2)))
print("1278B^2 - 2331B + 1061 at the fixed point: %.3e" % float(1278 * BB**2 - 2331 * BB + 1061))
