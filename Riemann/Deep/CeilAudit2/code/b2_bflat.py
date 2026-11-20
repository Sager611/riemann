# B2 float test: Liu's B-flat stage map (BFlat/bflat_sens.py definitions, executed read-only) with the M2 prime cost 6*kappa replaced by
# the constant costs of 2611.00082's two moment inductions (n1+n2+3z <= M, then n1+n2+2z <= M), everything else at NF3's constants.
import os, sys, io, contextlib
src = open(os.path.join(os.path.dirname(__file__), '..', '..', 'BFlat', 'bflat_sens.py')).read()
cut = src.index("BF = mm.mpf")
g = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src[:cut], 'bflat_sens_defs', 'exec'), g)
g['BF'] = g['mm'].mpf('0.8748866470399641309446140')
Q, pub, lam, k = g['Q'], g['pub'], g['lam'], g['k']
orig = {nm: g[nm] for nm in ('E', 'TS', 'RS', 'coef_th')}
g['report']("published: M2 cost 6*kappa", pub)
for c in (Q(9, 2), Q(4), Q(3), Q(5, 2), Q(2)):
    for nm in orig: g[nm] = orig[nm].subs(lam, c / k)
    try:
        g['report']("M2 cost constant %s" % c, pub)
    except Exception as e:
        print("M2 cost constant %s: %s" % (c, e))
    sys.stdout.flush()

# lever map below B_2 in the same system (M2 cost < 2), and the other NF3 constants at M2 cost 2
mu, al = g['mu'], g['al']
for c in (Q(7, 4), Q(3, 2), Q(1)):
    for nm in orig: g[nm] = orig[nm].subs(lam, c / k)
    try: g['report']("M2 cost constant %s" % c, pub, guess=(0.41, 0.745, 0.11))
    except Exception as e: print("M2 cost constant %s: %s" % (c, e))
    sys.stdout.flush()
for nm in orig: g[nm] = orig[nm].subs(lam, Q(2) / k)
for nm_, prm in (("M2 cost 2, M1 cost 2 -> 19/10", {**pub, mu: Q(19, 10)}), ("M2 cost 2, amplification 5/6 -> 4/5", {**pub, al: Q(4, 5)})):
    try: g['report'](nm_, prm, guess=(0.42, 0.746, 0.11))
    except Exception as e: print(nm_, e)
    sys.stdout.flush()
