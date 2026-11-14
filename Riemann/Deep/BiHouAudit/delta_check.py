# Close the reviewer's point "a zero could sit at ordinate T0": Arb counts around T0 and an mpmath check of gamma_{N+1}.
from flint import arb, ctx
from mpmath import mp, mpf, siegelz
ctx.prec = 80
for T in ('1299999.53', '1299999.54', '1299999.8171', '1300000', '1300000.1828', '1300000.1829'):
    print('Arb N(%s) = %s' % (T, arb(T).zeta_nzeros().str(10)))
mp.dps = 50
g = mpf('1300000.18286524889819717690167211334359833')   # Gallego row n = 2325573
t = g
for _ in range(2): t = t - siegelz(t) / siegelz(t, derivative=1)
print('mpmath: zero at', mp.nstr(t, 30), ' |t - data| =', mp.nstr(abs(t - g), 3), ' delta = min(T0 - g_N, g_{N+1} - T0) =', mp.nstr(min(mpf(1300000) - mpf('1299999.53632749461783640283430027729230148'), t - 1300000), 15))
