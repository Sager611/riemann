# Sensitivity of Table 1's B1 to the fourth-moment inputs of Cor 2.4 (M1(1e10) = 34.842, M2(1e10) = 1.3897e16, M2(1e10,2e10) <= 9.5079e15).
import mpmath as mp, table1_repro as t
rows = [t.TAB[i] for i in (0, 8, 11, 12)]
base = {r[0]: t.row_consts(r, mp.mpf('0.5501'))['B1'] for r in rows}
print('B1 at b4 = 0.5501; factor f applied to M1(1e10) and to both M2 constants')
for f in ['2', '5', '10']:
    t.M1_10, t.M2_10, t.MM2_10 = mp.mpf('34.842')*mp.mpf(f), mp.mpf('1.3897e16')*mp.mpf(f), mp.mpf('9.5079e15')*mp.mpf(f)
    out = []
    for r in rows:
        b = t.row_consts(r, mp.mpf('0.5501'))['B1']; out.append('%s-%s: %.3f (x%.4f)' % (r[0], r[1], b, b/base[r[0]]))
    print('f = %-3s ' % f + ' | '.join(out))
