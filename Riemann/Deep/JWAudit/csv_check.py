# JWAudit: independent check of the committed per-cell witness files (coverage of 0 <= eps <= 1/L*, and the maxima quoted).
import csv
from fractions import Fraction as Fr
for fn, Ls, Kpage in (('cells_main.csv', 80, None), ('cells_self.csv', 24, None)):
    rows = list(csv.DictReader(open('/Users/adriansager/playground/riemann/Riemann/Deep/SavvaK/' + fn)))
    iv = sorted((Fr(r['eps_lo']), Fr(r['eps_hi'])) for r in rows)
    gaps = [(a, b) for (x, a), (b, y) in zip(iv, iv[1:]) if b > a]
    cov = iv[0][0] <= 0 and max(h for _, h in iv) >= Fr(1, Ls)
    mx = {k: max(float(r[k]) for r in rows) for k in ('C6_up', 'Bn_up', 'Cfin_up', 'radius_budget_up', 'Cloc_up')}
    arg = max(rows, key=lambda r: float(r['Cfin_up']))
    print("%s: %d rows (%d bisected), covers [0, 1/%d]: %s, gaps: %d; max %s; max-Cfin cell %s eps [%s, %s] C6 %s B/n %s" % (
        fn, len(rows), sum(r['cell'] == '-1' for r in rows), Ls, cov, len(gaps),
        ", ".join("%s %.6g" % (k, v) for k, v in mx.items()), arg['cell'], arg['eps_lo'], arg['eps_hi'], arg['C6_up'], arg['Bn_up']))
    print("   all radius budgets < 1: %s; all Cfin <= %s" % (mx['radius_budget_up'] < 1, mx['Cfin_up']))
