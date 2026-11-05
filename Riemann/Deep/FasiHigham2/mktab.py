
# build Good-lemma tables at nu0, save JSON: per family, list of [r_lo, r_hi, e, bh, lw0, kind]
import good, goodsmall, json, sys, time
def fup(x):
    v = float(x.upper()); return v * (1 + 1e-9) + 1e-300 if v > 0 else v * (1 - 1e-9)
from flint import arb
nu0 = arb(sys.argv[1]); out = sys.argv[2]
NMIN = 8
fams = [('0', 1, 1), ('1', 3, 3), ('2', 5, 5), ('3', 7, 7), ('4', 9, 9), ('>=5', 0, 11)]
RANGES = {'0': (0.15, 1.6), '1': (0.03, 1.6), '2': (0.05, 1.6), '3': (0.07, 1.6), '4': (0.08, 1.6), '>=5': (0.1, 1.6)}
W = 0.03
tab = {}
t0 = time.time()
# l = 0 small r
rows = []
edges = [0, 1e-3, 3e-3, 0.01, 0.02, 0.03, 0.045, 0.06, 0.08, 0.1, 0.125, 0.15]
for i in range(len(edges) - 1):
    res = goodsmall.good_bound_l0small(arb(edges[i]), arb(edges[i + 1]), NMIN, nu0)
    # beta/nu <= t(r_hi) * bh_coef + (rmax^2 + lam c2)/nu0  (stored as upper bound)
    if res:
        rows.append([edges[i], edges[i + 1], fup(res['e']), fup(res['bnu']), None, 'small'])
    else:
        rows.append([edges[i], edges[i + 1], None, None, None, 'small'])
tab['0small'] = rows
print('0small done', ['%.3f' % r[2] if r[2] is not None else None for r in rows], '%.0fs' % (time.time() - t0), flush=True)
for (name, a, adiv) in fams:
    if name == '>=5':
        al_lo, al_hi = arb(0), arb(1) / 11
    else:
        al_lo = al_hi = arb(1) / a
    rlo, rhi = RANGES[name]
    r = rlo; rows = []
    while r < rhi - 1e-12:
        r2 = min(rhi, r + W * max(r, 0.05))
        # do not straddle yJ badly: split at 0.85
        if r < 0.85 < r2: r2 = 0.85
        res = good.good_bound(al_lo, al_hi, arb(r), arb(r2), nu0)
        if res is None:
            rows.append([r, r2, None, None, None, ''])
        else:
            rows.append([r, r2, fup(res['e']), fup(res['bh']), fup(res['lw0']), 'c' if res['custom'] else 'u'])
        r = r2
    tab[name] = rows
    bad = [x[:2] for x in rows if x[2] is None or x[2] > 0.5]
    print(name, 'boxes', len(rows), 'max e %.4f' % max(x[2] for x in rows if x[2] is not None), 'bad', bad[:6], '%.0fs' % (time.time() - t0), flush=True)
json.dump(tab, open(out, 'w'))
