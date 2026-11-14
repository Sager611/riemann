# Validate Gallego's rows (gamma, c_rho = 1/(rho zeta'(rho))) with mpmath at 50 digits (independent of Arb/FLINT).
# usage: python3 val_mpmath.py DATA_DIR SAMPLE_JSON N_EXTRA OUT_LOG
import sys, json, random, time, os
from mpmath import mp, mpf, mpc, siegelz, zeta, fabs
from crho_format import CrhoFile
mp.dps = 50
DATA, SAMPLE, NEXTRA, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
FILES = ['cc_00000001.crho', 'cc_00240000.crho', 'cc_00305000.crho', 'cc_00371494.crho', 'cc_00742987.crho',
         'cc_01114480.crho', 'cc_01485973.crho', 'cc_01857466.crho', 'cc_02228959.crho']
rows = [r for r in json.load(open(SAMPLE)) if 're' in r]
have = {r['n'] for r in rows}
# extra rows: top-weight zeros n=1..15, plus NEXTRA uniform random indices in [1, 2325572] (seeded)
rng = random.Random(777)
want = set(range(1, 16)) | {rng.randint(1, 2325572) for _ in range(NEXTRA)}
want -= have
for fn in FILES:
    cf = CrhoFile(os.path.join(DATA, fn))
    for n in sorted(w for w in want if cf.n0 <= w < cf.n0 + cf.count):
        i = n - cf.n0; rows.append({'n': n, 'gamma': cf.txt('gamma', i), 're': cf.txt('re', i), 'im': cf.txt('im', i),
                                    'radio_c': float(cf.col['radio_c'][i]), 'radio_gamma': cf.radius_gamma_txt}); want.discard(n)
    del cf
rows.sort(key=lambda r: r['n'])
out = open(OUT, 'w'); t0 = time.time()
worst = {'dg_over_rad': 0, 'dg_abs': 0, 'dc_over_rad': 0, 'dc_rel': 0, 'step2': 0}; fails = 0
for r in rows:
    g0 = mpf(r['gamma']); t = g0
    steps = []
    for _ in range(2):
        st = siegelz(t) / siegelz(t, derivative=1); t = t - st; steps.append(abs(st))
    dg = fabs(t - g0)
    rho = mpc(mpf(1) / 2, t); c = 1 / (rho * zeta(rho, derivative=1))
    dre, dim = fabs(c.real - mpf(r['re'])), fabs(c.imag - mpf(r['im']))
    dc = max(dre, dim); radg = mpf(r['radio_gamma']); radc = mpf(r['radio_c'])
    ok = dg <= radg and dc <= radc
    fails += (not ok)
    worst['dg_over_rad'] = max(worst['dg_over_rad'], float(dg / radg)); worst['dg_abs'] = max(worst['dg_abs'], float(dg))
    worst['dc_over_rad'] = max(worst['dc_over_rad'], float(dc / radc)); worst['dc_rel'] = max(worst['dc_rel'], float(dc / abs(c)))
    worst['step2'] = max(worst['step2'], float(steps[1]))
    out.write('%d %s dgamma=%s step2=%s |c|=%s dc=%s radc=%s %s\n' % (r['n'], r['gamma'][:20], mp.nstr(dg, 3), mp.nstr(steps[1], 3),
              mp.nstr(abs(c), 6), mp.nstr(dc, 3), mp.nstr(radc, 3), 'OK' if ok else 'FAIL')); out.flush()
summary = 'rows %d (sample %d + extra), fails %d, worst %s, %.0fs' % (len(rows), len(have), fails, worst, time.time() - t0)
out.write(summary + '\n'); print(summary)
