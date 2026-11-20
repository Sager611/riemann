# ceilaudit6: summarise delta = L^2 (d_N^2 - B) and gapL = L (d_N^2 - B) for zeta, chi_-4, DH from run2.log; windowed means; fits.
import re, sys, numpy as np
rows = {'Z': [], 'CHI': [], 'DH': []}
pat = re.compile(r'DELTA (?:(CHI|DH) )?N= *(\d+) L=([\d.]+) d2=([\d.e+-]+) B[^=]*=([\d.e+-]+)')
for l in open(sys.argv[1] if len(sys.argv) > 1 else 'logs/run2.log'):
    m = pat.match(l)
    if m: rows[m.group(1) or 'Z'].append((int(m.group(2)), float(m.group(3)), float(m.group(4)), float(m.group(5))))
for k, v in rows.items():
    a = np.array(v); N, L, d2, B = a.T; dl = L * L * (d2 - B); g = L * (d2 - B)
    print(f'--- {k}: N {int(N.min())}..{int(N.max())} ({len(N)} pts)')
    for lo, hi in [(5, 50), (50, 200), (200, 600), (600, 1000), (1000, 1500), (1500, 2000), (2000, 2500), (2500, 3000), (3000, 3500), (3500, 4001)]:
        s = (N >= lo) & (N < hi)
        if s.sum(): print(f'  [{lo},{hi}): n={s.sum():2d} delta {dl[s].mean():.4f} sd {dl[s].std():.4f} gapL {g[s].mean():.5f} d2L {np.mean(d2[s]*L[s]):.5f} BL {np.mean(B[s]*L[s]):.5f}')
    for lo in (600, 1000):
        s = N >= lo
        if s.sum() > 4:
            A1 = np.vstack([np.ones(s.sum()), L[s]]).T; c1 = np.linalg.lstsq(A1, dl[s], rcond=None)[0]
            A2 = np.vstack([np.ones(s.sum()), 1 / L[s]]).T; c2 = np.linalg.lstsq(A2, dl[s], rcond=None)[0]
            print(f'  fit N>={lo}: delta = {c1[0]:.4f} + {c1[1]:.4f} L (sd {(dl[s]-A1@c1).std():.5f}) | {c2[0]:.4f} + {c2[1]:.4f}/L (sd {(dl[s]-A2@c2).std():.5f})')
