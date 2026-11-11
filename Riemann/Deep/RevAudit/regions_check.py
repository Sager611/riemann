# Ratio (a1 t^{1/6} log t + a2 t^{1/6} + a3 t^{phi/2} + a4)/(t^{1/6} log t) of the 10 printed regions of 2602.05614 (TeX l. 962-1081),
# with the printed coefficients, on a fine log-grid of [t0, t1] (mpmath 30 digits): does any exceed 0.611?
# TeX: curl -s -L -o /tmp/zr.gz https://arxiv.org/e-print/2602.05614 && gunzip -c /tmp/zr.gz > /tmp/ZetaRevers.tex ; argv[1] overrides the path.
import re, sys, mpmath as mp
mp.mp.dps = 30
s = open(sys.argv[1] if len(sys.argv) > 1 else '/tmp/ZetaRevers.tex', encoding='utf8', errors='replace').read()
blk = s[s.index('Bottleneck region $3.80'):s.index('Proof for the region $8.97')]
pat = re.compile(r't_\{0\}=([\d.]+)\\cdot10\^\{(\d+)\}\$,\s*\$t_\{1\}\s*=\s*([\d.]+)\s*\\cdot\s*10\^\{(\d+)\}\$,\s*\$c=([\d.]+)\$\s*and\s*\$\\phi=([\d.]+)\$.*?\\leq\s*([\d.]+)t\^\{1/6\}\\log\s*t([-+][\d.]+)t\^\{1/6\}([-+][\d.]+)t\^\{\\phi/2\}([-+][\d.]+)', re.S)
rows = pat.findall(blk)
print('regions parsed:', len(rows))
worst_all = mp.mpf(0)
for (a, ea, b, eb, c, phi, a1, a2, a3, a4) in rows:
    t0 = mp.mpf(a) * 10 ** int(ea); t1 = mp.mpf(b) * 10 ** int(eb); ph = mp.mpf(phi)
    A1, A2, A3, A4 = map(mp.mpf, (a1, a2, a3, a4))
    f = lambda t: (A1 * t ** (mp.mpf(1) / 6) * mp.log(t) + A2 * t ** (mp.mpf(1) / 6) + A3 * t ** (ph / 2) + A4) / (t ** (mp.mpf(1) / 6) * mp.log(t))
    grid = [t0 * (t1 / t0) ** (mp.mpf(k) / 4000) for k in range(4001)]
    vals = [f(t) for t in grid]; k = max(range(len(vals)), key=lambda i: vals[i])
    worst_all = max(worst_all, vals[k])
    print('[%s, %s] c=%s phi=%s: ratio(t0)=%s ratio(t1)=%s max=%s at t=%s %s' % (mp.nstr(t0, 4), mp.nstr(t1, 4), c, phi, mp.nstr(vals[0], 8), mp.nstr(vals[-1], 8), mp.nstr(vals[k], 8), mp.nstr(grid[k], 6), 'EXCEEDS 0.611' if vals[k] > mp.mpf('0.611') else 'ok'))
print('max over all regions:', mp.nstr(worst_all, 9))
