# analysis.py (explore48): merge first-off-line-zero data; Airy test.
# Birth law tested: Im rho_1 = k T_s + |a'_j| (k T_s / 2)^{1/3}, T_s = 2 pi y, a'_j = zeros of Ai' (-1.0188, -3.2482, -4.8201, -6.1633).
import json, math, glob, os
A = [1.018793, 3.248198, 4.820099, 6.163307, 7.372177]
files = ['logs/t4_val.json', 'logs/t4_h2_0.json', 'logs/t4_rect_0.json', 'logs/t4_rect_7.json', 'logs/t4_rhomb_0.json',
         'logs/t4_redo.json', 'logs/t4_big.json', 'logs/t4_h1.json']
rows = {}
for f in files:
    if not os.path.exists(f): continue
    for r in json.load(open(f)):
        key = (r['label'], round(r['x'], 4), round(r['y'], 4))
        rows[key] = r  # later files (redo, big, h1 refined) override
print('%-20s %6s %8s %8s %s' % ('lattice', 'x', 'y', 'T_s', 'first off-line zero / ratio / k / Airy xi = (t - kT_s)(2/kT_s)^{1/3} vs |a_j\'|'))
stats = []
for key, r in sorted(rows.items(), key=lambda kv: (kv[0][0][:4], kv[0][2])):
    lab, x, y = key; Ts = r['Ts']
    f = r.get('first'); z = f['zero'] if f else None
    if z is None:
        print('%-20s %6.2f %8.4f %8.3f  none up to %.2f T_s (T1/T_s = %s)' % (lab, x, y, Ts, r['Tscan'] / Ts, ('%.3f' % r['T1r']) if r['T1r'] else '-'))
        continue
    t = z[1]; ratio = t / Ts
    k = max(1, int(math.floor(ratio + 0.06)))
    xi = (t - k * Ts) * (2 / (k * Ts)) ** (1 / 3)
    j = min(range(len(A)), key=lambda i: abs(A[i] - xi))
    stats.append((lab, y, ratio, k, xi, j + 1, A[j]))
    print('%-20s %6.2f %8.4f %8.3f  %.6f%+.6fi  %.4f  k=%d  xi=%+.3f  nearest |a_%d|=%.3f  Re-1/2=%.4f (x y = %.2f)' % (
        lab, x, y, Ts, z[0], t, ratio, k, xi, j + 1, A[j], z[0] - 0.5, (z[0] - 0.5) * y))
big = [s for s in stats if s[1] >= 8]
if big:
    print('y >= 8: n=%d, |xi - |a_j|| mean %.3f max %.3f' % (len(big), sum(abs(s[4] - s[6]) for s in big) / len(big), max(abs(s[4] - s[6]) for s in big)))
print('all: n=%d, ratio range %.3f..%.3f for k=1; k values: %s' % (len(stats), min(s[2] for s in stats if s[3] == 1), max(s[2] for s in stats if s[3] == 1), sorted(set(s[3] for s in stats))))
