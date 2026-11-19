import json
N = json.load(open('../IntCheb/CERTIFICATE_LOWER_new.json')); P = json.load(open('../../../.scratch/intchebaudit/tex/anc/CERTIFICATE_integer_Chebyshev_LOWER_2026-08-09.json'))
pe, pm = P['cell_edges'], P['masses']; ne = N['cell_edges']
# paper support intervals = maximal runs of positive-mass cells
runs = []; k = 0
while k < len(pm):
    if pm[k] > 0:
        j = k
        while j + 1 < len(pm) and pm[j + 1] > 0: j += 1
        runs.append((pe[k], pe[j + 1], j - k + 1)); k = j + 1
    else: k += 1
half = [r for r in runs if r[1] <= 0.5 + 1e-15]
print('paper support intervals: %d total, %d inside [0,1/2], straddling 1/2: %s' % (len(runs), len(half), [r for r in runs if r[0] < 0.5 < r[1]]))
ratios = []
for s, t, n in half:
    m = sum(1 for i in range(len(ne) - 1) if ne[i] >= s - 1e-12 and ne[i + 1] <= t + 1e-12)
    ratios.append((round(s, 6), round(t, 6), n, m))
print('per interval (start, end, paper cells, our cells inside):')
for r in ratios: print('  ', r)
print('exactly doubled:', sum(1 for r in ratios if r[3] == 2 * r[2]), 'of', len(ratios))
