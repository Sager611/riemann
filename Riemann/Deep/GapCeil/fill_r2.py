import json, glob
rows = []
for f in glob.glob('logs/certR_r2_*.json'):
    j = json.load(open(f))
    k = f[len('logs/certR_r2_'):-5]
    try:
        mxCd = max(float(c['C'].strip('[').split()[0]) / float(c['deltaD'].strip('[').split()[0]) for c in j['certs'])
    except Exception:
        mxCd = float('nan')
    rows.append((float(k), k, j['float_lam_over_k2'], len(j['certs']), j['ALL_OK'], len(j['cells']), j['D'], j['prec'], j['peak_MB'], j['seconds'], mxCd))
rows.sort()
for r in rows:
    print('| %s | %.7f | %s: %d operators (%d v-cells), D = %d, %d bits, max C/δ = %.2f, %.0f MB, %.0f s |' % (r[1], r[2], 'yes' if r[4] else 'NO', r[3], r[5], r[6], r[7], r[10], r[8], r[9]))
