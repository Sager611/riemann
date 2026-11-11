# table.py -- REPORT table: record row, published bounds combined, routine baseline (T' = 1000), new (7/8), new (B')
import json, math, sys
sys.argv = [sys.argv[0]]
exec(open('compare.py').read().split("cols = []")[0].replace("files = [f for f", "_x = [f for f"))
def load(f):
    try: return {round(r['theta'], 4): (r['log10c'], r['c']) for r in json.load(open(f))}
    except FileNotFoundError: return {}
N78, NBP, RT = load('logs/cert_ch_results.json'), load('logs/cert_ch_results_bp.json'), load('logs/cert_ch_results_bp_T1000.json')
def fmt(c):
    m, e = c.split('e'); m = int(m); e = int(e)
    while m >= 1000: m = (m + 9)//10; e += 1
    k = len(str(m)) - 1; return "%s.%s·10^%d" % (str(m)[0], str(m)[1:], e + k)
print("| θ | record row | published combined | routine (no RH strip) | new, 7/8 certificates | new, + B′ certificates | gain vs record | gain vs routine |")
print("|---|---|---|---|---|---|---|---|")
for th in sorted(NBP, reverse=True):
    re = rowenv(th); b = base(th); n = NBP[th]
    rt = RT.get(th)
    print("| %.4f | %s | 10^%.2f | %s | %s | %s | %s | %s |" % (th, ("10^%.2f" % re) if re is not None else "-", b,
          fmt(rt[1]) if rt else "-", fmt(N78[th][1]) if th in N78 else "-", "**" + fmt(n[1]) + "**",
          ("%.3gx" % 10**(re - n[0])) if re is not None else "-", ("%.3gx" % 10**(rt[0] - n[0])) if rt else "-"))
