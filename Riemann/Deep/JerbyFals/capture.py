# capture rate of the resonance prefilter: L-window CANDs inside the R-scan range that the R-scan also found
import re, sys
def cands(fn):
    out = {}
    for l in open(fn):
        if l.startswith('CAND '):
            n = int(re.search(r' n=(\d+)', l).group(1)); mu = float(re.search(r'mu=(\S+)', l).group(1)); out[n] = mu
    return out
R = {}
for fn in sys.argv[2:]: R.update(cands(fn))
L = cands(sys.argv[1])
lo = min(R); hi = max(R)
inr = {n: mu for n, mu in L.items() if lo <= n <= hi}
for thr in (4.6, 5.0, 6.0):
    s = [n for n, mu in inr.items() if mu < thr]; f = [n for n in s if n in R]
    print(f'L-window CANDs with mu<{thr} inside R range: {len(s)}, found by R scan: {len(f)}')
