# assemble.py: check the production sieve logs and report coverage.
# For each chunk log runs/p_<name>.log (chunks file lines "name lo hi"): the DONE line must report the chunk's N, LO, HI,
# lanes = certF + cert1 + t3, t4 = 0 (else the t4 file is listed), bad = 0, no FAIL / MISMATCH line; the lane count
# must equal the number of d = 3, 7, 11, 15, 4, 8 mod 16 in [LO, HI] computed here independently. Then it reports
# the covered interval [lo, top] made of consecutive DONE chunks starting at the given lower end, and for each
# threshold N' (thresholds file) whether e^{T_N'} <= top.
# Usage: python3 assemble.py CHUNKFILE N D_LOW THRESHOLDS
import sys, re, os
chunkf, N, DLOW, thrf = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
S = {3, 7, 11, 15, 4, 8}
def cnt(x):  # number of 0 <= d <= x with d mod 16 in S
    if x < 0: return 0
    q, r = divmod(x, 16)
    return q * len(S) + sum(1 for s in S if s <= r)
rows = []
for l in open(chunkf):
    a = l.split()
    if len(a) != 3: continue
    rows.append((a[0], int(a[1]), int(a[2])))
done = {}
tot = dict(lanes=0, certF=0, cert1=0, t3=0, t3pass=0, nonfund=0, t4=0, tested=0, bad=0, secs=0.0)
print('# chunk  [lo, hi]  lanes  certF  cert1  t3 (t3pass nonfund t4)  seconds  ns/lane  count_exact')
for name, lo, hi in rows:
    p = 'runs/p_%s.log' % name
    if not os.path.exists(p): continue
    txt = open(p).read()
    m = re.search(r'^DONE N=(\d+) \[(\d+),(\d+)\]: lanes (\d+) = certF (\d+) \+ cert1 (\d+) \+ t3 (\d+) \(t3pass (\d+), nonfund (\d+), t4 (\d+)\); tested (\d+), tier3 rechecked (\d+), bad (\d+); ([\d.]+) s', txt, re.M)
    if not m: print(name, 'NOT DONE'); continue
    v = [int(x) for x in m.groups()[:13]] + [float(m.group(14))]
    n_, lo_, hi_, lanes, cF, c1, t3, t3p, nf, t4, tested, rech, bad, secs = v
    ok = (n_ == N and lo_ == lo and hi_ == hi and lanes == cF + c1 + t3 and t3 == t3p + nf + t4 and bad == 0
          and 'FAIL' not in txt and 'MISMATCH' not in txt and lanes == cnt(hi) - cnt(lo - 1))
    if t4: print('  t4 entries in runs/p_%s.t4: %d (decide with t4.gp)' % (name, t4))
    print('%s [%d, %d] %d %d %d %d (%d %d %d) %.0f %.2f %s' % (name, lo, hi, lanes, cF, c1, t3, t3p, nf, t4, secs, secs * 1e9 / lanes, 'yes' if ok else 'NO'))
    if ok: done[lo] = (hi, t4)
    for k, x in zip(['lanes', 'certF', 'cert1', 't3', 't3pass', 'nonfund', 't4', 'tested', 'bad', 'secs'], [lanes, cF, c1, t3, t3p, nf, t4, tested, bad, secs]):
        tot[k] += x
print('# totals', tot)
top = DLOW - 1; t4tot = 0
while top + 1 in done:
    hi, t4 = done[top + 1]; t4tot += t4; top = hi
print('# covered by consecutive DONE chunks: [%d, %d], t4 entries %d' % (DLOW, top, t4tot))
for l in open(thrf):
    if l.startswith('#'): continue
    a = l.split()
    if len(a) < 3: continue
    Np, e0, up = int(a[0]), a[1], float(a[2])
    print('# N=%d: threshold e^%s <= %.4g: %s' % (Np, e0, up, 'COVERED' if up <= top else 'not covered'))
