# nbaudit c2: analyse the exact certified intervals of the lean N <= 1280 rerun (/tmp/nbaudit/loh_V1_lean.txt):
# widths, sign margins (mpmath, 50 digits), shipped-table rows, containment of the b1 mpmath values (N <= 100), agreement with b3 double sweep.
# usage: uv run --no-project --with mpmath --with numpy python3 c2_check.py
import mpmath as mp, numpy as np
from fractions import Fraction as Fr
mp.mp.dps = 50
C = 2 + mp.euler - mp.log(4 * mp.pi)
I = {}
for line in open('/tmp/nbaudit/loh_V1_lean.txt'):
    N, a, b, c, d = map(int, line.split()); I[N] = (Fr(a, b), Fr(c, d))
NM = max(I); assert sorted(I) == list(range(1, NM + 1)), 'gap in N'
W = {N: float(I[N][1] - I[N][0]) for N in I}
big = [N for N in I if W[N] >= 1e-23]
print(f'N covered: 1..{NM}, no gaps; max width {max(W.values()):.3e} at N={max(W, key=W.get)}; #N with width >= 1e-23: {len(big)} (first {big[:8]}, last {big[-3:]})')
print('width quantiles (50/90/99%):', [f'{np.quantile(list(W.values()), q):.2e}' for q in (0.5, 0.9, 0.99)])
marg = {}
for N in range(2, NM + 1):
    lo = mp.mpf(I[N][0].numerator) / I[N][0].denominator; hi = mp.mpf(I[N][1].numerator) / I[N][1].denominator
    fl, fh = lo * mp.log(N) - C, hi * mp.log(N) - C
    assert (fl > 0) == (fh > 0), f'undecided N={N}'
    marg[N] = fl if fl > 0 else fh   # signed value nearest to zero
s = [1 if marg[N] > 0 else -1 for N in range(2, NM + 1)]
ch = [N for N in range(3, NM + 1) if s[N - 2] != s[N - 3]]
print('sign changes (count, at N):', len(ch), ch)
print('tightest margins:', sorted(((abs(float(marg[N])), N) for N in marg))[:4])
for lo_, hi_ in ((2, 640), (641, NM)):
    print(f'below C on [{lo_},{hi_}]: {sum(1 for N in range(lo_, hi_ + 1) if marg[N] < 0)} of {hi_ - lo_ + 1}; min |margin| {min(abs(float(marg[N])) for N in range(lo_, hi_ + 1)):.3e}')
print('all N in [584, NM] below C:', all(marg[N] < 0 for N in range(584, NM + 1)), '; first below C:', next(N for N in range(2, NM + 1) if marg[N] < 0))
tab = {40: '0.012715197195337', 41: '0.012326491948694', 48: '0.011873207147379', 49: '0.011872614669304', 80: '0.010983481552538', 160: '0.0093455910085642',
       205: '0.0086793734308363', 206: '0.0086522007981348', 320: '0.0079537374311735', 411: '0.0076744466976446', 524: '0.0073767301615984',
       583: '0.0072572423228183', 584: '0.0072490880608171', 640: '0.0070170593529348', 695: '0.0067526049812340', 1280: '0.0063611026375655'}
bad = [N for N, v in tab.items() if abs(mp.mpf(v) - mp.mpf(I[N][0].numerator) / I[N][0].denominator) > 10 ** -(len(v.split('.')[1]) ) ]
print('REPORT table d_N^2 rows agree with certified intervals to the printed digits:', not bad, bad)
L640 = {}
for line in open('/tmp/nbaudit/loh_V1.txt'):
    N, a, b, c, d = map(int, line.split()); L640[N] = (Fr(a, b), Fr(c, d))
ov = all(max(I[N][0], L640[N][0]) <= min(I[N][1], L640[N][1]) for N in L640)
print('NMAX=1280 intervals intersect the NMAX=640 rerun intervals for every N <= 640:', ov)
d1 = np.load('/tmp/nbaudit/d1_b3.npy')
dev = max(abs(d1[N] - float((I[N][0] + I[N][1]) / 2)) for N in range(1, NM + 1))
print(f'max |b3 double - certified mid| over N <= {NM}: {dev:.3e}')
mono = all(I[N + 1][0] <= I[N][1] for N in range(1, NM))
print('interval consistency with monotonicity (lo_{N+1} <= hi_N) for all N:', mono)
