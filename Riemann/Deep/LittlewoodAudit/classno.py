# Falsification (i): exact class numbers h(-d) of all fundamental -d with 4 < d <= N by counting reduced forms (numpy
# strided slices), against Corollary 2: h(-d) >= 0.0350 sqrt(d)/max(3, loglog d); and L(1,chi_-d) = pi h/sqrt(d) against
# Theorem 1(a): >= 0.11/max(3, loglog d).  Run: uv run --no-project --with numpy python3 classno.py N
import sys, math, time
import numpy as np
N = int(sys.argv[1]); t0 = time.time()
cnt = np.zeros(N + 1, dtype=np.int32)
a = 1
while 3*a*a <= N:
    for b in range(-a + 1, a + 1):
        cmin = a if b >= 0 else a + 1                     # c = a needs b >= 0; |b| = a needs b = a
        cmin = max(cmin, (b*b + 4*a) // (4*a))            # d = 4ac - b^2 > 0 (automatic for c >= a)
        d0 = 4*a*cmin - b*b
        if d0 <= N: cnt[d0::4*a] += 1
    a += 1
# fundamental discriminants -d: d = 3 mod 4 squarefree, or d = 4m, m = 1, 2 mod 4 squarefree
sf = np.ones(N + 1, dtype=bool); p = 2
while p*p <= N:
    sf[p*p::p*p] = False; p += 1
d = np.arange(N + 1, dtype=np.int32)
fund = ((d % 4 == 3) & sf) | ((d % 4 == 0) & np.isin((d//4) % 4, [1, 2]) & sf[np.minimum(d//4, N)])
fund[:5] = False
D = d[fund]; H = cnt[fund].astype(np.float64)
known = {7: 1, 8: 1, 11: 1, 15: 2, 20: 2, 23: 3, 24: 2, 47: 5, 71: 7, 84: 4, 163: 1, 199: 9, 1555: 4, 5923: 7, 2383747: 98}
bad = [(k, int(cnt[k]), v) for k, v in known.items() if k <= N and cnt[k] != v]
print('N', N, 'fundamental d:', len(D), 'known-value check mismatches:', bad, f't={time.time()-t0:.1f}s')
sq = np.sqrt(D.astype(np.float64)); den = np.maximum(3.0, np.log(np.log(D.astype(np.float64))))
ratio = H*den/sq                                   # Corollary 2 needs ratio >= 0.0350
i = int(np.argmin(ratio)); print(f'min h*max(3,loglog d)/sqrt(d) = {ratio[i]:.5f} at d = {D[i]} (h = {int(H[i])}); claim 0.0350: {"OK" if ratio[i] >= 0.0350 else "VIOLATED"}')
L1 = math.pi*H/sq
j = int(np.argmin(L1*den)); print(f'min L(1,chi_-d)*max(3,loglog d) = {L1[j]*den[j]:.5f} at d = {D[j]}; claim 0.11: {"OK" if L1[j]*den[j] >= 0.11 else "VIOLATED"}')
for lo, hi in ((5, 10**3), (10**3, 10**4), (10**4, 10**5), (10**5, 10**6), (10**6, 10**7)):
    m = (D > lo) & (D <= hi)
    if m.any():
        k = int(np.argmin(L1[m])); print(f'  d in ({lo},{hi}]: min L(1) = {L1[m][k]:.4f} at d = {D[m][k]}, min h = {int(H[m].min())}, max d with h <= 100: {int(D[m][H[m] <= 100].max()) if (H[m] <= 100).any() else None}')
print('largest d <= N with h(-d) <= 100:', int(D[H <= 100].max()), ' (Watkins: 2383747)' , f't={time.time()-t0:.1f}s')
