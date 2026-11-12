from mpmath import mp, mpf, cot, sqrt
mp.dps = 30
c = lambda t: 2 - t/2 - cot(t/sqrt(2))/sqrt(2)
h = lambda t, k: (3 + k - sqrt((1-k)*(9 - k - 8*c(t))))/4
k = mpf('0.08447790796974840320324474')   # nu=3/25, R=177/20, (P,Q) of frozen nu=999/10000 set (nu12.log)
for t in ['0.5', '0.51', '0.527', '0.55', '0.6']:
    print('theta', t, 'h =', mp.nstr(h(mpf(t), k), 12), ' (3theta-1)/4 =', mp.nstr((3*mpf(t)-1)/4, 6))
