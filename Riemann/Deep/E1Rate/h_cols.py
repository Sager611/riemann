# e1rate: h(theta;kappa) itself (old nu vs new nu) next to Wang's c(theta) and (c+2kappa)/3; kappa from h_table.log (mpmath values).
import mpmath as mp, re
mp.mp.dps = 25
K = {}
for line in open('h_table.log'):
    m = re.match(r'nu=([0-9.]+) .*kappa_mp\(sinh P\)=([0-9.]+)', line)
    if m: K[round(float(m.group(1)), 6)] = mp.mpf(m.group(2))
s2 = mp.sqrt(2)
c = lambda th: 2-th/2-mp.cot(th/s2)/s2
h = lambda th, k: (3+k-mp.sqrt((1-k)*(9-k-8*c(th))))/4
print("theta | nu_old kappa_old h_old | nu_new kappa_new h_new | c(theta) | (c+2k_new)/3 | max (LB)")
for th in [0.51, 0.6, 0.75, 0.82, 0.85, 0.9, 0.95, 0.99]:
    nE = (3*th-1)/4; nB = min(0.5, 17/33*(2*th-1)); no = max(nE, nB); nn = min(0.5, max(nE, 4/7*(2*th-1)))
    t = mp.mpf(th); ko = K[round(no, 6)]; kn = K[round(nn, 6)]
    print(f"{th} | {no:.6f} {mp.nstr(ko,9)} {mp.nstr(h(t,ko),9)} | {nn:.6f} {mp.nstr(kn,9)} {mp.nstr(h(t,kn),9)} | {mp.nstr(c(t),9)} | {mp.nstr((c(t)+2*kn)/3,9)} | {mp.nstr(max(0,c(t),(c(t)+2*kn)/3,h(t,kn)),9)}")
