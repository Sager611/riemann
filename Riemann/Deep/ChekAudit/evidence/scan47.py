# Independent float64 scan of the lattice symbol at c=47 (odd nodes threshold 1, even nodes threshold 0),
# phases k*beta_n mod 1 by exact uint64 fixed-point products (error < 1e-10 cycles), safety margin 1e-7.
# Compares the flagged nodes with Chekmachev's retained node lists. Usage: scan47.py KLO KHI
import sys, csv, time
import numpy as np, mpmath as mp
mp.mp.dps = 40
c = 47; L = mp.log(c)
pp = []
for n in range(2, c):
    f = [p for p in range(2, n+1) if n % p == 0 and all(p % q for q in range(2, p))]
    if len(f) == 1:
        pp.append((n, mp.log(f[0])/mp.sqrt(n)))
B = np.array([int(mp.nint(mp.log(n)/(2*L)*mp.mpf(2)**64)) for n, _ in pp], dtype=np.uint64)
W = np.array([float(w) for _, w in pp]); sig2 = 2*W.sum()
Lf = float(L); lnpi = float(mp.log(mp.pi)); TWO64 = float(2**64)
an = np.array([2*n + 0.5 for n in range(6)]); can = np.array([float(mp.mpf(c)**(-(2*n + mp.mpf(1)/2))) for n in range(6)])
def mu(k):
    k = k.astype(np.uint64)
    r = np.pi*k.astype(np.float64)/Lf
    w = 0.25 + 0.5j*r
    psi = np.log(w) - 1/(2*w) - 1/(12*w**2) + 1/(120*w**4) - 1/(252*w**6) + 1/(240*w**8)
    hp = psi.real - lnpi
    sgn = np.where(k % 2 == 0, 1.0, -1.0)
    T = -2*sgn*sum(can[i]*an[i]/(an[i]**2 + r*r) for i in range(6))
    s = np.zeros_like(r)
    for i in range(len(B)):
        ph = (k*B[i]).astype(np.float64)/TWO64
        s += W[i]*np.cos(2*np.pi*ph)
    return hp - T - 2*s
K2 = 3683; KLO = int(float(sys.argv[1])); KHI = int(float(sys.argv[2])); marg = 1e-7
t0 = time.time()
flo, fle = [], []
CH = 2_000_000
start = max(KLO, K2 + 1)
k0 = start + (1 - start % 2)   # first odd
while k0 <= KHI:
    k = np.arange(k0, min(k0 + 2*CH, KHI + 1), 2, dtype=np.int64)
    m = mu(k); sel = k[m < 1 + marg]; flo.extend(sel.tolist())
    k0 += 2*CH
j0 = max(KLO//2, K2 + 1)
while 2*j0 <= KHI:
    j = np.arange(j0, min(j0 + CH, KHI//2 + 1), dtype=np.int64)
    m = mu(2*j); sel = j[m < 0 + marg]; fle.extend(sel.tolist())
    j0 += CH
t1 = time.time()
base = 'ex/truncated-weil-form-certificates-en/certificates/c47_even_eps/'
odd = {int(r['k']): float(r['mu_mid_40']) for r in csv.DictReader(open(base+'nodes_odd.csv'))}
ev = {int(r['j']): float(r['mu_mid_40']) for r in csv.DictReader(open(base+'nodes_even.csv'))}
their_o = {k for k, v in odd.items() if K2 < k and KLO <= k <= KHI and v < 1}
their_e = {j for j, v in ev.items() if j > K2 and KLO <= 2*j <= KHI and v < 0}
mine_o = set(flo); mine_e = set(fle)
miss_o = sorted(k for k in mine_o if k not in odd)          # flagged by me, not retained by them
miss_e = sorted(j for j in mine_e if j not in ev)
extra_o = sorted(their_o - mine_o); extra_e = sorted(their_e - mine_e)
print(f"c=47 scan k in [{KLO},{KHI}]: odd nodes {(KHI-K2)//2}, even nodes {max(0, KHI//2-K2)}; {t1-t0:.1f}s")
print(f"  odd flagged (mu<1+1e-7): {len(mine_o)}; their retained with mu<1 in range: {len(their_o)}; flagged-but-not-retained: {len(miss_o)} {miss_o[:10]}; theirs-not-flagged: {len(extra_o)} {extra_o[:10]}")
print(f"  even flagged (mu<0+1e-7): {len(mine_e)}; their retained with mu<0 in range: {len(their_e)}; flagged-but-not-retained: {len(miss_e)} {miss_e[:10]}; theirs-not-flagged: {len(extra_e)} {extra_e[:10]}")
# value agreement on their retained nodes in range
ko = np.array(sorted(k for k in odd if KLO <= k <= KHI and k > K2), dtype=np.int64)
if len(ko):
    d = np.abs(mu(ko) - np.array([odd[k] for k in ko])); print(f"  max |mu_float - mu_theirs| on {len(ko)} retained odd nodes: {d.max():.2e}")
je = np.array(sorted(j for j in ev if KLO <= 2*j <= KHI and j > K2), dtype=np.int64)
if len(je):
    d = np.abs(mu(2*je) - np.array([ev[j] for j in je])); print(f"  max |mu_float - mu_theirs| on {len(je)} retained even nodes: {d.max():.2e}")
