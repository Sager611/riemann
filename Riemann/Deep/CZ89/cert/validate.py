
import sys, time, numpy as np
import bnb as Bm
import explore as E
rng = np.random.default_rng(5)
T = np.load("tables.npz")
HT = Bm.HT
# 1) tables vs float evaluation at grid points and random points
m = rng.integers(0, len(T["qlo"]), 20000)
x = m*HT
qf = E.K(x)**2; qpf = 2*E.K(x)*E.Kp(x)
print("grid: q in [qlo,qhi] violations", int(((qf < T["qlo"][m]-1e-14) | (qf > T["qhi"][m]+1e-14)).sum()),
      " q' violations", int(((qpf < T["q1lo"][m]-1e-13) | (qpf > T["q1hi"][m]+1e-13)).sum()))
c = rng.integers(0, len(T["cqmin"]), 200000); u = rng.random(200000)
xr = (c+u)*HT
qr = E.K(xr)**2
print("cells: cqmin > q(x) violations", int((T["cqmin"][c] > qr + 1e-15).sum()), " max(cqmin - q)", float((T["cqmin"][c]-qr).max()))
hh = 1e-4
q2 = (E.K(xr+hh)**2 - 2*E.K(xr)**2 + E.K(xr-hh)**2)/hh**2
print("cells: |q''| > cq2 violations", int((np.abs(q2) > T["cq2"][c] + 1e-5).sum()), " max ratio", float((np.abs(q2)/T["cq2"][c]).max()))
# 2) accepted-box bounds vs sampled S
def S(g): return E.S_exact(g)
nb = 0; viol = 0; worst = -1
for trial in range(4000):
    # random box somewhere relevant
    if trial % 2 == 0:
        ctr = np.array([1.0837,1.8939,1.0414,1.0249,1.0289,1.0256]) + rng.normal(0, 0.05, 6)
    else:
        ctr = rng.uniform(0.0, 3.0, 6)
    wdt = rng.integers(1, 400, 6)
    L = np.floor(ctr/Bm.HB).astype(np.int64)[None,:]; U = L + wdt[None,:]
    if L.sum() >= Bm.SUMB or (L < 0).any(): continue
    U = Bm.tighten(L, U)
    if (U < L).any(): continue
    a1, a2, L1, L2 = Bm.evaluate(L, U)
    lb = max(L1[0], L2[0])
    G = (L + (U-L)*rng.random((3000,6)))*Bm.HB
    G = G[G.sum(axis=1) <= 8.8]
    corners = (L + (U-L)*rng.integers(0,2,(64,6)))*Bm.HB
    corners = corners[corners.sum(axis=1) <= 8.8]
    G = np.concatenate([G, corners])
    if len(G) == 0: continue
    smin = S(G).min()
    nb += 1
    if smin < lb - 1e-13: viol += 1
    worst = max(worst, lb - smin)
print("boxes tested", nb, "bound violations", viol, "max(lb - sampled min)", worst)

