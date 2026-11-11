# explore30 qmc2: (a) delta-dependence of the adversary's needed gap (delta strictly inside the B_2 half-plane);
# (b) real zeros near t0 = 15000 (ExploreHelio zeta_zeros.npy): true configuration, one merged off-line pair,
#     and adversarial blocks of K off-line zeros (count-preserving) seen by the c = 47 window.
# usage: uv run --no-project --with numpy python3 qmc2_real.py      (imports qmc1_adversary; chunked, < 100 MB)
import numpy as np, math, sys, time
sys.path.insert(0, '.')
from qmc1_adversary import TwoBump, adversary, Sbounds, memcheck
t0 = time.time()
Sb, DT = Sbounds(3e12 + 1e4); Sallow = Sb - 1.0
PART = sys.argv[1] if len(sys.argv) > 1 else 'ab'
if 'a' in PART:
  print(f'# (a) needed gap vs delta, w = 0.06a, q = 1 (T = 3e12: D = {DT:.4f}, allowed |S| {Sallow:.3f}, allowed gap {2*Sallow/DT:.3f}; t -> inf: 1.304)')
  for d in (0.37, 0.3, 0.2, 0.1, 0.05):
    row = []
    for a in (1.2, 1.4, 1.6, 1.8, 1.925, 2.131, 2.303):
        r = adversary(a, 0.06 * a, DT, d, 1.0)
        ok = r['W'] < 0 and r['Sabs'] <= Sallow
        row.append(f'a={a}: gap {r["gap"]:.3f} W {r["W"]:+.3f}{"*" if ok else ""}{"+" if r["gap"] <= 1.304 and r["W"] < 0 else ""}')
    print(f'  d = {d}: ' + '; '.join(row), flush=True)
  print('  (* admissible at T = 3e12 with |S| <= allowed; + needed gap <= 1.304, admissible at every height)')
memcheck('a')
z = np.load('../ExploreHelio/results/zeta_zeros.npy')
z = np.sort(z[(z > 0)])
T0 = 15000.0; a = 1.925; w = 0.06 * a
UMAX = 60.0 / w                      # |u w| <= 60: 160-point Gauss-Legendre exact to ~1e-90, tail |Phi|^2 < 1e-17
nall = len(z); z = z[np.abs(z - T0) <= UMAX]
print(f'# (b) zeros loaded: {nall}; kept {len(z)} with |gamma - t0| <= {UMAX:.1f}, range {z[0]:.3f} .. {z[-1]:.3f}')
D = math.log(T0 / (2 * math.pi)) / (2 * math.pi)
dd = 0.3
f = TwoBump(a, w)
u = z - T0
base = D * 2 * math.pi * f.norm2
Wtrue = float(np.sum(np.abs(f.F(u))**2)) / base
print(f'  t0 = {T0}, D = {D:.4f}, a = {a}, w = {w:.4f}: true zeros W/(D 2pi||f||^2) = {Wtrue:.4f}')
# one merged pair: two consecutive zeros nearest t0 -> off-line pair at their mean, delta = 0.3; scan phase theta
k = int(np.searchsorted(z, T0))
best = None
for th in np.linspace(0, 2 * math.pi, 73)[:-1]:
    g = TwoBump(a, w, th)
    W0 = np.abs(g.F(u))**2
    for (i, j) in ((k - 1, k), (k, k + 1)):
        um = (u[i] + u[j]) / 2
        Wm = (float(W0.sum()) - W0[i] - W0[j] + 2 * float(g.h(np.array([um]), dd)[0])) / base
        if best is None or Wm < best[0]:
            best = (Wm, th, i, j)
print(f'  one merged pair (delta = 0.3): min over 72 phases W = {best[0]:.4f} at theta = {best[1]:.3f} (true {Wtrue:.4f}): invisible')
# adversarial blocks: replace zeros with |u| <= R by pairs in the negative intervals of h_d (count-preserving per period)
def block(R, th=0.0, d=dd):
    g = TwoBump(a, w, th)
    P = math.pi / g.ap
    grid = np.arange(-R, R, P / 200)
    neg = g.h(grid, d) < 0
    e = np.flatnonzero(np.diff(neg.astype(np.int8)))
    st = [q + 1 for q in e if not neg[q]]; en = [q for q in e if neg[q]]
    if neg[0]: st = [0] + st
    if neg[-1]: en = en + [len(grid) - 1]
    ints = [(grid[s_], grid[e_]) for s_, e_ in zip(st, en) if e_ > s_]
    if len(ints) < 3:
        return None
    lo = ints[0][0] - (ints[1][0] - ints[0][1]) / 2; hi = ints[-1][1] + (ints[-1][0] - ints[-2][1]) / 2
    inside = (u >= lo) & (u < hi)
    bounds = [lo] + [(ints[q][1] + ints[q + 1][0]) / 2 for q in range(len(ints) - 1)] + [hi]
    cum = lambda x: int(np.searchsorted(u, x)) - int(np.searchsorted(u, lo))
    ntot = int(inside.sum())
    placed = 0; offs = []
    for q, (x, y) in enumerate(ints):
        tgt = cum(bounds[q + 1]) if q < len(ints) - 1 else ntot
        npair = max(0, (tgt - placed) // 2)
        if q == len(ints) - 1:
            npair = (ntot - placed) // 2
        placed += 2 * npair
        c = (x + y) / 2; r = (y - x) / 2
        if npair == 1: offs.append(c)
        elif npair > 1: offs.extend((c + r * np.linspace(-1, 1, npair)).tolist())
    offs = np.array(offs)
    rest = u[~inside]
    if placed < ntot:                       # odd leftover: keep one on-line zero at the block centre
        rest = np.concatenate([rest, [0.0]])
    Won = float(np.sum(np.abs(g.F(rest))**2)); Woff = 2 * float(np.sum(g.h(offs, d)))
    ev = np.sort(np.concatenate([np.repeat(offs, 2), [0.0] * (ntot - placed)]))
    gu = np.linspace(lo, hi, 20001)
    dev = (np.searchsorted(ev, gu, side='right') - (np.searchsorted(u, gu) - np.searchsorted(u, lo)))
    return (Won + Woff) / base, ntot, float(np.max(np.abs(dev))), hi - lo
print('  adversarial blocks (delta = 0.3; K = off-line zeros; dev = max |N_mod - N_true| inside the block):')
for R in (10, 20, 40, 60, 80, 100, 130, 170, 220, 300):
    out = block(R)
    if out is None: continue
    Wb, K, dev, ln = out
    print(f'    R = {R:4d}: block length {ln:7.2f}, K = {K:4d}: W/(D 2pi||f||^2) = {Wb:+.4f}, dev {dev:.1f}', flush=True)
print(f'# maxrss {memcheck("end")/1e6:.0f} MB, {time.time()-t0:.0f}s')
