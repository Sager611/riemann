# chek_spot.py: (1) Chekmachev's retained sets contain every odd k <= K2 and every j <= J (premise of Thm 5.10 (i));
# (2) our symbol mu (cert_zeta.SymZeta, Arb, 256 bits) equals his published 40-digit midpoints at random retained nodes.
# Usage: chek_spot.py c K2 nodes_odd.csv nodes_even.csv nsample
import sys, csv, random
from decimal import Decimal
from flint import arb, ctx, fmpq
from cert_zeta import SymZeta
c = int(sys.argv[1]); K2 = int(sys.argv[2]); ns = int(sys.argv[5])
def load(fn):
    out = {}
    with open(fn) as f:
        r = csv.reader(f); next(r)
        for row in r: out[int(row[0])] = (row[1], row[2])
    return out
no = load(sys.argv[3]); ne = load(sys.argv[4])
miss_o = [k for k in range(1, K2 + 1, 2) if k not in no]; miss_e = [j for j in range(0, K2 + 1) if j not in ne]
print(f"c={c}: odd k <= {K2} missing from retained: {len(miss_o)}; j <= {K2} missing: {len(miss_e)} {miss_e[:3]}")
ctx.prec = 256; S = SymZeta(fmpq(c))
random.seed(1)
so = random.sample(sorted(no), min(ns, len(no))); se = random.sample(sorted(ne), min(ns, len(ne)))
worst = Decimal(0); bad = 0; wk = None
for kind, sample, tab in (('odd', so, no), ('even', se, ne)):
    for i in sample:
        k = i if kind == 'odd' else 2*i
        m = S.mu(k); theirs = Decimal(tab[i][0])
        ours = Decimal(m.mid().str(45, radius=False))
        d = abs(ours - theirs)
        if d > worst: worst = d; wk = (kind, i)
        if not m.overlaps(arb(tab[i][0]) + arb(0, 1e-37)): bad += 1     # their midpoints: 40 digits, truncated
print(f"spot check: {len(so)} odd + {len(se)} even retained nodes (incl. k up to {max(so)} / j up to {max(se)}): max |our mid - their mid| = {float(worst):.2e} at {wk}; balls disjoint from their value +- 1e-37: {bad}")
