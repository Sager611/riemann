# hkcheck.py PERH CGO: our per-h counts (perh_h10000.txt) vs Holmin et al.'s GRH counts (class_group_orders.txt, odd h).
import sys
ours = {int(a[0]): int(a[1]) for a in (l.split() for l in open(sys.argv[1])) if a and a[0].isdigit()}
hk = {}
for l in open(sys.argv[2]):
    a = l.split()
    if len(a) >= 2 and a[0].isdigit(): hk[int(a[0])] = int(a[1])
odd = range(1, 10000, 2)
diff = [h for h in odd if ours.get(h) != hk.get(h)]
print('odd h compared', len(odd), 'differences', len(diff), diff[:5], 'sum ours', sum(ours[h] for h in odd), 'sum HK', sum(hk.get(h, 0) for h in odd))
