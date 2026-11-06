# hk_compare.py: compare our unconditional per-h counts F(h) (odd h <= HMAX) with Holmin-Jones-Kurlberg-McLeman-
# Petersen's F(h) for odd h < 10^6 (conditional on GRH; https://people.kth.se/~kurlberg/class_group_data/class_group_orders.txt).
# Usage: python3 hk_compare.py perh.txt cgo.txt HMAX
import sys
ours = {}
for l in open(sys.argv[1]):
    if l.startswith('#'): continue
    h, c, m = l.split(); ours[int(h)] = int(c)
hk = {}
for l in open(sys.argv[2]):
    if l.startswith('#') or not l.strip(): continue
    a = l.split(); hk[int(a[0])] = int(a[1])
H = int(sys.argv[3])
n = bad = so = sh = 0
for h in range(1, H + 1, 2):
    n += 1; so += ours.get(h, 0); sh += hk.get(h, -10**9)
    if ours.get(h) != hk.get(h): bad += 1; print('DIFF h=%d ours %s HK %s' % (h, ours.get(h), hk.get(h)))
print('odd h <= %d compared: %d, differences: %d, total fields ours %d, HK %d' % (H, n, bad, so, sh))
