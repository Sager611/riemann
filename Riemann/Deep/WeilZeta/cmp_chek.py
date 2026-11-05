# cmp_chek.py: compare our rigorous scan (scan cache JSON) with Chekmachev's retained nodes (certificate data, Zenodo 23847003).
# Usage: cmp_chek.py <our scan json> <their nodes_odd.csv> <their nodes_even.csv> K2
import sys, json, csv
from decimal import Decimal
sc = json.load(open(sys.argv[1])); K2 = int(sys.argv[4]); J = K2
def load(fn):
    out = {}
    with open(fn) as f:
        r = csv.reader(f); next(r)
        for row in r: out[int(row[0])] = (Decimal(row[1]), Decimal(row[2]))
    return out
no = load(sys.argv[2]); ne = load(sys.argv[3])
selo = set(sc['selo']); sele = set(sc['sele'])
th_o = {k for k, (m, r) in no.items() if k > K2 and m - r < 1}
th_e = {j for j, (m, r) in ne.items() if j > J and m - r < 0}
print(f"ours: kscan {sc['ks']} kscan_e {sc['kse']}, selected odd {len(selo)} even {len(sele)}; scan info {sc.get('scan')}")
print(f"theirs (retained, beyond K2 = J = {K2}): odd with mu < 1: {len(th_o)}; even with mu < 0: {len(th_e)}; retained odd {len(no)} even {len(ne)}")
print(f"odd: ours-not-theirs {len(selo - th_o)} {sorted(selo - th_o)[:10]}; theirs-not-ours {len(th_o - selo)} {sorted(th_o - selo)[:10]}; ours not retained by them {len(selo - set(no))}")
print(f"even: ours-not-theirs {len(sele - th_e)} {sorted(sele - th_e)[:10]}; theirs-not-ours {len(th_e - sele)} {sorted(th_e - sele)[:10]}; ours not retained by them {len(sele - set(ne))}")
if selo: print(f"largest selected odd k: ours {max(selo)}, theirs {max(th_o) if th_o else None}; even j: ours {max(sele) if sele else None}, theirs {max(th_e) if th_e else None}")
