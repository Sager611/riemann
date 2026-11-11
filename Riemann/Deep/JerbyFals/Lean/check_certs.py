# Check that every certificate in Main.lean is exactly the Arb ball [mid - rad, mid + rad] printed in the logs.
# Run from the repo root: python3 Riemann/Deep/JerbyFals/Lean/check_certs.py
import re
from fractions import Fraction as F
D = 'Riemann/Deep/JerbyFals/'
lean = open(D + 'Lean/Main.lean').read()
defs = {m[0]: (m[1], m[2], [F(x) for x in m[3].replace('(', ' ').replace(')', ' ').split()])
        for m in re.findall(r'def cert(\w\w) : Prop := (GramCertD?) (\d+)\s+([\s\S]*?)\n\n', lean)}
ok = True
for c, f in [('A', 'certify_cx.log'), ('B', 'certify_cx2.log'), ('C', 'certify_cx3.log')]:
    lines = open(D + 'logs/' + f).read().split('\n')
    recs, dz = {}, None
    for i, ln in enumerate(lines, 1):
        m = re.match(r'k=(\d+) g=\[(\S+) \+/- (\S+)\] Z=\[(\S+) \+/- (\S+)\]', ln)
        if m: recs[m[1]] = (i, F(m[2]) - F(m[3]), F(m[2]) + F(m[3]), F(m[4]) - F(m[5]), F(m[4]) + F(m[5]))
        m = re.search(r"Z'\(g_(\d+)\) in \[(\S+) \+/- (\S+)\]", ln)
        if m: dz = (i, int(m[1]), F(m[2]) - F(m[3]), F(m[2]) + F(m[3]))
    n = dz[1]
    for pos, k in [('L', n - 1), ('M', n), ('R', n + 1)]:
        kind, kk, nums = defs[c + pos]
        i, a, b, zl, zu = recs[str(k)]
        want = [a, b, zl, zu] + ([dz[2], dz[3]] if pos == 'M' else [])
        good = int(kk) == k and nums == want and (kind == 'GramCertD') == (pos == 'M')
        ok &= good
        print(f'cert{c}{pos}: k={k} {f} line {i}' + (f", Z' line {dz[0]}" if pos == 'M' else '') + (' MATCH' if good else ' MISMATCH'))
print('ALL MATCH' if ok else 'MISMATCH FOUND')
