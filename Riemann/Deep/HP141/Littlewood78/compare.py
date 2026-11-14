# compare.py -- key lines of the original Littlewood78/cert_output.txt (B = 0.87338196661320752) vs this folder's (theta141 rounded up)
import re, os
HERE = os.path.dirname(os.path.abspath(__file__))
old = open(os.path.join(HERE, '..', '..', 'Littlewood78', 'cert_output.txt')).read().split('\n')
new = open(os.path.join(HERE, 'cert_output.txt')).read().split('\n')
pats = [r'^limit constant', r'LamB = ', r'LamB2 = ', r'h\(-d\) >', r'=> \|L\(1,chi\)\| >=', r'small q', r'\(A\)|\(B\)|\(B2\)', r'h\(-d\) >= ',
        r'limit constant e\^gamma', r'\|L\(1,chi\)\| <= 1[56] loglog']
def pick(lines):
    return [l.strip() for l in lines if any(re.search(p, l) for p in pats)]
po, pn = pick(old), pick(new)
same = 0
for a, b in zip(po, pn):
    if a == b: same += 1; continue
    print('OLD:', a[:230]); print('NEW:', b[:230]); print()
print(f'{same} key lines identical; {len(po)} old / {len(pn)} new key lines')
