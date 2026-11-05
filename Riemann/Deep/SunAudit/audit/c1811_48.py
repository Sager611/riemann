import sys
sys.path.insert(0, '.')
from lib import *
N = int(sys.argv[1])
for n in range(2, N + 1):
    full = (1 << n) - 1
    f = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        f[1 << v][v] = 1  # sum 0
    for mask in range(1, 1 << n):
        row = f[mask]
        for last in range(n):
            b = row[last]
            if not b: continue
            for w in range(n):
                if mask >> w & 1: continue
                f[mask | (1 << w)][w] |= b << ((last + 1) * (w + 1))
    allsums = 0
    for last in range(n): allsums |= f[full][last]
    endn = f[full][n - 1]
    # (i): some 2^m+1 in allsums
    m = 0; hit = []
    while (1 << m) + 1 < allsums.bit_length():
        if allsums >> ((1 << m) + 1) & 1: hit.append((1 << m) + 1)
        m += 1
    pw = [1 << m for m in range(0, endn.bit_length()) if (1 << m) < endn.bit_length() and endn >> (1 << m) & 1]
    print(n, '4.8(i) hits 2^m+1:', hit, '| 4.8(ii) powers of 2 with pi(n)=n:', pw, flush=True)
