#!/usr/bin/env python3
"""Check of littlewood.gp output: the reducible monic Littlewood polynomials with a root +-1 must be exactly
#{f(1) = 0 or f(-1) = 0} = 2 C(n,(n-1)/2) - [(n+1)/2 even] C((n+1)/2,(n+1)/4) C((n-1)/2,(n-3)/4) for odd n, 0 for even n.
Also prints the asymptotic 2 sqrt(2/(pi n)) (odd n) next to the observed reducible fraction."""
import math, sys
rows = [list(map(float, l.split())) for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
print("# n  2^n  reducible  red_with_root_pm1  exact_pm1_count  other_red  frac_red  2sqrt(2/(pi n))  n*other/2^n")
ok = True
for r in rows:
    n = int(r[0]); N, red, pm, other = int(r[1]), int(r[2]), int(r[3]), int(r[4])
    if n % 2:
        ex = 2 * math.comb(n, (n - 1) // 2)
        if ((n + 1) // 2) % 2 == 0:
            ex -= math.comb((n + 1) // 2, (n + 1) // 4) * math.comb((n - 1) // 2, (n - 3) // 4)
        asym = 2 * math.sqrt(2 / (math.pi * n))
    else:
        ex, asym = 0, 0.0
    ok &= (ex == pm)
    print("%3d %8d %7d %7d %7d %7d %.4f %.4f %.3f" % (n, N, red, pm, ex, other, red / N, asym, n * other / N))
print("# root-(+-1) counts match the exact formula:", ok)
