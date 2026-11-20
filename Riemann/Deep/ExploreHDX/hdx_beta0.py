# hdx_beta0.py (explore64): beta0 = min over integers 2 <= N' <= 100 of the best integer-cycle bound (hdx_path), all factorizations
# into <= 4 prime powers; for every window c <= 101 and link c' with an integer N' in (c', c) (N' > c' strictly), Lam*(c, c') >= beta0.
import sys, math
sys.path.insert(0, 'Riemann/Deep/ExploreHDX')
from hdx_path import factorizations, best, certify, ppow
worst = []
for Np in range(2, 101):
    vals = []
    for fac in factorizations(Np, kmax=4):
        ws = [math.log(ppow(q)[0])/math.sqrt(q) for q in fac]
        v, (perm, al) = best(ws)
        pool = list(fac); qs = []
        for w in perm:
            j = min(range(len(pool)), key=lambda t: abs(math.log(ppow(pool[t])[0])/math.sqrt(pool[t]) - w)); qs.append(pool.pop(j))
        vals.append((float(certify(perm, al, qs).lower()), v, fac))
    vals.sort(reverse=True); worst.append((vals[0][0], Np, vals[0][2]))
worst.sort()
print('# five weakest N\' <= 100 (arb lower bound of the best cycle, factorization):')
for v, Np, fac in worst[:5]: print(f"  N'={Np:3d} = {'*'.join(map(str, fac))}: {v:.4f}")
print(f'beta0 >= {worst[0][0]:.4f}')
