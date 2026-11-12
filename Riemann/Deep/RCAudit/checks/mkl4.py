import math
def fund(D):
    n = -D; best = 1
    for g in range(1, math.isqrt(n) + 1):
        if n % (g*g) == 0 and (D // (g*g)) % 4 in (0, 1): best = g
    return D // (best*best), best
fd = [D for D in range(-3, -10001, -1) if D % 4 in (0, 1) and fund(D)[1] == 1]
S = [D for D in fd if -D <= 200]
targets = [200 * 50 ** (k / 70) for k in range(1, 71)]
for t in targets:
    d = min(fd, key=lambda D: abs(-D - t))
    if d not in S: S.append(d)
F = [(D, 1) for D in S]
nf = []
for dK in [-3, -4, -7, -8, -11, -15, -20, -23, -24, -31, -39, -47, -71, -95, -119]:
    for f in [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 42, 60]:
        if -dK * f * f <= 10000: nf.append((dK, f))
nf = nf[::2] + [x for x in nf[1::2] if x[1] >= 20]
lines = ['read("l4g.gp");'] + [f'lem4g({dK}, 1, [1,5,20], [0,10,30], "l4g_out.txt");' for (dK, _) in F] + [f'lem4g({dK}, {f}, [1,5,20], [0,10,30], "l4g_out.txt");' for (dK, f) in nf]
open('l4run.gp', 'w').write('\n'.join(lines) + '\n')
print(len(F), 'fundamental,', len(nf), 'non-fundamental; max |D| nonfund', max(-d*f*f for d, f in nf))
