# decomp.py -- jerbyconf: decomposition of the normalised viscosity s*Z'/Z(g_n) (s = 2 pi / log(g_n/2pi)) into the confined-pair term
# P = s sum_pair 1/(g_n - gamma) (Jerby 2610.21952 l. 2022-2041: |1/a_1 + 1/a_2| >= 2), the other zeros in [g_{n-6}, g_{n+6}], the drift s*pi/4
# (Z'/Z = pi/4 + sum 1/(t - gamma) + O(1/t)), and the remainder (zeros beyond 6 Gram intervals).  Input: gram_win -zeros 6 64 output.
import re, math, sys
L = open(sys.argv[1]).read().split('\n')
print('| n | t | s*mu | side | a_1, a_2 | pair term P | other zeros within 6 Gram intervals | drift s*pi/4 | far remainder | nearest outside zeros (Gram units beyond the walls) |')
print('|---|---|---|---|---|---|---|---|---|---|')
for i, line in enumerate(L):
    if not line.startswith('ZG'): continue
    n = int(re.search(r'n=(\d+)', line).group(1)); tc = float(re.search(r'tc=([\d.e+]+)', line).group(1)); W = int(re.search(r'W=(\d+)', line).group(1))
    Zp = float(re.search(r'Zp=([-\d.e+]+)', line).group(1))
    gz = [tuple(map(float, x.split(':'))) for x in line.split()[5:]]
    dg = [a for a, b in gz]; zg = [b for a, b in gz]
    zz = [float(x) for x in L[i+1].split()[2:]]
    dn = dg[W]; t = tc + dn; s = 2*math.pi/math.log(t/(2*math.pi)); Z = zg[W]
    smu = s*abs(Zp/Z); sgn = 1 if Zp/Z > 0 else -1
    left = [z for z in zz if dg[W-1] < z < dn]; right = [z for z in zz if dn < z < dg[W+1]]
    pair = left if len(left) == 2 else right if len(right) == 2 else []
    side = 'left' if pair is left and pair else 'right' if pair else 'zero-free' if not left and not right else f'{len(left)}+{len(right)}'
    P = s*sum(1/(dn - z) for z in pair); others = [z for z in zz if z not in pair]
    near = s*sum(1/(dn - z) for z in others); drift = s*math.pi/4
    far = s*Zp/Z - P - near - drift
    a = [(dn - z)/s for z in pair]
    outL = max([z for z in zz if z < dg[W-1]], default=None); outR = min([z for z in zz if z > dg[W+1]], default=None)
    ol = (dg[W-1] - outL)/s if outL is not None else float('nan'); orr = (outR - dg[W+1])/s if outR is not None else float('nan')
    # signs relative to the sign of Z'/Z, so that the columns add up to s*mu
    f = lambda x: sgn*x
    print(f"| {n} | {t:.6e} | {smu:.4f} | {side} | {', '.join(f'{x:.3f}' for x in a) if a else '–'} | {f(P):.4f} | {f(near):.4f} | {f(drift):.4f} | {f(far):.4f} | {ol:.3f}, {orr:.3f} |")
