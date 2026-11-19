# nbaudit c3: write the certified V1 table for every N <= 1280 from the exact intervals of the lean rerun, outward-rounded to 25 significant digits
# usage: python3 c3_table.py   (stdlib only)
from fractions import Fraction as Fr
from decimal import Decimal, getcontext, ROUND_FLOOR, ROUND_CEILING
import math
getcontext().prec = 60
def dec(fr, mode):
    q = Decimal(fr.numerator) / Decimal(fr.denominator)          # 60 digits, then outward to 25 significant
    e = q.adjusted() - 24
    return q.quantize(Decimal(1).scaleb(e), rounding=mode)
out = open('/Users/adriansager/playground/riemann/Riemann/Deep/NBAudit/logs/dN2_certified_1280.txt', 'w')
out.write('# N  lo  hi : rigorous enclosure lo <= d_N^2 <= hi (Baez-Duarte V1, L^2(0,inf)), outward-rounded from the exact intervals of code/c1_s2lean.py 1280 -21\n')
w41 = None
for line in open('/tmp/nbaudit/loh_V1_lean.txt'):
    N, a, b, c, d = map(int, line.split()); lo, hi = Fr(a, b), Fr(c, d)
    out.write(f'{N} {dec(lo, ROUND_FLOOR)} {dec(hi, ROUND_CEILING)}\n')
    if N == 41: w41 = float(hi - lo); m41 = float(lo) * math.log(41)
out.close()
print('width at N=41:', w41, ' (d^2 log 41 approx', m41, ')')
