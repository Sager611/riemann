# t3_check.py (explore48): validate count.py: N(T) integer, N(T) - N0(T) = 0 below and 2 just above known off-line zeros.
import math, sys, time
sys.path.insert(0, '.')
from count import NE, N0_scan, Etilde_real_excess
cases = [('D-20 principal', 0.0, math.sqrt(5), [10.0, 15.0, 16.5, 20.0]),
         ('D-15 principal', -0.5, math.sqrt(15) / 2, [11.0, 12.5, 20.0]),
         ('D-163 h1', -0.5, math.sqrt(163) / 2, [30.0, 62.2, 84.6, 97.2])]
for lab, x, y, Ts in cases:
    print(lab, 'Etilde(3)-1 <= %.4f' % Etilde_real_excess(x, y, 3.0))
    for T in Ts:
        t0 = time.time(); v, a = NE(T, x, y); n0, pos = N0_scan(T, x, y)
        print('  T=%.2f  N(T)=%.4f  N0=%d  offline=%.4f  argEt=%.3f  [%.1fs]' % (T, v, n0, v - n0, a, time.time() - t0), flush=True)
