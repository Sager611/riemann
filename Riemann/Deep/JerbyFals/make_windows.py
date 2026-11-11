# Build the log-uniform window list: J windows, t_j = 1e8 * 5000^(j/(J-1)), size clamp(1e10/N(t), 1e5, 1e6) Gram points.
import sys, mpmath
mpmath.mp.dps = 30
J = int(sys.argv[1]); out = sys.argv[2]
tot = 0
with open(out, 'w') as f:
    for j in range(J):
        t = mpmath.mpf(10)**8 * mpmath.mpf(5000)**(mpmath.mpf(j)/(J-1))
        N = float(mpmath.sqrt(t/(2*mpmath.pi)))
        cnt = int(max(1e5, min(1e6, 1e10/N)))
        n0 = int(mpmath.floor(mpmath.siegeltheta(t)/mpmath.pi)) - cnt//2
        f.write(f'L{j:02d} {n0} {cnt}\n'); tot += cnt*N/3.3e8
print('est secs', round(tot))
