# zhp.py (lostlaw): zeta zeros at high precision (mpmath zetazero), output "k value" lines. Usage: zhp.py dps k1,k2,...
import sys, mpmath as mp
dps = int(sys.argv[1]); mp.mp.dps = dps + 10
for k in [int(t) for t in sys.argv[2].split(',')]:
    print(k, mp.nstr(mp.zetazero(k).imag, dps), flush=True)
