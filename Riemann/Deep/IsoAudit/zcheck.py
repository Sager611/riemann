# zcheck.py: fraction of c-isolated zeros among Odlyzko's first 100,000 zeros (zeros1), definition of [La] Thm 1.1 with T = gamma_100000.
# Run: uv run --no-project --with numpy python3 zcheck.py ../../../.scratch/isoaudit/zeros1
import numpy as np, sys
g = np.loadtxt(sys.argv[1]); n = len(g); T = g[-1]; L = np.log(T)
d = np.diff(g)
print('zeros: %d, T = %.6f, min gap %.6f (normalised by 2pi/log T: %.4f)' % (n, T, d.min(), d.min()*L/(2*np.pi)))
for c in (0.001, 0.00476, 0.01, 0.02485, 0.05495, 0.06427, 0.1, 1/6, 0.5, 1.0):
    thr = 2*np.pi*c/L
    iso = (d[1:] > thr) & (d[:-1] > thr)       # zeros n = 2..n-1 (1-based): both neighbours more than thr away
    # local version 2 pi c / log gamma_n (abstract's definition)
    thr_loc = 2*np.pi*c/np.log(g[1:-1])
    iso_loc = (d[1:] > thr_loc) & (d[:-1] > thr_loc)
    print('c = %-8.5f B_0/N = %.6f   (log gamma_n version %.6f)' % (c, iso.sum()/n, iso_loc.sum()/n))
