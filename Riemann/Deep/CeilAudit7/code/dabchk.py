import numpy as np, sys
sys.argv = ['x', 'none']
exec(open('dab.py').read().split("mode = sys.argv[1]")[0])
for (a, b) in [(0.5, -0.25), (0.8, -0.4), (0.25, 0.5), (1.3, 0.2), (0.7, -1.1), (2.0, -0.5)]:
    v = [Dab(N, a, b)[0] for N in (400, 800)]
    ext = 2*v[1]-v[0]
    print(f"D({a},{b}): N=400 {v[0]:+.6f} N=800 {v[1]:+.6f} Richardson {ext:+.6f} formula {np.sinc(a-b)-np.sinc(a)*np.sinc(b):+.6f}")
