"""k = 0 check of zk.forms against ExploreShape crosscheck (r = 1, theta = 1/4, kappa = 3.18, degree 5, v = 1.26, unshifted
mollifier): h1 = 1.04053590, <(v-S)^k>/<1> = 1, 0.030934125, 0.095036118, 0.0088906238, 0.020059772 (k = 0..4)."""
import sys, numpy as np, scipy.linalg as sl, zk
nodes = tuple(int(x) for x in sys.argv[1].split(','))
M, herm, imag = zk.forms(0, 3.18, 0.0, 1.26, 5, 4, th=0.25, nodes=nodes)
# note zk uses (S - v0)^n; BM's frequency is (v - S): sign (-1)^n
lam = 3.18
ev, U = sl.eigh(M[0], lam**2 * M[2])
b = U[:, -1]; n0 = b @ M[0] @ b
print('nodes', nodes, 'h1 = %.8f' % ev[-1], 'ratios', ['%.10g' % ((-1)**n * (b @ M[n] @ b) / n0) for n in range(5)], 'herm %.1e imag %.1e' % (herm, imag))
