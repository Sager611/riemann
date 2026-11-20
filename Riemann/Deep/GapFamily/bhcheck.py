"""Reproduce Bui-Hall's printed h_k (Wirtinger) from zk.forms: h = c0/(lam^2 c1), P = 1 - a x, V = v+1."""
import sys, json, time
import numpy as np
import zk
k = int(sys.argv[1]); lam = float(sys.argv[2]); a = float(sys.argv[3]); v = float(sys.argv[4]); eta = float(sys.argv[5])
nodes = tuple(int(x) for x in sys.argv[6].split(',')) if len(sys.argv) > 6 else (11, 12, 10, 12, 10, 16, 16)
t = time.time()
M, herm, imag = zk.forms(k, lam, eta, v + 1, 1, 2, nodes=nodes)
b = np.array([1.0, -a])
c0 = b @ M[0] @ b; c1 = b @ M[2] @ b
print('k=%d lam=%.4f P=1-%.3gx v=%.3f eta=%.4f nodes=%s: h=%.8f  c0=%.10g c1=%.10g herm=%.1e imag=%.1e (%.1fs)' % (k, lam, a, v, eta, nodes, c0 / (lam**2 * c1), c0, c1, herm, imag, time.time() - t))
