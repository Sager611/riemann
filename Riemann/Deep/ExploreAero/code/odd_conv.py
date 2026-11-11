import numpy as np, sys
sys.path.insert(0, '.')
from weilleg import WeilLeg
from scipy.special import spherical_in
class O(WeilLeg):
    def __init__(s, a, K):
        super().__init__(a, K); s.deg = 2*np.arange(K) + 1; s.nrm = np.sqrt((2*s.deg + 1)/(2*a))
        from numpy.polynomial import legendre as Lg; s.nin = s.deg[-1] + 40; s.nout = 2*s.deg[-1] + 60; s.ti, s.wi = Lg.leggauss(s.nin)
    def pole(s): return s.nrm*2*s.a*spherical_in(s.deg, s.a/2)
for K in [10, 20, 40, 80]:
    W = O(np.log(2)/2, K).build(); Q = -2*np.outer(W.l, W.l) + W.A - W.P
    print('odd c=2 K', K, np.linalg.eigvalsh(Q)[:2])
