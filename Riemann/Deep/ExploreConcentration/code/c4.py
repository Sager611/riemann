# explore44 c4: Littlewood's lemma for DH on [sigma, 3] x [10, 420] with the known off-line zeros (explore23 list): zero sum = mean shift of log|F|.
import numpy as np, json, time, sys
sys.path.insert(0, '/Users/adriansager/playground/riemann/.scratch/explore44')
from dhlib import hurwitz, f as fdh
t0 = time.time()
Z = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/dh_zeros.json'))['offline']
T1, T2, SI = 10.0, 420.0, 3.0
for sig in [0.6, 0.7]:
    t = np.arange(T1, T2, 0.002); tm = 0.5*(t[1:] + t[:-1])
    lv = np.log(np.abs(fdh(sig + 1j*tm, M=300))).sum()*0.002
    li = np.log(np.abs(fdh(SI + 1j*tm, M=300))).sum()*0.002
    def harg(T):
        u = np.arange(SI, sig - 1e-12, -0.0005)
        a = np.unwrap(np.angle(fdh(u + 1j*T, M=300)))
        a = a - 2*np.pi*np.round(a[0]/(2*np.pi))
        um = 0.5*(u[1:] + u[:-1]); am = 0.5*(a[1:] + a[:-1])
        return float(np.sum(am*np.abs(np.diff(u))))   # integral over [sig, SI] of arg F(u + iT) du
    a1, a2 = harg(T1), harg(T2)
    zs = sum(b - sig for b, g in Z if b > sig and T1 < g < T2)
    lhs1 = lv - li + a1 - a2; lhs2 = lv - li - a1 + a2
    print('sigma', sig, 'int log|F(sig+it)| %.5f  int log|F(3+it)| %.5f  harg(T1) %.5f harg(T2) %.5f' % (lv, li, a1, a2))
    print('   Littlewood LHS (+a1 - a2) %.5f   (-a1 + a2) %.5f   2 pi sum(beta - sigma) %.5f  (%d zeros)' % (lhs1, lhs2, 2*np.pi*zs, sum(1 for b, g in Z if b > sig)), flush=True)
    print('   mean shift E log|F(sig+it)| over [10,420] = %.5f ; zero part 2pi sum/(T2-T1) = %.5f' % (lv/(T2 - T1), 2*np.pi*zs/(T2 - T1)))
print('done', round(time.time() - t0, 1))
