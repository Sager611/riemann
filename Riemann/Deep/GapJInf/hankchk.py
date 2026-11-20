import sys, numpy as np
sys.path.insert(0, '.')
from ampmom import load
sv, am, kap, v = load(sys.argv[1], sys.argv[2])
m = sv[0::2]
for J in range(4, len(m)):
    n0 = J // 2; n1 = (J - 1) // 2
    H0 = np.array([[m[i + k] for k in range(n0 + 1)] for i in range(n0 + 1)])
    H1 = np.array([[m[i + k + 1] for k in range(n1 + 1)] for i in range(n1 + 1)])
    e0 = np.linalg.eigvalsh(H0); e1 = np.linalg.eigvalsh(H1)
    print('J=%2d: Stieltjes H0 min/max %.2e  H1 min/max %.2e  ratio m_J/m_{J-1} %.4f' % (J, e0[0]/e0[-1], e1[0]/e1[-1], m[J]/m[J-1]))
