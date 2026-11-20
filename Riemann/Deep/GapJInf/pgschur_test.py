import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pgschur as PS
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ExploreShape'))
import exact4 as E4
a2 = [-5899179108889603/6134092019850000, 1, 14789027537/948575124080, -4037122085/176350902122, 1982999731/771179571628]
p12, p32 = 3691572375143311/104390384155500000, 4037122085/176350902122
mine = ([-0.95767, 1, 0.011945, -0.023295, 0.0025546], 0.03403, 0.02335)
for name, (a, q12, q32) in [('cert4', (a2, p12, p32)), ('mineLP', mine)]:
    x = np.array(list(a) + [q12, q32])
    for g in [np.pi, 0.9*np.pi]:
        out = []
        for par in (1, 0):
            D = PS.mode_data(g, 4, par, 20000)
            Amin, R, Bt = PS.reduce(D, x)
            out.append(R[0, 0])
        # exact4 S_odd, S_even + (a3+phi32)/2 times (2/g)*(2/g)*(g/2)... compare signs/ratios
        k = np.arange(1, 20001); xk = (k*np.pi/g)**2
        rho, n3 = E4.rho_stable(xk, a, q12, q32)
        so = rho[0::2].sum(); se = rho[1::2].sum() + 0.5*(a[3]+q32)
        print(name, 'g/pi=%.2f' % (g/np.pi), 'R_odd=%.4e R_even=%.4e | exact4 S_odd=%.4e S_even+c=%.4e | ratio %.6f %.6f' % (out[0], out[1], so, se, out[0]/so, out[1]/se))
