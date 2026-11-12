import numpy as np, math, sys, flint
from flint import acb, ctx
ctx.prec = 80
sys.argv = ['levin.py','zeta','10000','2','0.01']
exec(open('levin.py').read().split("J = int(round")[0])
for (sg_, t) in [(0.37, 10001.3), (0.5, 10001.0), (1.2, 10001.7)]:
    s = np.array([sg_ + 1j*t])
    ph = np.exp(-1j*t*logm)
    F = (coef*np.exp(-sg_*logm)) @ ph + tail(s)[0][0]
    Fp = (-coef*logm*np.exp(-sg_*logm)) @ ph + tail(s)[1][0]
    ser = flint.acb_series([acb(sg_, t), 1], prec=2).zeta()
    z0 = complex(ser.coeffs()[0]); z1 = complex(ser.coeffs()[1])
    print(sg_, t, abs(F - z0)/abs(z0), abs(Fp - z1)/abs(z1))
# DH check at known zero
sys.argv = ['levin.py','dh','60','30','0.01']
exec(open('levin.py').read().split("J = int(round")[0])
for (sg_, t) in [(0.808517, 85.699348), (0.5, 85.0)]:
    s = np.array([sg_ + 1j*t]); ph = np.exp(-1j*t*logm)
    F = (coef*np.exp(-sg_*logm)) @ ph + tail(s)[0][0]
    print('DH', sg_, t, abs(F))
