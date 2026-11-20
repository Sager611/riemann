# validate fz against mpmath at a few points (zeta and DH, F and F')
import numpy as np, mpmath as mp, math, sys
sys.path.insert(0, '.'); from fz import Fn, KAP
mp.mp.dps = 20
pts = [0.5+14.134725j, 0.5+85.0j, 0.8085+85.699j, 0.3+200.3j, 2.0+1000.5j, 0.5+10000.7j]
for kind in ('zeta', 'dh'):
    f = Fn(kind, 10200); F, Fp = f(np.array(pts)); err = 0
    for k, s in enumerate(pts):
        if kind == 'zeta': a = mp.zeta(s); ap = mp.zeta(s, derivative=1)
        else: a = mp.dirichlet(s, [0, 1, KAP, -KAP, -1]); ap = mp.dirichlet(s, [0, 1, KAP, -KAP, -1], derivative=1)
        e = max(abs(complex(a) - F[k]), abs(complex(ap) - Fp[k]))/max(1, abs(complex(a)), abs(complex(ap))); err = max(err, e)
    # functional equation check on the line: |chi|=1, Z real
    t = np.array([100.3, 1234.5]); s = 0.5 + 1j*t; F2, _ = f(s); Z = np.exp(1j*f.theta(t))*F2
    print(kind, 'max rel err vs mpmath', f'{err:.2e}', ' Im Z/|Z|', np.abs(np.imag(Z))/np.abs(Z))
