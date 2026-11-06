# Uncertified estimate (floats): asymptotic LP bound for a double real zero with smoothed weights. See REPORT item 6.

import numpy as np, math
# asymptotic LP for a double real zero: need 2*Phi(a+b) - Phi(a) > f(0)/2 = 1/2 for some admissible f, a >= 0
def Phi_tri(z, T):  # f(t) = (1 - t/T)_+, Phi(z) = 1/z - (1 - e^{-zT})/(z^2 T)
    if abs(z) < 1e-12: return T/2
    return 1/z - (1 - math.exp(-z*T))/(z*z*T)
def best_b(Phi, params):
    lo, hi = 0.0, 3.0
    for _ in range(50):
        b = (lo+hi)/2
        ok = False
        for T in params:
            for a in np.linspace(0, 3, 301):
                if 2*Phi(a+b, T) - Phi(a, T) > 0.5: ok = True; break
            if ok: break
        if ok: lo = b
        else: hi = b
    return lo
Ts = list(np.linspace(0.2, 20, 199))
print("triangle kernel: b* =", best_b(Phi_tri, Ts))
# plain kernel check: Phi = 1/z
print("plain 1/z: b* =", best_b(lambda z,T: 1/z if z>0 else 1e9, [1.0]), " (6-4sqrt2 =", 6-4*math.sqrt(2), ")")
# Fejer-type f(t) = (1-t/T)_+ * (cos(pi t/T)... ) skip; also f = (1 - t/T)^2_+ (still PD? test only)
def Phi_num(z, T, f):
    t = np.linspace(0, T, 4001); y = f(t/T)*np.exp(-z*t); return np.trapezoid(y, t)
for name, f in [("(1-u)^2", lambda u: (1-u)**2), ("cos^2(pi u/2)", lambda u: np.cos(np.pi*u/2)**2), ("(1-u)cos(pi u)+sin(pi u)/pi", lambda u: (1-u)*np.cos(np.pi*u)+np.sin(np.pi*u)/np.pi)]:
    Ts2 = list(np.linspace(0.5, 12, 47))
    print(name, "b* ~", best_b(lambda z,T: Phi_num(z,T,f), Ts2))
