#!/usr/bin/env python3
"""Weight exponent of an off-line zero in DPC's moments (toy, float64).
At depth y (zero at beta = 1/2 + d, slit point x0 - iy, 0 <= y <= d), the hairpin weight is
  W(x0, y) = R(x0 - iy) R*(x0 - iy) = M(1/2 + y + i x0) * conj(M(1/2 - y + i x0)),  M(s) = sum c_q q^{1/2 - s}.
Mean squares: E|R(x-iy)|^2 = S_-(y) = sum |c_q|^2 q^{-2y},  E|R*(x-iy)|^2 = S_+(y) = sum |c_q|^2 q^{2y}.
Claim: E|W| ~ sqrt(S_+ S_-) = S(0) sinh(a)/a with a = y log Q  (i.e. T^{theta(beta-1/2)} up to 1/a), not S_+ ~ e^{2a}.
Also the 4th moment over x of R*(x - iY): E|R*|^4 / S_+^2 stays bounded (MVT regime, 2*len < range).
"""
import numpy as np, json, sys, time
rng = np.random.default_rng(1)
def mob_sieve(Q):
    mu = np.ones(Q+1, dtype=np.int64); mu[0] = 0
    isp = np.ones(Q+1, dtype=bool); isp[:2] = False
    for p in range(2, Q+1):
        if isp[p]:
            isp[2*p::p] = False
            mu[p::p] *= -1
            mu[p*p::p*p] = 0
    return mu
out = {}
t0 = time.time()
for Q in (100, 1000, 5000):
    mu = mob_sieve(Q)
    q = np.nonzero(mu)[0]; q = q[q >= 1]
    c = mu[q] / np.sqrt(q)
    lq = np.log(q)
    X = 1.0e6 + rng.uniform(0, 1.0e6, size=4000)   # range >> Q^2 so MVT holds for R^2
    rows = []
    for a in (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0):
        y = a / np.log(Q)
        Sm = np.sum(c**2 * q**(-2*y)); Sp = np.sum(c**2 * q**(2*y)); S0 = np.sum(c**2)
        absW = []; W2 = []; Rs4 = []
        for k in range(0, len(X), 500):
            x = X[k:k+500]
            E = np.exp(-1j * np.outer(x, lq))          # q^{-ix}
            Rv = E @ (c * q**(-y))                      # R(x - iy)
            Rsv = np.conj(E) @ (c * q**(y))             # R*(x - iy)
            W = Rv * Rsv
            absW.append(np.abs(W)); W2.append(np.abs(W)**2); Rs4.append(np.abs(Rsv)**4)
        absW = np.concatenate(absW); W2 = np.concatenate(W2); Rs4 = np.concatenate(Rs4)
        rows.append(dict(a=a, y=y, S0=S0, Sm=Sm, Sp=Sp,
                         EabsW_over_S0=float(absW.mean()/S0),
                         sqrtSpSm_over_S0=float(np.sqrt(Sp*Sm)/S0),
                         sinh_over_a=float(np.sinh(a)/a) if a > 0 else 1.0,
                         Sp_over_S0=float(Sp/S0),
                         ratio_EabsW_to_sqrtSpSm=float(absW.mean()/np.sqrt(Sp*Sm)),
                         ratio_EabsW_to_Sp=float(absW.mean()/Sp),
                         EW2_over_SpSm=float(W2.mean()/(Sp*Sm)),
                         ERs4_over_Sp2=float(Rs4.mean()/Sp**2)))
        r = rows[-1]
        print("Q=%5d a=%.2f  E|W|/S0=%9.4f  sqrt(S+S-)/S0=%9.4f  sinh(a)/a=%8.4f  S+/S0=%9.3f  E|W|/sqrt(S+S-)=%.3f  E|W|/S+=%.4f  E|W|^2/(S+S-)=%.3f  E|R*|^4/S+^2=%.3f" % (
            Q, a, r["EabsW_over_S0"], r["sqrtSpSm_over_S0"], r["sinh_over_a"], r["Sp_over_S0"], r["ratio_EabsW_to_sqrtSpSm"], r["ratio_EabsW_to_Sp"], r["EW2_over_SpSm"], r["ERs4_over_Sp2"]), flush=True)
    out[str(Q)] = rows
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "toy_weight.json", "w"), indent=1)
print("time %.1f s" % (time.time() - t0))
