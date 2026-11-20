# Contrapositive test of the windowed Landau theorem (zeroFree_of_rpow_window, ZeroFreeOfRpowWindow.lean:11-14):
# toy Euler products F(s) = zeta(s) * P(2^-s) with known zeros off the line; G_F = proper Riesz function
# (ramp - y log y + (1 + gamma + kappa) y, kappa = -(-P'/P)(1), the analogue of + y L'/L(1,chi) in Gchi).
# The window bound  G_F(y) <= C y^theta on [Y, Y^(1+delta)]  must fail for all large Y when Re rho0 > theta.
import numpy as np, math, sys
YMAX = 4_000_000
EG = 0.5772156649015329
L2 = math.log(2.0)

def lam_array(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if s[p]: s[p*p::p] = False
    pr = np.nonzero(s)[0]
    lam = np.zeros(N + 1); lam[pr] = np.log(pr)
    for p in pr[pr <= int(N**0.5)]:
        q = p * p
        while q <= N: lam[q] = math.log(p); q *= p
    return lam

lam = lam_array(YMAX)
n = np.arange(YMAX + 1, dtype=np.float64)
S1 = np.cumsum(np.where(n > 0, lam / np.maximum(n, 1), 0.0)); S0 = np.cumsum(lam)
del lam
y = n.copy(); y[0] = 1.0
Gz = y * S1 - S0 - y * np.log(y) + (1 + EG) * y
del S1, S0

def T_part(coef):  # coef(k) = extra Lambda_F(2^k) / log 2 (beyond Lambda), returns ramp part + kappa*y
    T = np.zeros_like(y); kap = 0.0
    for k in range(1, 200):
        ck = coef(k)
        kap += -L2 * ck * 2.0**(-k)      # -(-P'/P)(1) contribution: Lambda_extra(2^k)/2^k summed with sign
        if 2**k <= YMAX:
            m = y >= 2**k
            T[m] += L2 * ck * (y[m] / 2**k - 1)
    return T + kap * y, kap

def window_test(G, theta, C, delta=0.5, label=""):
    r = G / y**theta
    res = []
    j = 16
    while True:
        Y = 2 ** (j / 4); Yt = Y ** (1 + delta)
        if Yt > YMAX: break
        lo, hi = int(math.ceil(Y)), int(math.floor(Yt))
        i = lo + int(np.argmax(r[lo:hi + 1]))
        res.append((Y, r[i], i))
        j += 1
    viol = np.nonzero(r[16:] > C)[0]
    yfirst = int(viol[0]) + 16 if len(viol) else None
    fails = [Y for (Y, mx, i) in res if mx > C]
    first = fails[0] if fails else None
    allafter = first is not None and all(mx > C for (Y, mx, i) in res if Y >= first)
    print(f"{label} theta={theta} C={C}: windows {len(res)} (Y=2^4..{res[-1][0]:.0f}); first failing Y = {first if first is None else round(first)}; all later fail: {allafter}; first violating y = {yfirst}; max G/y^theta over all = {max(m for _,m,_ in res):.3f}")
    return first, res

out = []
# zeta itself: no zero off the line -> the bound should hold on every window
for th, C in [(0.6, 0.5), (0.55, 0.5)]:
    window_test(Gz, th, C, label="zeta alone:")
print("zeta alone: max |Gz(y)|/sqrt(y) for 1e3<=y<=4e6 =", float(np.max(np.abs(Gz[1000:]) / np.sqrt(y[1000:]))))

for c in (0.8, 0.9):
    # Toy A: zeta(s)(1 + 4^(c-s)): zeros c + i pi(2k+1)/(2 log 2), no real zero.
    # -P'/P(s) = -log2 sum_j 2*(-1)^j 4^{jc} 4^{-js}  -> extra Lambda at 2^(2j): -2(-1)^j 2^{2jc} log 2
    coefA = lambda k, c=c: (-2.0 * (-1) ** (k // 2) * 2.0 ** (k * c)) if k % 2 == 0 else 0.0
    TA, kapA = T_part(coefA)
    GA = Gz + TA
    rho0 = complex(c, math.pi / (2 * L2))
    A0 = 2 / abs(rho0 * (rho0 - 1))
    env = float(np.max(TA[10**5:] / y[10**5:] ** c))
    print(f"\nToy A c={c}: kappa={kapA:.6f}; leading zero pair rho0={rho0:.4f}, explicit-formula amplitude 2/|rho0(rho0-1)| = {A0:.4f}; empirical max T/y^c (y>=1e5) = {env:.4f}")
    for th in (0.6,):
        for C in (0.5, 1, 2, 4):
            pred_y = (C / A0) ** (1 / (c - th))   # first y with A0 y^(c-theta) > C
            first, res = window_test(GA, th, C, label=f"Toy A c={c}:")
            pred_env = (C / env) ** (1 / (c - th))
            print(f"   predicted first violating y: (C/A0)^(1/(c-theta)) = {pred_y:.3g} (leading pair only); (C/envelope)^(1/(c-theta)) = {pred_env:.3g} (all zeros of P); first failing window then Y ~ y^(2/3) = {pred_env**(2/3):.3g}")
    # Toy B: zeta(s)(1 - 2^(c-s)): real zero at c (+ c + 2 pi i k/log 2)
    coefB = lambda k, c=c: -(2.0 ** (k * c))
    TB, kapB = T_part(coefB)
    GB = Gz + TB
    print(f"Toy B c={c} (real zero at c): min T/y^c (y>=1e5) = {float(np.min(TB[10**5:]/y[10**5:]**c)):.4f} > 0 (real-zero term y^c/(c(1-c)) = {1/(c*(1-c)):.3f} y^c)")
    for C in (1, 4):
        window_test(GB, 0.6, C, label=f"Toy B c={c}:")
