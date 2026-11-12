#!/usr/bin/env python3
"""uaudit N1: independent falsification test of the slit identity (H) of DPCUncond PROOF Step 3,
with coincident, clustered and same-ordinate cuts (cases the dpcuncond toy did not contain).

X(s) = exp(alpha*s) * prod_rho (s-rho)^m ; log X = alpha*s + sum m*Log(s-rho) (principal Log = horizontal
continuation from Re s=+inf: every cut is {Im s = gamma, Re s < beta}).
V(t) = (log X(1/2+i(t+h/2)) - log X(1/2+i(t-h/2)))/pi,  F = g R Rstar P(V),  P complex cubic.
Slits in -1 <= Im t <= 0: x = gamma -+ h/2, depth beta-1/2 (beta>1/2). Merged by x.
(H): int_a^b F dx = int_{a-i}^{b-i} F dt + int_{b-i}^{b} F dt - int_{a-i}^{a} F dt - i sum_x int_0^{dep(x)} (F_R-F_L)(x-iy) dy.
Also evaluated: W0 = no hairpins; W1 = one hairpin per zero and sign with its own increment J_rho
(the 'separate increments' rule PROOF Step 3 rejects).
"""
import sys, json, time, cmath, math
import numpy as np
from scipy.integrate import quad

H = 0.8                      # window h
ALPHA = 0.07 + 0.03j         # zero-free entire factor exp(alpha s)
PC = [0.3, 0.5j, 1 - 2j, 1]  # P(w) = 0.3 + 0.5i w + (1-2i) w^2 + w^3
DP = [0.5j, 2 * (1 - 2j), 3]  # P'
Q = 30

def mob(n):
    r, p, m = 1, 2, n
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            r = -r
        p += 1
    return -r if m > 1 else r

QS = np.array([q for q in range(1, Q + 1) if mob(q) != 0], dtype=float)
CS = np.array([mob(int(q)) / math.sqrt(q) * (1 - math.log(q) / math.log(Q + 1)) for q in QS])
LQ = np.log(QS)

def R(t):   # sum c_q q^{-it}
    return complex(np.sum(CS * np.exp(-1j * t * LQ)))

def Rs(t):  # sum conj(c_q) q^{it}
    return complex(np.sum(CS * np.exp(1j * t * LQ)))

TC, W = 0.0, 2.5
def g(t):
    return cmath.exp(-(t - TC) ** 2 / (2 * W * W))

def P(w):
    return PC[0] + PC[1] * w + PC[2] * w * w + PC[3] * w ** 3

def logX(s, zeros, side=0, tol=1e-9):
    """side=0 principal; side=-1: factors exactly on their cut take arg -pi (Im s -> gamma from below);
    side=+1: arg +pi."""
    v = ALPHA * s
    for (beta, gam, m) in zeros:
        z = s - complex(beta, gam)
        if z == 0:
            z = complex(1e-300, 0)
        if side != 0 and abs(z.imag) < tol and z.real < 0:
            v += m * complex(math.log(-z.real), side * math.pi)
        else:
            v += m * cmath.log(z)
    return v

def V(t, zeros, side=0):
    sp = 0.5 + 1j * (t + H / 2)
    sm = 0.5 + 1j * (t - H / 2)
    return (logX(sp, zeros, side) - logX(sm, zeros, side)) / math.pi

def F(t, zeros, side=0):
    return g(t) * R(t) * Rs(t) * P(V(t, zeros, side))

USE_MP = False

def cquad(f, a, b, pts=None):
    if USE_MP:   # tanh-sinh: endpoint log^3 singularities (corner case)
        import mpmath as mp
        mp.mp.dps = 18
        nodes = [a] + sorted(p for p in (pts or []) if a < p < b) + [b]
        v = mp.quad(lambda x: mp.mpc(f(float(x))), nodes, maxdegree=10, error=True)
        return complex(v[0]), float(v[1])
    kw = dict(limit=800, epsabs=1e-13, epsrel=1e-11)
    if pts:
        pts = sorted(p for p in pts if a < p < b)
        kw['points'] = pts if pts else None
        if not pts:
            kw.pop('points')
    re = quad(lambda x: f(x).real, a, b, **kw)
    im = quad(lambda x: f(x).imag, a, b, **kw)
    return complex(re[0], im[0]), abs(re[1]) + abs(im[1])

def slits(zeros):
    """merged slits: x -> list of (depth, sign J/2i, mult)"""
    S = {}
    for (beta, gam, m) in zeros:
        if beta > 0.5:
            for (x, sg) in ((gam - H / 2, +1), (gam + H / 2, -1)):
                key = round(x, 9)
                S.setdefault(key, []).append((beta - 0.5, sg, m))
    return S

def run(zeros, a=-12.0, b=12.0):
    sing = []
    for (beta, gam, m) in zeros:
        sing += [gam - H / 2, gam + H / 2]
    I_real, e1 = cquad(lambda x: F(x, zeros), a, b, sing)
    I_low, e2 = cquad(lambda x: F(x - 1j, zeros), a, b)
    sR, e3 = cquad(lambda y: 1j * F(b - 1j * y, zeros), 0, 1)
    sL, e4 = cquad(lambda y: 1j * F(a - 1j * y, zeros), 0, 1)
    S = slits(zeros)
    hp, w1 = 0j, 0j
    det = []
    for x, lst in S.items():
        dep = max(d for (d, _, _) in lst)
        tips = [d for (d, _, _) in lst]
        # merged hairpin: true jump from the two one-sided branches
        def fj(y):
            t = x - 1j * y
            return g(t) * R(t) * Rs(t) * (P(V(t, zeros, +1)) - P(V(t, zeros, -1)))
        val, e = cquad(fj, 0, dep, tips + [0.5 * dep])
        hp += -1j * val
        # jump check at mid-depth of the shallowest tip
        ych = 0.5 * min(tips)
        Jnum = V(x - 1j * ych, zeros, +1) - V(x - 1j * ych, zeros, -1)
        Jth = 2j * sum(sg * m for (d, sg, m) in lst if d > ych)
        det.append(dict(x=x, depth=dep, n=len(lst), J_num=[Jnum.real, Jnum.imag], J_th=[Jth.real, Jth.imag],
                        hairpin=[(-1j * val).real, (-1j * val).imag]))
        # W1: separate increments per zero/sign, each over its own depth, V_L = full left value
        for (d, sg, m) in lst:
            def fw(y, sg=sg, m=m):
                t = x - 1j * y
                VL = V(t, zeros, -1)
                return g(t) * R(t) * Rs(t) * (P(VL + 2j * sg * m) - P(VL))
            v1, _ = cquad(fw, 0, d, [0.5 * d])
            w1 += -1j * v1
    rhs = I_low + sR - sL + hp
    return dict(I_real=[I_real.real, I_real.imag], rhs=[rhs.real, rhs.imag], resid=abs(I_real - rhs),
                resid_W0=abs(I_real - (I_low + sR - sL)), resid_W1=abs(I_real - (I_low + sR - sL + w1)),
                quad_err=e1 + e2 + e3 + e4, hairpins=[hp.real, hp.imag], slits=det)

ON = [(0.5, g0, 1) for g0 in (-7.3, -4.1, -1.9, 0.35, 2.7, 5.2, 8.8)]
CASES = {
    "single":            [(0.75, 0.9, 1)],
    "coinc_opp_equal":   [(0.70, -0.5, 1), (0.70, -0.5 + H, 1)],          # gamma2-gamma1 = h, equal depth: J = 0
    "coinc_opp_unequal": [(0.80, -0.5, 1), (0.62, -0.5 + H, 1)],          # J = 2i(m+ - m-) changes with y
    "same_ordinate":     [(0.85, 1.3, 1), (0.60, 1.3, 1)],                # same-sign coincident, two depths
    "double_offline":    [(0.72, -2.6, 2)],                               # multiplicity 2: J = +-4i
    "online_at_same_ord":[(0.78, 3.6, 1), (0.5, 3.6, 1)],                 # on-line zero at the off-line ordinate
    "online_at_slit_top":[(0.66, 4.4, 1), (0.5, 4.4 - H, 1)],             # on-line shifted ordinate = slit top
    "cluster_mixed":     [(0.80, -0.5, 1), (0.62, -0.5 + H, 1), (0.70, -0.5 + 2 * H, 1), (0.55, -0.5 + H, 1),
                          (0.9, 6.0, 1), (0.53, 6.0, 1)],
    "deep_near_edge":    [(0.97, 10.5, 1), (0.95, -10.9, 1)],
}

if __name__ == "__main__":
    t0 = time.time()
    out = {}
    names = sys.argv[2:] if len(sys.argv) > 2 else list(CASES)
    if names and names[0] == "--mp":
        USE_MP = True
        names = names[1:]
    for name in names:
        offl = CASES[name]
        zeros = []
        for (beta, gam, m) in offl:
            zeros.append((beta, gam, m))
            if beta != 0.5:
                zeros.append((1 - beta, gam, m))
        zeros += ON
        r = run(zeros)
        out[name] = r
        print("%-20s resid=%.2e  noHairpin=%.2e  separateIncr=%.2e  quadErr=%.1e  |I|=%.3g  slits=%d  Jcheck=%s" % (
            name, r['resid'], r['resid_W0'], r['resid_W1'], r['quad_err'], abs(complex(*r['I_real'])), len(r['slits']),
            all(abs(complex(*s['J_num']) - complex(*s['J_th'])) < 1e-9 for s in r['slits'])), flush=True)
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "slit_identity.json", "w"), indent=1)
    print("time %.1f s" % (time.time() - t0))
