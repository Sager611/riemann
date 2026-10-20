#!/usr/bin/env python3
"""Do de Faveri-Pandey (arXiv 2610.01035), Guth-Maynard (arXiv 2405.20552) or Heath-Brown's twelfth moment
lower the refined quasi-RH exponent B_flat = 0.874886647040 (Liu, arXiv 2610.26835, Sec. 14)?

The three inputs are t-aspect statements about zeta alone. The refined system needs bounds for the sextic
conductor family {psi_u : u sixth-power-free, N(u) ~ U} over F, at heights <= Z^tau with tau -> 0
(OpenAI 7/8 paper, Lemma 8.1). So no published input fills a slot. Below, each input is plugged in as a
HYPOTHETICAL transfer T -> U at its exact published exponent (optimistically: the incumbent's numerator
shift theta is kept for every new numerator route), and the optimisation is rerun.

Slots: numerator first moment over a class B, |B| = U^R: K = (3R+1)/4 (Liu (14.5), Xu-Sang);
       class count R = R_kappa(delta,x) (Liu (14.3)), at most the family zero count U^{A(a)(1-a)}, a=(1+delta)/2.
Run: python3 new_inputs_check.py        (numpy only; about 100 s, about 55 MB)
"""
import math
import time
import numpy as np

T0 = time.time()
# ----------------------------------------------------------------------------------------------------------
# Liu (14.3) and (14.8), re-implemented from the paper (not imported).
def R143(d, x, k):
    c = 1/(3*k); D = 3-(1+2*c)*x; P = (2-2*c*x)*(1-x); J = (5/6-d)*D+d*P
    return 1-d+(5/6-d)*d*P/(2*J)

def quartic_root():
    rts = np.roots([133488, -376596, 384534, -167541, 26231])
    return [r.real for r in rts if abs(r.imag) < 1e-12 and 0.87 < r.real < 0.88][0]

# ---------------------------------------------- exact published exponents --------------------------------
ETA6 = 1/4-1/60                       # DFP Thm 1.1:  M6(T)  << T^{5/4-1/60}
ETA12 = 1-520/445121                  # DFP Thm 10.1: M12(T) << T^{2-520/445121}  (Heath-Brown 1978: 2)
def eta(k):                           # excess of M_k: convex envelope of (4,0), (6,ETA6), (12,ETA12)
    return (k-4)/2*ETA6 if k <= 6 else ETA6+(k-6)/6*(ETA12-ETA6)
KS = (4.5, 5, 5.5, 6, 7, 8, 9, 10, 11, 12)

# Large-value exponent n(v) of N(T,T^v) on [0, 13/84] (Bourgain's pointwise bound: -inf beyond, optimistic):
# min of pieces A - S v on [lo,hi]: (2) [4th, 12th moment] and DFP Thm 9.1 (k=2 is Thm 1.3).
VMAX = 13/84
PIECES = [(0.0, VMAX, 1.0, 4.0), (0.0, VMAX, 2.0, 12.0), (7/60, 5/38, 24/17, 128/17)]
for k in range(3, 101):
    q = 16*k*k+10*k-17
    PIECES.append(((16*k*k-22*k+7)/(8*(16*k*k-20*k+5)), (4*k-3)/(30*k-22),
                   (32*k*k-12*k-10)/q, (192*k*k-120*k-28)/q))
PLO, PHI, PA, PS = (np.array(t) for t in zip(*PIECES))

def n_at(v, right=False):             # n at v (closed ranges) or its right limit
    act = (PLO <= v) & ((PHI > v) if right else (PHI >= v))
    return np.min(np.where(act, PA-PS*v, np.inf))

BPTS = np.unique(np.concatenate([PLO, PHI, [0.0, VMAX]]))
BVAL = np.array([[n_at(b), n_at(b, True) if b < VMAX else -np.inf] for b in BPTS])

def K_lv(R):                          # exact max_v [v + min(R, n(v))] (level-set first moment exponent)
    best = np.max(BPTS[:, None]+np.minimum(R, BVAL))             # breakpoints and right limits
    vx = (PA-R)/PS                                                # crossings n_p(v) = R inside ranges
    for v in vx[(vx >= PLO) & (vx <= PHI)]:
        best = max(best, v+min(R, n_at(v)))
    return best

RTAB = np.linspace(0, 1, 4001)
KLV = np.array([K_lv(r) for r in RTAB])
RLV = 8/15     # for R >= 1-4*(7/60) the crossing lies below 7/60, where n = 1-4v: K_lv = (3R+1)/4 exactly

TTY = [(0.5, 0.7, lambda s: 3/(2-s)), (0.7, 19/25, lambda s: 15/(3+5*s)),          # GM Thm 1.2
       (19/25, 127/167, lambda s: 9/(8*s-2)), (127/167, 13/17, lambda s: 15/(13*s-3)),
       (13/17, 17/22, lambda s: 6/(5*s-1)), (17/22, 41/53, lambda s: 2/(9*s-6)),
       (41/53, 7/9, lambda s: 9/(7*s-1)), (7/9, 1867/2347, lambda s: 9/(8*(2*s-1))),
       (1867/2347, 7/8, lambda s: 3/(2*s)), (7/8, 279/314, lambda s: 3/(10*s-7)),
       (279/314, 9/10, lambda s: 24/(30*s-11)), (9/10, 31/34, lambda s: 3/(10*s-7)),
       (31/34, 14/15, lambda s: 11/(48*s-36))]
def A_best(s):                        # best known zeta density exponent (Tao-Trudgian-Yang Table 2)
    return min(f(s) for lo, hi, f in TTY if lo <= s <= hi)

# ------------------------------------------- A. independent baseline --------------------------------------
def beta_b(d, k):                     # d E / d b at x = 1/2  (E is affine in b)
    return 1/6-d/(16*3*k)-3*(1-R143(d, 0.5, k))/8

def E_half(d, b, l, k, R=None):       # Liu (14.8) at x = 1/2
    c = 1/(3*k); h = (1+b+3*l)/2; ly = (1+b-l)/2
    R = R143(d, 0.5, k) if R is None else R
    return -1/4+b/6+5*l/4+d/2+(d/2)*(l+c*(h/4-ly/2))-(3*h/4)*(1-R)

def bisect(f, lo, hi, n=200):
    flo = f(lo)
    for _ in range(n):
        mid = (lo+hi)/2; fm = f(mid)
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi = mid
    return (lo+hi)/2

def tangent(k):                       # critical class delta_c(kappa), l*(kappa), b*(kappa)
    dc = bisect(lambda d: beta_b(d, k), 0.3, 0.6)
    e0 = E_half(dc, 0, 0, k); e1 = E_half(dc, 0, 1, k)-e0
    ls = -e0/e1
    hd = 1e-6
    dE = (E_half(dc+hd, 0, ls, k)-E_half(dc-hd, 0, ls, k))/(2*hd)
    dB = (beta_b(dc+hd, k)-beta_b(dc-hd, k))/(2*hd)
    return dc, ls, -dE/dB

kf = bisect(lambda k: 11/12-tangent(k)[1]/4-(1+k)/2, 0.70, 0.76)
dc, lstar, bstar = tangent(kf)
Bt = 11/12-lstar/4; Bq = quartic_root()
print(f"A) tangent fixed point: B={Bt:.12f}  quartic root B_flat={Bq:.12f}  diff={Bt-Bq:+.1e}")
print(f"   kappa={kf:.9f} (2B-1={2*Bq-1:.9f})  delta_c={dc:.6f}  l*={lstar:.6f}  b*={bstar:.6f}")

# full-rectangle check of the tangent choice (Liu Lemma 14.10), x in [0,1/2], delta in [1/8,5/6]
dg = np.linspace(1/8, 5/6, 1201); xg = np.linspace(0, 0.5, 241)
DD, XX = np.meshgrid(dg, xg, indexing="ij")

def E_full(d, x, b, l, k, R=None, K=None):    # Liu (14.8) on the rectangle
    c = 1/(3*k); h = (1+b+3*l)/2; ly = (1+b-l)/2
    R = R143(d, x, k) if R is None else R
    K = (3*R+1)/4 if K is None else K
    return -1/4+b/6+5*l/4+d/2+d*x*l+d*x*c*(h/4-ly/2)-h+h*K   # = C0+d/2+dx(l+c(h/4-ly/2))-(3h/4)(1-R)

Eb = E_full(DD, XX, bstar, lstar, kf)
i = np.unravel_index(Eb.argmax(), Eb.shape)
print(f"   max E over rectangle at (b*,l*) = {Eb.max():+.2e} at delta={dg[i[0]]:.4f}, x={xg[i[1]]:.3f}")

# ------------------------------------------- B. critical-class bookkeeping --------------------------------
a = (1+dc)/2; Rc = R143(dc, 0.5, kf); K4 = (3*Rc+1)/4
h = (1+bstar+3*lstar)/2
print(f"B) critical class: delta={dc:.5f} a={a:.5f} x=1/2  R={Rc:.5f}  K4=(3R+1)/4={K4:.5f}  h={h:.5f}")
print("   numerator slot, k-th moment Hölder K_k=((k-1)R+1+eta_k)/k beats K4 iff eta_k < (k-4)(1-R)/4:")
for kk in (5, 6, 8, 12):
    need = (kk-4)*(1-Rc)/4; Kk = ((kk-1)*Rc+1+eta(kk))/kk
    print(f"     k={kk:>2}: need eta<{need:.4f}, DFP eta={eta(kk):.4f} (HB/interp {(kk-4)/8:.4f}) -> K_k-K4={Kk-K4:+.5f}")
vstar = (1-Rc)/4
print(f"   numerator slot, level sets: critical level v*=(1-R)/4={vstar:.5f} < 7/60={7/60:.5f} (DFP range start);"
      f" exact K_LV-K4={K_lv(Rc)-K4:+.1e}")
Aneed = 2*Rc/(1-dc)
print(f"   count slot, zero density U^(A(a)(1-a)) beats R iff A(a) < 2R/(1-delta) = {Aneed:.5f}; known at a={a:.4f}:")
for name, A in (("Ingham 3/(2-s)", 3/(2-a)), ("Guth-Maynard 30/13", 30/13),
                ("Guth-Maynard 15/(3+5s) = TTY best", A_best(a)), ("Chen-Gupta-Li (chars mod q) 7/3", 7/3)):
    print(f"     {name:<36} A={A:.5f}  count exponent {A*(1-a):.5f} vs R={Rc:.5f}")

# ------------------------------------------- C. reruns with hypothetical transfers ------------------------
RD = np.array([A_best((1+d)/2)*(1-(1+d)/2) for d in dg])

def E_scen(d, x, b, l, k, scen):
    R = R143(d, x, k)
    if "count" in scen:
        R = np.minimum(R, np.interp(d, dg, RD))
    K = (3*R+1)/4
    if "mom" in scen:
        for kk in KS:
            K = np.minimum(K, ((kk-1)*R+1+eta(kk))/kk)
    if "lv" in scen:
        K = np.minimum(K, np.where(R >= RLV, (3*R+1)/4, np.interp(R, RTAB, KLV)))
    return E_full(d, x, b, l, k, R=R, K=K)

def maxE(b, l, k, scen):
    E = E_scen(DD, XX, b, l, k, scen)
    j = np.unravel_index(E.argmax(), E.shape)
    xs = xg[j[1]]; lo, hi = max(1/8, dg[j[0]]-0.002), min(5/6, dg[j[0]]+0.002)
    f = lambda d: float(E_scen(d, xs, b, l, k, scen)); g = (math.sqrt(5)-1)/2
    for _ in range(50):
        m1, m2 = hi-g*(hi-lo), lo+g*(hi-lo)
        if f(m1) > f(m2): hi = m2
        else: lo = m1
    return max(E.max(), f((lo+hi)/2))

def golden_min(f, a0, b0, tol=1e-8):
    g = (math.sqrt(5)-1)/2; c = b0-g*(b0-a0); d = a0+g*(b0-a0); fc, fd = f(c), f(d)
    while b0-a0 > tol:
        if fc < fd: b0, d, fd = d, c, fc; c = b0-g*(b0-a0); fc = f(c)
        else: a0, c, fc = c, d, fd; d = a0+g*(b0-a0); fd = f(d)
    return (a0+b0)/2, min(fc, fd)

def solve(scen, k=None):
    k = 2*Bq-1 if k is None else k
    B, lo = 11/12-0.160/4, 0.160
    for _ in range(2):                                   # licensing fixed point kappa = 2B-1
        lo, hi = 0.160, 0.200
        for _ in range(32):
            l = (lo+hi)/2
            _, v = golden_min(lambda b: maxE(b, l, k, scen), 0.0, 0.2, 1e-7)
            lo, hi = (l, hi) if v <= 0 else (lo, l)
        B = 11/12-lo/4; k = 2*B-1
    return B, lo

print("C) reruns (balanced lengths, b free, l bisected, kappa=2B-1):")
for tag, scen in (("baseline (Liu Sec. 14)", ()),
                  ("numerator: DFP/HB k-th moments", ("mom",)),
                  ("numerator: DFP large values (Thm 9.1)", ("lv",)),
                  ("count: GM/TTY zero density", ("count",)),
                  ("all of the above at once", ("mom", "lv", "count"))):
    B, l = solve(scen)
    print(f"   {tag:<40} B={B:.9f}  (B-B_flat={B-Bq:+.1e})  l={l:.6f}  [{time.time()-T0:.0f}s]", flush=True)
