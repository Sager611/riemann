#!/usr/bin/env python3
"""dfp_moments.py -- exact rational certificate for the optimisation of de Faveri-Pandey (arXiv 2610.01035).
Part 1 reproduces every printed exponent. Part 2 optimises the free parameters (tau = log_T W, k, the
choice of N-range endpoint, the exponent pair in Ivic's Lemma) and recomputes the moment exponents.
Pure Python, fractions only. Run: python3 dfp_moments.py
"""
from fractions import Fraction as Fr
import sys
F = Fr
OK = True
def check(cond, msg):
    global OK
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond: OK = False

# ---------------------------------------------------------------- the four terms of (eq:curyjytduj), general N
# exponents of T; s = sigma (V = T^s), nu = log_T N, om = log_T W, G = V^2 = T^{2s}
def terms(k, s, nu, om):
    k = F(k)
    t1 = F(1,2) + F(3,2)*nu - 2*s - om
    t2 = 1/(2*k) + (2 - 3/(2*k))*nu - 2*s - (k-1)*om/(2*k)
    t3 = F(3,4) - k/2 + (k/2 + F(1,4))*nu + (k-1)*om/2 - 2*s + max(F(0), nu + om - 1)/4
    t4 = (6 + (16*k-12)*nu - 16*k*s - (4*k-4)*om - 12*s)/(8*k+1)
    return [t1, t2, t3, t4]

def constraints(k, s, nu, om, strict_C6=True):
    """returns list of (name, slack) with slack >= 0 meaning satisfied (C6 needs > 0)."""
    c = []
    c.append(("C1 N<=T^(3/5)", F(3,5) - nu))
    c.append(("C2 W>=N", om - nu))
    c.append(("C3 W<=N^(4/3)", F(4,3)*nu - om))
    c.append(("C4 N^(1/2)W^2<=T^(3/2)", F(3,2) - nu/2 - 2*om))
    c.append(("C5 W<=sqrt(TN)", (1+nu)/2 - om))
    c.append(("C6 N^3/(TW)=o(1) [strict]", 1 + om - 3*nu))
    c.append(("C7a sqrt(T/N)<=V^4", 4*s - (1-nu)/2))
    c.append(("C7b W sqrt(N/T)<=V^4", 4*s - (om + nu/2 - F(1,2))))
    c.append(("C0 N>=T^(1/3)", nu - F(1,3)))
    return c

def feasible(k, s, nu, om):
    for name, sl in constraints(k, s, nu, om):
        if name.startswith("C6"):
            if sl <= 0: return False
        elif sl < 0: return False
    return True

# ---------------------------------------------------------------- PART 1: reproduce the printed exponents
def E_paper(k, s):
    k = F(k)
    if k == 2: return (24 - 128*s)/17
    return (32*k*k - 12*k - 10 - (192*k*k - 120*k - 28)*s)/(16*k*k + 10*k - 17)
def tau_paper(k, s):
    k = F(k)
    if k == 2: return F(1,2) + s
    return 4*((8*k - 7) + (16*k*k - 46*k + 40)*s)/(16*k*k + 10*k - 17)
def sig_m(k):
    k = F(k)
    return F(7,60) if k == 2 else (16*k*k - 22*k + 7)/(8*(16*k*k - 20*k + 5))
def sig_p(k):
    k = F(k)
    return F(5,38) if k == 2 else (4*k - 3)/(30*k - 22)

def part1():
    print("=== PART 1: reproduction of de Faveri-Pandey's printed exponents ===")
    for k in range(2, 13):
        a, b = sig_m(k), sig_p(k)
        check(a < b, f"k={k}: range [{a}, {b}] nonempty")
        for s in (a, b, (a+b)/2):
            nu, om = 1 - 4*s, tau_paper(k, s)
            t = terms(k, s, nu, om)
            E = E_paper(k, s)
            check(t[3] == E, f"k={k} s={s}: T4 == E_k")
            if k >= 3:
                check(t[2] == E, f"k={k} s={s}: T3 == E_k (tau balances T3=T4)")
            else:
                check(t[1] == t[2] == F(11,8) - F(29,4)*s, f"k=2 s={s}: T2 == T3 == 11/8-29s/4")
            check(t[1] <= E, f"k={k} s={s}: T2 <= E_k")
            check(t[0] < E, f"k={k} s={s}: T1 < E_k")
            rho = om/nu
            check(1 <= rho <= F(4,3), f"k={k} s={s}: rho=tau/(1-4s)={float(rho):.5f} in [1,4/3] (terms increase in N)")
            for name, sl in constraints(k, s, nu, om):
                cond = sl > 0 if name.startswith("C6") else sl >= 0
                check(cond, f"k={k} s={s}: {name} slack {sl} ({float(sl):.5f})")
            # small-N end
            nu0, om0 = F(1,3), rho*F(1,3)
            for name, sl in constraints(k, s, nu0, om0):
                cond = sl > 0 if name.startswith("C6") else sl >= 0
                check(cond, f"k={k} s={s} N=T^(1/3): {name} slack {float(sl):.5f}")
        # endpoint identities
        check(E_paper(k, a) == 1 - 4*a, f"k={k}: E_k(sigma_-) = 1-4sigma_- (meets 4th moment)")
        if k >= 3:
            check(E_paper(k, b) == 2 - 12*b, f"k={k}: E_k(sigma_+) = 2-12sigma_+ (meets 12th moment)")
            nu, om = 1 - 4*b, tau_paper(k, b)
            check(terms(k, b, nu, om)[1] == E_paper(k, b), f"k={k}: T2 = E_k exactly at sigma_+")
    # Theorem 1.3
    check(E_paper(2, F(0)) == F(24,17) and E_paper(2, F(1)) - E_paper(2, F(0)) == -F(128,17), "Thm 1.3: N(T,V) << T^(24/17) V^(-128/17) on [T^(7/60),T^(5/38)]")
    # Theorem 1.1 (sixth moment)
    V0, V1 = F(7,60), F(5,38)
    six = max(1 + 2*V0, F(24,17) - F(26,17)*V0, 2 - 6*V1)
    check(six == F(37,30) and F(37,30) == F(5,4) - F(1,60), f"Thm 1.1: sixth moment exponent {six} = 5/4 - 1/60")
    # Theorem 9.1 (twelfth moment)
    check(E_paper(6, F(0)) == F(1070,619) and E_paper(6, F(1)) - E_paper(6, F(0)) == -F(6164,619), "k=6: E_6 = (1070-6164s)/619")
    check(sig_m(6) == F(451,3688) and sig_p(6) == F(21,158), "k=6 range [451/3688, 21/158]")
    m = 10
    kap = F(1,2) - F(3*m-2, m*(m-1)*(m+2)); lam = F(1,2) + F(2, (m-1)**2*(m+2))
    check((kap, lam) == (F(64,135), F(122,243)), f"TY m=10 pair B(p_10,q_10) = ({kap},{lam})")
    A_ty, B_ty = 1 + lam/kap, 2*(1 + 2*kap + 2*lam)/kap
    check((A_ty, B_ty) == (F(593,288), F(3587,288)), f"Ivic/TY bound T^{A_ty} V^-{B_ty}")
    s1 = sig_m(6)
    s2 = (A_ty - F(1070,619))/(F(1264,619) + B_ty - 12)
    check(s2 == F(58907,445121), f"sigma_2 = {s2}")
    val = F(1070,619) + F(1264,619)*s2
    check(val == 2 - F(520,445121) and val == A_ty + (12 - B_ty)*s2, f"twelfth moment exponent {val} = 2 - 520/445121")
    check(s1 <= s2 <= sig_p(6), "sigma_1 <= sigma_2 <= sigma_+(6)")
    check(1 + 8*s1 <= val and 1 + 6*F(13,84) <= val and A_ty + (12-B_ty)*F(13,84) <= val, "other ranges below the max")
    return val


# ---------------------------------------------------------------- PART 2: optimise tau (=om) for each N, k
def om_interval(k, s, nu):
    """feasible om-interval from C2..C7b; None if a om-free constraint fails. C6 is strict: lo_strict flag."""
    for name, sl in constraints(k, s, nu, nu):   # om-free constraints: C1, C7a, C0
        if name[:3] in ("C1 ", "C7a", "C0 ") and sl < 0: return None
    lo, hi = nu, min(F(4,3)*nu, (F(3,2) - nu/2)/2, (1 + nu)/2, 4*s - nu/2 + F(1,2))
    lo6 = 3*nu - 1
    strict = False
    if lo6 >= lo: lo, strict = lo6, True
    if hi < lo or (strict and hi == lo): return None
    return lo, hi, strict

def lin_pieces(k, s, nu):
    """terms as affine functions a + b*om (T3 split at om = 1 - nu)."""
    def aff(f):
        y0 = f(F(0)); y1 = f(F(1)); return (y0, y1 - y0)
    k_ = F(k)
    p = []
    p.append(aff(lambda om: terms(k, s, nu, om)[0]))
    p.append(aff(lambda om: terms(k, s, nu, om)[1]))
    p.append(aff(lambda om: terms(k, s, nu, om)[3]))
    base3 = lambda om: F(3,4) - k_/2 + (k_/2 + F(1,4))*nu + (k_-1)*om/2 - 2*s
    p.append(aff(base3))                                   # T3 when nu+om <= 1
    p.append(aff(lambda om: base3(om) + (nu + om - 1)/4))  # T3 when nu+om >= 1 (max of both = T3)
    return p

def g(k, s, nu):
    """min over feasible om of max_i T_i(nu, om); returns (value, om*, at_C6_boundary) or None."""
    iv = om_interval(k, s, nu)
    if iv is None: return None
    lo, hi, strict = iv
    pcs = lin_pieces(k, s, nu)
    cands = {lo, hi}
    for i in range(len(pcs)):
        for j in range(i+1, len(pcs)):
            (a1, b1), (a2, b2) = pcs[i], pcs[j]
            if b1 != b2:
                x = (a2 - a1)/(b1 - b2)
                if lo <= x <= hi: cands.add(x)
    best = None
    for om in cands:
        v = max(a + b*om for a, b in pcs)
        if best is None or v < best[0] or (v == best[0] and om > best[1]): best = (v, om)
    return best[0], best[1], (strict and best[1] == lo)

def Fk(k, s):
    """optimised large-value exponent for given k: max over the two N-endpoints (g convex in nu)."""
    r1 = g(k, s, F(1,3)); r2 = g(k, s, 1 - 4*s)
    if r1 is None or r2 is None: return None
    return max(r1[0], r2[0]), r2[1], r1[0], (r1[2] or r2[2])

# ---------------------------------------------------------------- exponent pairs (Ivic Thm 8.2 / TY Lemma 5.1)
def hb_pair(m):   # Heath-Brown 2017 (p_m, q_m)
    return F(2, (m-1)**2*(m+2)), 1 - F(3*m-2, m*(m-1)*(m+2))
def Bp(p):  return (p[1] - F(1,2), p[0] + F(1,2))
def Ap(p):  return (p[0]/(2*p[0] + 2), p[1]/(2*p[0] + 2) + F(1,2))
def ivic_line(pair):
    kap, lam = pair
    return 1 + lam/kap, 2*(1 + 2*kap + 2*lam)/kap    # N << T^a V^-b  (+ T V^-6)

def explore():
    print("=== exploration: F_k(s) vs paper E_k(s) and classical ===")
    for s in [F(1,10), F(11,100), F(7,60), F(23,200), F(3,25), F(1,8), F(13,100), F(132339,1000000), F(133,1000), F(135,1000), F(14,100), F(145,1000), F(15,100), F(155,1000)]:
        row = []
        for k in range(2, 16):
            r = Fk(k, s)
            row.append((k, None if r is None else float(r[0])))
        best = min((v, k) for k, v in row if v is not None) if any(v is not None for _, v in row) else None
        print(f"s={float(s):.6f} 1-4s={float(1-4*s):.6f} 2-12s={float(2-12*s):.6f} best={best} ", " ".join(f"{k}:{v:.6f}" if v is not None else f"{k}:-" for k, v in row))

# ---------------------------------------------------------------- PART 3: exponent pairs and moment exponents
def primitive_pairs():
    P = {}
    P["Bourgain (13/84,55/84)"] = (F(13,84), F(55,84))
    for i, p in enumerate([(F(4742,38463), F(35731,51284)), (F(18,199), F(593,796)), (F(2779,38033), F(58699,76066)), (F(715,10238), F(7955,10238))]):
        P[f"TY24 Lemma 1.1 #{i+1}"] = p
    for i, p in enumerate([(F(89,1282), F(997,1282)), (F(652397,9713986), F(7599781,9713986)), (F(10769,351096), F(609317,702192)), (F(89,3478), F(15327,17390))]):
        P[f"TTY25 new-exp-pair #{i+1}"] = p
    for i, p in enumerate([(F(2,13), F(35,52)), (F(516247,6629696), F(5080955,6629696)), (F(6299,43860), F(29507,43860)), (F(771,8116), F(1499,2029)), (F(21,232), F(173,232)), (F(1959,21656), F(16135,21656))]):
        P[f"Huxley/Sargos #{i+1}"] = p
    for m in range(1, 12):
        d = 1424*2**m - 338
        P[f"Huxley T17.3 m={m}"] = (F(169, d), 1 - F(169, d)*F(712*m + 1577, 712))
    for (mm, th) in [(4, F(1,13)), (8, F(1,204)), (9, F(1,360)), (9, F(7,2640)), (10, F(1,615)), (11, F(1,915))]:
        P[f"deriv test m={mm} th={th}"] = (th, 1 - (mm-1)*th)
    for m in range(3, 121):
        P[f"HB17 (p_{m},q_{m})"] = hb_pair(m)
    P["(1/2,1/2)=B(0,1)"] = (F(1,2), F(1,2))
    return P

def Cp(p):  # Sargos C-process (TY24 eq. (C_process_defn)): (k, l+eps) -> (k/(12(1+4k)), (11(1+4k)+l)/(12(1+4k)) + eps)
    k, l = p
    d = 12*(1 + 4*k)
    return (k/d, (11*(1 + 4*k) + l)/d)

def closure(P, depth=3):
    out = dict(P)
    seen = set(P.values())
    frontier = dict(P)
    for d in range(depth):
        new = {}
        for name, p in frontier.items():
            for op, f in (("A", Ap), ("B", Bp), ("C", Cp)):
                if op == "B" and name.startswith("B"):  # B is an involution
                    continue
                q = f(p)
                nm = op + "." + name
                if q not in seen and q[0] > 0:
                    new[nm] = q; seen.add(q)
        out.update(new); frontier = new
    return out

def ivic(pair):
    """Ivic Thm 8.2 (TY24 Lemma 5.1): N(T,T^s) << T^{max(1-6s, a-b s)}."""
    a, b = ivic_line(pair)
    return a, b

def explore2():
    P = closure(primitive_pairs(), depth=4)
    print("pairs:", len(P))
    # float envelope of Ivic lines on s in [0.12, 13/84]
    lines = [(float(ivic(p)[0]), float(ivic(p)[1]), nm, p) for nm, p in P.items()]
    best = {}
    N = 3000
    for i in range(N+1):
        s = 0.12 + (13/84 - 0.12)*i/N
        v, nm = min((max(1-6*s, a - b*s), nm) for a, b, nm, p in lines)
        best.setdefault(nm, []).append(s)
    for nm, ss in sorted(best.items(), key=lambda x: x[1][0]):
        p = P[nm]; a, b = ivic(p)
        print(f"{nm:40s} s in [{ss[0]:.5f},{ss[-1]:.5f}] pair=({float(p[0]):.6f},{float(p[1]):.6f}) a={float(a):.6f} b={float(b):.6f}")
    # twelfth moment with each candidate line vs E_6 (and min_k E_k)
    res = []
    for nm, p in P.items():
        a, b = ivic(p)
        for k in range(3, 12):
            Ek0, Ek1 = E_paper(k, F(0)), E_paper(k, F(1)) - E_paper(k, F(0))
            # 12s + Ek0 + Ek1 s = 12 s + a - b s  -> s* = (a - Ek0)/(Ek1 + b)
            if Ek1 + b == 0: continue
            s = (a - Ek0)/(Ek1 + b)
            if sig_m(k) <= s <= sig_p(k):
                res.append((float(12*s + a - b*s), k, nm, s))
    res.sort()
    for r in res[:12]: print("12th:", r[0], "k=", r[1], r[2], float(r[3]))

# ---------------------------------------------------------------- PART 4: exact moment exponents mu(A) = max(1, sup_s A s + E*(s))
S_MAX = F(13,84)    # Bourgain: |zeta(1/2+it)| << t^{13/84+eps}
class Bound:
    """a valid large-value bound N(T,T^s) << T^{B(s)} on [lo,hi]; B = max of affine lines (c0 + c1 s)."""
    def __init__(self, name, lines, lo=F(0), hi=S_MAX):
        self.name, self.lines, self.lo, self.hi = name, lines, lo, hi
    def val(self, s):
        if s < self.lo or s > self.hi: return None
        return max(c0 + c1*s for c0, c1 in self.lines)

def classical_bounds(pairs):
    B = [Bound("4th moment (Hardy-Littlewood/Ingham): T V^-4", [(F(1), F(-4))]),
         Bound("12th moment (Heath-Brown 1978): T^2 V^-12", [(F(2), F(-12))]),
         Bound("Heath-Brown auxiliary (ANTEDB hb-12-aux): T V^-6 + T^5 V^-32", [(F(1), F(-6)), (F(5), F(-32))])]
    for nm, p in pairs.items():
        a, b = ivic(p)
        B.append(Bound(f"Ivic Thm 8.2 with {nm} = ({p[0]},{p[1]})", [(F(1), F(-6)), (a, -b)]))
    return B

def dfp_bounds(K=60):
    B = []
    for k in range(2, K+1):
        e0 = E_paper(k, F(0)); e1 = E_paper(k, F(1)) - e0
        B.append(Bound(f"DFP Thm 8.1 k={k}: T^({e0}) V^({e1})", [(e0, e1)], sig_m(k), sig_p(k)))
    return B

def Estar(bounds, s):
    vals = [(b.val(s), b) for b in bounds]
    vals = [(v, b) for v, b in vals if v is not None]
    return min(vals, key=lambda x: x[0])

def prune(bounds, A_list=None, grid=3000, chunk=800):
    """keep bounds that attain (within 1e-9) the float envelope somewhere on [0, 13/84] (plus the classical 0,1,2).
    Pruning only removes bounds, so every exact result below is a valid upper bound; optimality over the full set
    is then certified exactly at the maximising vertices (vertex_certificate)."""
    import numpy as np
    s = np.linspace(0.0, 13/84, grid+1)
    def block(i0, i1):
        out = np.empty((i1 - i0, grid+1))
        for j, b in enumerate(bounds[i0:i1]):
            v = np.max(np.array([float(c0) + float(c1)*s for c0, c1 in b.lines]), axis=0)
            v[(s < float(b.lo) - 1e-15) | (s > float(b.hi) + 1e-15)] = np.inf
            out[j] = v
        return out
    env = np.full(grid+1, np.inf); arg = np.zeros(grid+1, dtype=int)
    for i0 in range(0, len(bounds), chunk):
        V = block(i0, min(len(bounds), i0+chunk))
        m = V.min(axis=0); a = V.argmin(axis=0) + i0
        better = m < env
        env = np.where(better, m, env); arg = np.where(better, a, arg)
    keep = set([0, 1, 2]) | set(arg.tolist())
    for i0 in range(0, len(bounds), chunk):
        V = block(i0, min(len(bounds), i0+chunk))
        near = V <= env + 1e-9
        for j in np.nonzero(near.any(axis=1))[0].tolist():
            i = i0 + j; b = bounds[i]
            if len(b.lines) == 1: keep.add(i); continue
            c0, c1 = b.lines[1]
            lin = float(c0) + float(c1)*s
            if np.any(near[j] & (lin >= 1 - 6*s + 1e-12)): keep.add(i)
    return [bounds[i] for i in sorted(keep)]

def sup_moment(bounds, A):
    """exact max over s in [0, 13/84] of A s + E*(s) for the continuous piecewise-linear E* = min of bounds."""
    lines = []
    cands = {F(0), S_MAX}
    for b in bounds:
        cands.add(b.lo); cands.add(b.hi)
        for i in range(len(b.lines)):
            lines.append(b.lines[i])
            for j in range(i+1, len(b.lines)):
                (a0, a1), (c0, c1) = b.lines[i], b.lines[j]
                if a1 != c1: cands.add((c0 - a0)/(a1 - c1))
    L = list(set(lines))
    for i in range(len(L)):
        for j in range(i+1, len(L)):
            (a0, a1), (c0, c1) = L[i], L[j]
            if a1 != c1:
                x = (c0 - a0)/(a1 - c1)
                if 0 <= x <= S_MAX: cands.add(x)
    best = None
    for s in cands:
        if s < 0 or s > S_MAX: continue
        v, b = Estar(bounds, s)
        val = A*s + v
        if best is None or val > best[0]: best = (val, s, b.name)
    return max(best[0], F(1)), best

def before_interp(A):
    """best bound implied by the printed results (M_4, M_6 = 37/30, M_12 = 2-520/445121) and Hoelder (log-convexity)."""
    pts = [(F(4), F(1)), (F(6), F(37,30)), (F(12), 2 - F(520,445121))]
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= A <= x1: return y0 + (A - x0)*(y1 - y0)/(x1 - x0)
    return None

def part4():
    print("=== PART 4: optimised moment exponents ===")
    import time; t0 = time.time()
    pairs = closure(primitive_pairs(), depth=4)
    print("pairs", len(pairs), "t", round(time.time()-t0,1), flush=True)
    CL = classical_bounds(pairs)
    DF = dfp_bounds(60)
    ALL = prune(CL + DF, None)
    CLp = prune(CL, None)
    print("kept bounds (all):", len(ALL), flush=True); [print("   ", b.name) for b in ALL]
    print("kept bounds (classical only):", len(CLp))
    rows = []
    for A in [F(5), F(6), F(13,2), F(7), F(15,2), F(8), F(17,2), F(9), F(19,2), F(10), F(21,2), F(11), F(23,2), F(12)]:
        mu, (val, s, nm) = sup_moment(ALL, A)
        cl, _ = sup_moment(CLp, A)
        bf = before_interp(A)
        rows.append((A, mu, s, nm, bf, cl))
        print(f"A={A}: mu*={mu} ({float(mu):.7f}) at s*={s} ({float(s):.6f}) [{nm}]; before(printed+Hoelder)={float(bf):.7f}; classical={float(cl):.7f}; delta*={float(1+(A-4)/8-mu):.7f} vs {float(1+(A-4)/8-bf):.7f}", flush=True)
    return rows

def vertices(bounds):
    """exact candidate breakpoints s_j of E* and values e_j = E*(s_j)."""
    lines, cands = [], {F(0), S_MAX}
    for b in bounds:
        cands.add(b.lo); cands.add(b.hi); lines += b.lines
    L = list(set(lines))
    for i in range(len(L)):
        for j in range(i+1, len(L)):
            (a0, a1), (c0, c1) = L[i], L[j]
            if a1 != c1:
                x = (c0 - a0)/(a1 - c1)
                if 0 <= x <= S_MAX: cands.add(x)
    out = []
    for s in sorted(cands):
        v, b = Estar(bounds, s); out.append((s, v, b.name))
    return out

def curve(V, A_lo, A_hi, n=400):
    """maximising vertices of A s + e over A in [A_lo, A_hi], with exact switch points."""
    seq = []
    for i in range(n+1):
        A = A_lo + (A_hi - A_lo)*F(i, n)
        j = max(range(len(V)), key=lambda t: (A*V[t][0] + V[t][1], V[t][0]))
        if not seq or seq[-1] != j: seq.append(j)
    pieces = []
    for a, b in zip(seq, seq[1:]):
        (s1, e1, _), (s2, e2, _) = V[a], V[b]
        pieces.append((e1 - e2)/(s2 - s1))      # A where A s1 + e1 = A s2 + e2
    return seq, pieces

def part5(V=None, Vc=None):
    print("=== PART 5: the curve mu*(A), 4 < A <= 12, and the range above 12 ===", flush=True)
    if V is None:
        pairs = closure(primitive_pairs(), depth=4)
        CL = classical_bounds(pairs); DF = dfp_bounds(60)
        ALL = prune(CL + DF); CLp = prune(CL)
        V = vertices(ALL); Vc = vertices(CLp)
    seq, sw = curve(V, F(4), F(14))
    print("mu*(A) = max_j (A s_j + e_j); maximising vertices for A in [4,14]:")
    prev = F(4)
    for t, j in enumerate(seq):
        s, e, nm = V[j]
        nxt = sw[t] if t < len(sw) else F(14)
        print(f"  A in [{prev}, {nxt}] ~ [{float(prev):.5f},{float(nxt):.5f}]: s_j={s} ({float(s):.6f}), e_j=E*(s_j)={e} ({float(e):.6f}) [{nm}]", flush=True)
        prev = nxt
    seqc, swc = curve(Vc, F(12), F(14))
    print("classical mu_cl(A) for A in [12,14]:")
    prev = F(12)
    for t, j in enumerate(seqc):
        s, e, nm = Vc[j]
        nxt = swc[t] if t < len(swc) else F(14)
        print(f"  A in [{float(prev):.5f},{float(nxt):.5f}]: s={float(s):.6f} e={float(e):.6f} [{nm}]", flush=True)
        prev = nxt
    # crossing: smallest A > 12 with mu*(A) >= mu_cl(A)
    def mu_from(VV, A): return max(A*s + e for s, e, _ in VV)
    lo, hi = F(12), F(14)
    assert mu_from(V, lo) < mu_from(Vc, lo)
    for _ in range(60):
        mid = (lo + hi)/2
        if mu_from(V, mid) < mu_from(Vc, mid): lo = mid
        else: hi = mid
    # exact: identify the linear pieces at the crossing and solve
    j1 = max(range(len(V)), key=lambda t: lo*V[t][0] + V[t][1]); j2 = max(range(len(Vc)), key=lambda t: lo*Vc[t][0] + Vc[t][1])
    (s1, e1, n1), (s2, e2, n2) = V[j1], Vc[j2]
    Ax = (e2 - e1)/(s1 - s2)
    ok = mu_from(V, Ax) == mu_from(Vc, Ax) == Ax*s1 + e1
    print(f"range: DFP-based bound beats every classical bound for 12 <= A < {Ax} = 12 + {Ax-12} ({float(Ax):.6f}); exact check {ok}")
    print(f"   DFP side vertex s={s1} [{n1}]; classical side vertex s={s2} [{n2}]")
    return Ax

# ---------------------------------------------------------------- PART 6: refinements, sensitivity, other inputs
def pieces6(k, s, nu, eta, cap):
    """affine pieces (a, b) in om of all terms; T3 and the capped T2 split at om = 1 - nu."""
    k_ = F(k)
    def aff(f):
        y0 = f(F(0)); return (y0, f(F(1)) - y0)
    t = lambda om: terms(k, s, nu, om)
    base3 = lambda om: F(3,4) - k_/2 + (k_/2 + F(1,4))*nu + (k_-1)*om/2 - 2*s
    P = [aff(lambda om: t(om)[0]), aff(lambda om: t(om)[1]), aff(base3), aff(lambda om: base3(om) + (nu + om - 1)/4),
         aff(lambda om: t(om)[3] - 2*eta/(8*k_ + 1))]
    if cap: P.append(aff(lambda om: t(om)[1] + (nu + om - 1)/(4*k_)))
    return P

def obj6(k, s, nu, om, eta, cap):
    """exponent of R at (N, W) = (T^nu, T^om): eta = saving T^-eta in Iwaniec's Prop 6.1 (X-term T4);
    cap: sup over D <= Upsilon of min((2D+1)R, R(1+X'))/sqrt(D) = R sqrt(min(1+X', Upsilon)) (middle term)."""
    t1, t2, t3, t4 = terms(k, s, nu, om)
    t4 = t4 - 2*eta/(8*k + 1)
    mid = max(t2, t4)
    if cap: mid = min(mid, t2 + max(F(0), nu + om - 1)/(4*k))
    return max(t1, t3, mid)

def g6(k, s, nu, eta=F(0), cap=False):
    """exact min over feasible om (continuous piecewise-linear objective: candidates = endpoints, kinks, crossings)."""
    iv = om_interval(k, s, nu)
    if iv is None: return None
    lo, hi, strict = iv
    P = pieces6(k, s, nu, eta, cap)
    cands = {lo, hi, 1 - nu}
    for i in range(len(P)):
        for j in range(i+1, len(P)):
            (a1, b1), (a2, b2) = P[i], P[j]
            if b1 != b2: cands.add((a2 - a1)/(b1 - b2))
    best = None
    for om in cands:
        if lo <= om <= hi:
            v = obj6(k, s, nu, om, eta, cap)
            if best is None or v < best[0]: best = (v, om)
    return best

def F6(k, s, eta=F(0)):
    """valid exponent (cap off): objective jointly convex in (nu, om) => min over om convex in nu => max at endpoints."""
    a = g6(k, s, F(1,3), eta); b = g6(k, s, 1 - 4*s, eta)
    if a is None or b is None: return None
    return max(a[0], b[0])

def part6(verts, ALL):
    print("=== PART 6: refinements, sensitivity, other inputs ===", flush=True)
    # (a) N-dependence and tau: at every curve vertex the LP over W reproduces E_k, worst N is the top N = T^(1-4s)
    for s, e, nm in verts:
        for k in range(2, 13):
            if not (sig_m(k) <= s <= sig_p(k)): continue
            top, bot = g6(k, s, 1 - 4*s), g6(k, s, F(1,3))
            check(top[0] == E_paper(k, s) and bot[0] <= top[0], f"(a) s={s} k={k}: min_W bound at N=T^(1-4s) = E_k = {top[0]} (W=T^{top[1]}); at N=T^(1/3): {float(bot[0]):.6f} <= it")
            cap = g6(k, s, 1 - 4*s, cap=True)
            check(cap[0] == E_paper(k, s), f"(b) s={s} k={k}: Upsilon-capped middle term gives {cap[0]} = E_k (no gain)")
    # (b') Upsilon cap on a grid inside every range
    worst = None
    for k in range(2, 13):
        a, b = sig_m(k), sig_p(k)
        for i in range(25):
            s = a + (b - a)*F(i, 24)
            cap = g6(k, s, 1 - 4*s, cap=True)
            d = E_paper(k, s) - cap[0]
            worst = d if worst is None or d > worst else worst
    check(worst == 0, f"(b') Upsilon cap, 11 x 25 exact grid points: max gain {worst}")
    # (c) using the chain outside the stated ranges [sigma_-(k), sigma_+(k)] (k = 2..20, exact grid on [3/26, 13/84])
    gain = F(0); where = None
    for i in range(1, 121):
        s = F(3,26) + (S_MAX - F(3,26))*F(i, 120)
        ek = Estar(ALL, s)[0]
        for k in range(2, 21):
            v = F6(k, s)
            if v is not None and ek - v > gain: gain, where = ek - v, (s, k)
    check(gain == 0, f"(c) chain outside its stated ranges: max gain over the envelope {gain} {where}")
    # (d) sensitivity to the X-term input (Iwaniec, Prop 6.1 with an extra factor T^-eta)
    for eta in (F(1,100), F(1,50), F(1,30)):
        st = (7 - 2*eta)/60
        v = F6(2, st, eta)
        check(v == 1 - 4*st and v == (24 - 128*st - 2*eta)/17, f"(d) eta={eta}: k=2 bound at s*=(7-2eta)/60={st} equals 1-4s* -> sixth moment 1+2s* = {1+2*st} = 37/30 - eta/15")
    # twelfth moment: crossing of min_k (chain with eta) with the classical envelope; X-term removed entirely as eta = 1
    CLx = [b for b in ALL if not b.name.startswith("DFP")]
    for eta in (F(0), F(1,100), F(1,26), F(1)):
        def f(s):
            fk = min(v for v in (F6(k, s, eta) for k in range(3, 13)) if v is not None)
            return fk - Estar(CLx, s)[0], fk
        lo, hi = F(13,100), F(134,1000)
        assert f(lo)[0] < 0 < f(hi)[0]
        for _ in range(45):
            mid = (lo + hi)/2
            if f(mid)[0] < 0: lo = mid
            else: hi = mid
        tw = 12*hi + f(hi)[1]
        ks = [k for k in range(3, 13) if F6(k, hi, eta) == f(hi)[1]]
        t = terms(ks[0], hi, 1 - 4*hi, g6(ks[0], hi, 1 - 4*hi, eta)[1])
        print(f"(d) eta={eta}: twelfth moment exponent {float(tw):.7f} (crossing s={float(hi):.7f}, k={ks}, terms T1..T4 at top N: {[round(float(x),6) for x in t[:3]] + [round(float(t[3] - 2*eta/(8*ks[0]+1)),6)]})")
        if eta == 0: check(abs(tw - (2 - F(520,445121))) < F(1, 10**9), "(d) eta=0 reproduces 2 - 520/445121")
    for eta in (F(0), F(1,100), F(1,26), F(1,7)):
        for s in (F(1,9), F(3,26), F(1,8), F(5,38)):
            dd = (24 - 128*s - 2*eta)/17 - (F(11,8) - F(29,4)*s)
            check(dd == (5 - 38*s - 16*eta)/136, f"(d) algebra eta={eta} s={s}: E_2^eta - (11/8 - 29s/4) = (5 - 38s - 16eta)/136")
    check((7 - 2*F(1,26))/60 == (5 - 16*F(1,26))/38 == F(3,26), "(d) s*=(7-2eta)/60 <= (5-16eta)/38 iff eta <= 1/26 (equal at 1/26: 3/26)")
    s0 = F(3,26)
    check(F(1,2) + s0 == 2 - 12*s0 == F(8,13), "(d) floor: at s=3/26, C4 cap tau=1/2+s equals C6 limit tau=2-12s (=8/13): chain empty for s<=3/26 at N=T^(1-4s)")
    t = terms(2, s0, 1 - 4*s0, F(8,13))
    check(t[1] == t[2] == 1 - 4*s0 == F(7,13) and 1 + 2*s0 == F(16,13), "(d) floor: there T2=T3=1-4s=7/13; sixth moment floor 1+2s = 16/13 (eta = 1/26)")
    # (e) Dirichlet-polynomial inputs through the approximate functional equation at nu = 1/2 (worst piece), s = 1/8
    s, nu = F(1,8), F(1,2)
    gm, hux = 1 + F(2,5)*nu - 4*s, 1 + nu - 6*s
    check(gm > 1 - 4*s and hux > 1 - 4*s, f"(e) Guth-Maynard third term T N^(12/5-4s') -> T^(1+2nu/5-4s) = T^{gm}, Huxley T N^(4-6s') -> T^(1+nu-6s) = T^{hux}; both > T^(1-4s) = T^{1-4*s}")

def vertex_certificate(verts, FULL):
    """exact: at every maximising vertex, no bound of the full set (and no tail k > 60, m > 120) is below the kept envelope."""
    print("=== vertex certificate: optimality over the full input set ===", flush=True)
    for s, e, nm in verts:
        bad = [b.name for b in FULL if b.val(s) is not None and b.val(s) < e]
        g = (2 - 12*s) - e
        check(not bad, f"s={s}: all {len(FULL)} bounds >= E*(s)={e} [{nm}]" + (f" VIOLATED by {bad[:3]}" if bad else ""))
        K = max(61, int((4 - 30*s)/(2*g)) + 2) if g > 0 else 61
        badk = [k for k in range(61, K+1) if sig_m(k) <= s <= sig_p(k) and E_paper(k, s) < e]
        check(g > 0 and s <= F(3,22) and not badk, f"s={s}: DFP k=61..{K} none below; k>{K}: 2-12s-E_k <= (4-30s)/(2k) < gap {float(g):.6f}")
        if s > F(1,8):
            tail = 6*(8*s - 1)/(121**2 - 6)
            check(tail <= g, f"s={s}: Ivic with B(p_m,q_m), m>120: v_m >= 2-12s-6(8s-1)/(m^2-6), deficit {float(tail):.2e} <= gap")
        else:
            check(True, f"s={s} <= 1/8: every pair (lambda >= 1/2) gives Ivic line >= 1-4s >= E*(s)")

def run_all():
    import time; t0 = time.time()
    part1()
    pairs = closure(primitive_pairs(), depth=4)
    CL = classical_bounds(pairs); DF = dfp_bounds(60)
    FULL = CL + DF
    ALL = prune(FULL); CLp = prune(CL)
    print(f"exponent pairs {len(pairs)}, bounds {len(FULL)}, kept {len(ALL)} (classical kept {len(CLp)}), t={time.time()-t0:.1f}s", flush=True)
    for b in ALL: print("   kept:", b.name)
    V = vertices(ALL); Vc = vertices(CLp)
    print(f"envelope candidate vertices: {len(V)} (classical {len(Vc)}), t={time.time()-t0:.1f}s", flush=True)
    print("=== TABLE: mu*(A) = max(1, sup_s A s + E*(s)); before = DFP printed (4th, 37/30, 2-520/445121) + Hoelder ===")
    rows = []
    for A in [F(n, 2) for n in range(9, 25)]:
        j = max(range(len(V)), key=lambda t: A*V[t][0] + V[t][1])
        s, e, nm = V[j]; mu = max(F(1), A*s + e)
        cl = max(F(1), max(A*x + y for x, y, _ in Vc))
        bf = before_interp(A)
        rows.append((A, mu, bf, cl))
        print(f"A={A}: mu*={mu} = {float(mu):.7f} (s*={s}, {nm}); before={float(bf):.7f}; classical={float(cl):.7f}; delta*={float(1+(A-4)/8-mu):.7f} vs {float(1+(A-4)/8-bf):.7f}", flush=True)
    seq, sw = curve(V, F(4), F(12))
    verts = [V[j] for j in seq]
    print("=== CURVE on [4,12]: mu*(A) = A s_j + e_j on the A-intervals below ===")
    prev = F(4)
    for t, j in enumerate(seq):
        s, e, nm = V[j]; nxt = sw[t] if t < len(sw) else F(12)
        print(f"  A in [{prev}, {nxt}]: mu*(A) = {s} A + {e}  [{nm}]")
        prev = nxt
    # gap before - mu* is concave on [6,12] (before linear, mu* convex): its max is at a breakpoint of mu*
    gaps = []
    for t in range(len(sw)):
        A = sw[t]
        if 6 < A < 12:
            j = max(range(len(V)), key=lambda u: A*V[u][0] + V[u][1])
            mu = A*V[j][0] + V[j][1]; gaps.append((before_interp(A) - mu, A, mu))
            print(f"breakpoint A={A}: mu*={mu} ({float(mu):.7f}), before={float(before_interp(A)):.7f}, gap={float(before_interp(A)-mu):.7f}")
    g, A, mu = max(gaps)
    check(A == F(128,17) and mu == F(24,17), f"largest gain {float(g):.7f} at A=128/17, where mu*=24/17")
    # Hoelder/convexity cross-checks
    for A, mu, bf, cl in rows:
        if A <= 6: check(mu == bf, f"A={A}: mu* equals Hoelder(4, 6) (no gain for A <= 6)")
        elif A < 12: check(mu < bf, f"A={A}: mu* < before by {float(bf-mu):.7f}")
        else: check(mu == bf, f"A={A}: mu* = printed twelfth moment")
    part5(V, Vc)
    a10, b10 = ivic(Bp(hb_pair(10)))
    mu12 = 2 - F(520,445121)
    s2 = F(58907,445121)
    check(mu12 + (b10 - 12)*s2 == a10, "for 12 <= A <= 3587/288 the bound (58907 A + 182838)/445121 is Hoelder between M_12 (printed) and M_{3587/288} <= 593/288 (Ivic, B(p_10,q_10))")
    vertex_certificate(verts + [(s2, a10 - b10*s2, "range end")], FULL)
    part6(verts, ALL)
    print(f"total time {time.time()-t0:.1f}s")
    print("ALL PASS" if OK else "SOME FAIL")

if __name__ == "__main__":
    if "--part1" in sys.argv:
        part1(); print("ALL PASS" if OK else "SOME FAIL")
    if "--explore" in sys.argv:
        explore()
    if "--explore2" in sys.argv:
        explore2()
    if "--part4" in sys.argv:
        part4()
    if "--part5" in sys.argv:
        part5()
    if "--all" in sys.argv:
        run_all()
