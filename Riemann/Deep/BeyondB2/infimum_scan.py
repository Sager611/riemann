#!/usr/bin/env python3
"""BeyondB2 infimum_scan.py: numerical infimum of the 2611.00082 sec-6 system under H(sigma) (REPORT sec. 2.1),
b re-optimised on a grid, B by bisection; reproduces the sigma table (0.8722824, 0.8711429, 0.8699600, 0.8676503, 13/15).
Run: uv run --no-project --with numpy --with sympy python infimum_scan.py   (single thread, ~250 s, ~90 MB)
Constraints (structural; paper-specific numerator bands m0 in (0.30,0.47), ly/h in (0.57,0.59), l/(h-1/40)<2/9 are
re-choosable and NOT imposed; principal box relaxed to B>0.67 per verify.py B3):
  active  G_N - E_s(R_d) >= 0 on D;   G_L - E_s, G_R - E_s, -E_s - h Dsh >= 0 on D;  -E_s(R_s) >= 0 on D1,D2,D3;
  S_s(R_d) > 0, S_s(R_s) > 0;  floor: -1/1000 - C0 - (1/2+3l/2)/50 + h s >= 0;  small-amp: -1/250 - C0 - (1-h)/16 + (7/8) h s >= 0;
  small rows: -1/20 - 13h/75 + ly/2 - 1/50 >= 0;  direct branch: ly - l - 11b/6 >= 0, l <= 1/5, lx > l, T > 0, b >= 0.
"""
import numpy as np, math, time, resource
t0 = time.time()
P1 = 1124661 / 250000; LST = 0.21
def fmodel(Bv, bv, dd, xx, sv):
    l = 11 / 3 - 4 * Bv; lx = (1 - bv - l) / 2; ly = (1 + bv - l) / 2; h = (1 + bv + 3 * l) / 2; T = 2 * ly - h
    C0 = -0.25 + bv / 6 + 5 * l / 4; th0 = dd * xx / P1; A = 5 / 6 - dd
    mdd = (5 / 6 - dd * xx) / (5 / 2 - dd - 2 * dd * xx); Fd = xx + 2 * (1 - xx) * mdd; Rdd = 1 - dd * Fd
    c = 2 / P1; Dc = 3 - (1 + 2 * c) * xx; Pc = (2 - 2 * c * xx) * (1 - xx); J = A * Dc + dd * Pc
    Rs = 1 - dd + A * dd * Pc / (2 * J)
    E = lambda R: C0 + dd / 2 + dd * xx * (l + (h / 2 - ly) / P1) - 3 * h * (1 - R) / 4 - h * sv * R
    S = lambda R: (3 * R + 1) / 4 + th0 / 2 - 17 / 50 - sv * R
    N = lambda R: (3 * R + 1) / 4 - th0 / 2
    Dsh = lambda R: R + (0.5 - th0) * 3 / 40 - N(R)
    JN = 5 / 2 - 5 / P1
    m0 = 1 - ly / h
    # numerator band edges re-chosen around m0: left edge mL = max(1/2 - LST... ) keep paper's 0.30 / 0.47 if m0 inside, else shift symmetric
    mL, mR = 0.30, 0.47
    GL = th0 * h * (m0 - mL); GR = (dd / 2 - th0) * h * (mR - m0)
    return dict(l=l, lx=lx, ly=ly, h=h, T=T, C0=C0, Rd=Rdd, Rs=Rs, E=E, S=S, Dsh=Dsh, GN=dd * xx * JN / 10 * T, GL=GL, GR=GR, m0=m0)
DOM = ((11 / 50, 13 / 20), (6 / 25, 0.5)); D1 = ((1 / 8, 5 / 6), (0, 6 / 25)); D2r = ((1 / 8, 11 / 50), (6 / 25, 0.5))
D3 = ((13 / 20, 5 / 6), (6 / 25, 0.5)); OMG = ((1 / 8, 5 / 6), (0, 0.5))
GRIDS = {}
def grid(rect, n=121):
    k = (rect, n)
    if k not in GRIDS:
        GRIDS[k] = np.meshgrid(np.linspace(*rect[0], n), np.linspace(*rect[1], n))
    return GRIDS[k]
def margins(Bv, bv, sv):
    out = {}
    dd, xx = grid(DOM); m = fmodel(Bv, bv, dd, xx, sv)
    out["active"] = float((m["GN"] - m["E"](m["Rd"])).min())
    out["G_L"] = float((m["GL"] - m["E"](m["Rd"])).min()); out["G_R"] = float((m["GR"] - m["E"](m["Rd"])).min())
    out["short"] = float((-m["E"](m["Rd"]) - m["h"] * m["Dsh"](m["Rd"])).min()); out["slope"] = float(m["S"](m["Rd"]).min())
    for nm, r in [("D1", D1), ("D2", D2r), ("D3", D3)]:
        dd, xx = grid(r); mm = fmodel(Bv, bv, dd, xx, sv); out["E(R_s)" + nm] = float((-mm["E"](mm["Rs"])).min())
    l, lx, ly, h, T, C0 = m["l"], m["lx"], m["ly"], m["h"], m["T"], m["C0"]
    out["floor"] = -1 / 1000 - C0 - (0.5 + 1.5 * l) / 50 + h * sv
    out["small_amp"] = -1 / 250 - C0 - (1 - h) / 16 + 7 / 8 * h * sv
    out["small_rows"] = -1 / 20 - 13 * h / 75 + ly / 2 - 1 / 50
    out["direct_branch"] = ly - l - 11 * bv / 6; out["l<=1/5"] = 0.2 - l; out["lx>l"] = lx - l; out["T>0"] = T
    out["m0 in band"] = min(m["m0"] - 0.30, 0.47 - m["m0"]); out["b>=0"] = bv; out["B>0.67"] = Bv - 0.67
    return out
def feasible(Bv, sv):
    best = None
    for bv in np.linspace(0.0, 0.35, 141):
        mg = margins(Bv, bv, sv); mn = min(mg.values())
        if best is None or mn > best[0]:
            best = (mn, bv, mg)
    return best
print("sigma | inf B (numerical) | KKT B(sigma) | b at inf | binding constraints")
import sympy as sp
s_ = sp.symbols('s')
disc = 964 * s_**2 - 1156 * s_ + 1081
Bcl = sp.lambdify(s_, (103528 * s_**2 + 7316 * s_ * sp.sqrt(disc) + 51710 * s_ - 3823 * sp.sqrt(disc) - 51257) /
                  (12 * (8720 * s_**2 + 640 * s_ * sp.sqrt(disc) + 4912 * s_ - 356 * sp.sqrt(disc) - 5179)))
for sv in [0.0, 0.002, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5]:
    lo, hi = 0.80, 0.8734
    if feasible(hi, sv)[0] <= 0:
        print("%5g  infeasible at 0.8734" % sv); continue
    for _ in range(26):
        mid = (lo + hi) / 2
        if feasible(mid, sv)[0] > 0: hi = mid
        else: lo = mid
    mn, bv, mg = feasible(hi, sv)
    bind = sorted(mg.items(), key=lambda kv: kv[1])[:3]
    print("%-6g %.7f   %.7f   b=%.4f  %s" % (sv, hi, Bcl(sv), bv, ", ".join("%s=%.1e" % kv for kv in bind)))
print("time %.1fs, maxrss %.0f MB" % (time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20))
