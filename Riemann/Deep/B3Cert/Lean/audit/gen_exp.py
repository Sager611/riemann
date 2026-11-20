
# b3lean gen_exp.py: exp chains for kappa_2 (exact Fractions) + Lean snippets.
from fractions import Fraction as Fr
import math, json
D = 30
P10 = 10**D
def rd(q): return Fr(math.floor(q*P10), P10)
def ru(q): return Fr(math.ceil(q*P10), P10)
def lit(q):
    q = Fr(q)
    return f"({q.numerator} / {q.denominator} : ℝ)" if q.denominator != 1 else f"({q.numerator} : ℝ)"
def raw(q):
    q = Fr(q)
    return f"{q.numerator} / {q.denominator}" if q.denominator != 1 else f"{q.numerator}"
clo = Fr(3138997965019943478277, 10**18); chi = Fr(3138997965019943478278, 10**18)
elo = Fr(27182818283, 10**10); ehi = Fr(27182818286, 10**10)
e8lo, e8hi = elo**8, ehi**8
logc = math.log(3138.99796501994347827760614809)  # float reference only
def lower_chain(f, m):
    y = f / 2**m; L = []
    l = rd(1 + y + y**2/2)
    L.append(f"  have q0 := F.quadratic_le_exp_of_nonneg ({raw(y)}) (by norm_num)")
    L.append(f"  have l0 : {lit(l)} ≤ F.exp ({raw(y)}) := le_trans (by norm_num) q0")
    prev = "l0"
    for i in range(1, m+1):
        l2 = rd(l*l)
        L.append(f"  have l{i} := sq_lo F (l' := {lit(l2)}) {prev} (by norm_num) (by norm_num)")
        L.append(f"  rw [show (2 : ℝ) * ({raw(y)}) = {raw(2*y)} by norm_num] at l{i}")
        y = 2*y; l = l2; prev = f"l{i}"
    assert y == f
    return l, L, prev
def upper_chain(g, m):
    y = g / 2**m; L = []
    u = ru(1/(1 - y))
    L.append(f"  have b0 := F.exp_bound_div_one_sub_of_interval ({raw(y)}) (by norm_num) (by norm_num)")
    L.append(f"  have u0 : F.exp ({raw(y)}) ≤ {lit(u)} := le_trans b0 (by norm_num)")
    prev = "u0"
    for i in range(1, m+1):
        u2 = ru(u*u)
        L.append(f"  have u{i} := sq_up F (u' := {lit(u2)}) {prev} (by norm_num)")
        L.append(f"  rw [show (2 : ℝ) * ({raw(y)}) = {raw(2*y)} by norm_num] at u{i}")
        y = 2*y; u = u2; prev = f"u{i}"
    assert y == g
    return u, L, prev
# L_hi: smallest 10-decimal L with chi <= e8lo * lower(L-8)
mL, mU = 8, 18
L0 = Fr(math.floor(logc*10**10), 10**10)
for j in range(0, 200):
    Lhi = L0 + Fr(j, 10**10)
    l, LA, pA = lower_chain(Lhi - 8, mL)
    if chi <= e8lo * l: break
for j in range(0, 200):
    Llo = L0 + Fr(1, 10**10) - Fr(j, 10**10)
    u, LB, pB = upper_chain(Llo - 8, mU)
    if e8hi * u <= clo: break
R = Fr(44, 5)
klo = 1 - Lhi / R; khi = 1 - Llo / R
print("Lhi", Lhi, float(Lhi) - logc, "Llo", Llo, float(Llo) - logc)
print("kappa enclosure", float(klo) - 0.0850387603917657602445, float(khi) - 0.0850387603917657602445)
json.dump({"Lhi": str(Lhi), "Llo": str(Llo), "klo": str(klo), "khi": str(khi), "lA": str(l), "uB": str(u), "mL": mL, "mU": mU}, open("Riemann/Deep/B3Cert/Lean/audit/exp.json", "w"))
open("Riemann/Deep/B3Cert/Lean/audit/chainL.txt", "w").write("\n".join(LA) + "\n")
open("Riemann/Deep/B3Cert/Lean/audit/chainU.txt", "w").write("\n".join(LB) + "\n")
print("lines", len(LA), len(LB), pA, pB)

