
# b3lean gen.py: exact-Fraction enclosures for Theorem E1' numerics + Lean snippet generator.
from fractions import Fraction as Fr
import math, sys, json
D = 30
P10 = 10**D
def rd(q): return Fr(math.floor(q*P10), P10)
def ru(q): return Fr(math.ceil(q*P10), P10)
def lit(q):
    q = Fr(q)
    if q.denominator == 1: return f"({q.numerator} : ℝ)"
    return f"({q.numerator} / {q.denominator} : ℝ)"
def raw(q):
    q = Fr(q)
    return f"{q.numerator} / {q.denominator}" if q.denominator != 1 else f"{q.numerator}"
th0 = Fr(371, 750)
slo, shi = Fr(141421356237, 10**11), Fr(141421356238, 10**11)
xlo = Fr(math.floor(th0/shi*10**12), 10**12)
xhi = Fr(math.ceil(th0/slo*10**12), 10**12)
assert xlo*shi <= th0 and th0 <= xhi*slo
def base_sin(a): return (a - a**3/6, a - a**3/6 + a**5/100)
def chain(X, k, name):
    """returns (sl, sh, cl, ch) at X and Lean proof lines"""
    a0 = X / 2**k; b = a0/2
    L = []
    sbl, sbh = base_sin(b); sbl, sbh = rd(sbl), ru(sbh)
    L.append(f"  have bb := base_sin F (a := {raw(b)}) (by norm_num) (by norm_num)")
    L.append(f"  have sb : {lit(sbl)} ≤ F.sin ({raw(b)}) ∧ F.sin ({raw(b)}) ≤ {lit(sbh)} := ⟨le_trans (by norm_num) bb.1, le_trans bb.2 (by norm_num)⟩")
    cl, ch = rd(1 - 2*sbh**2), ru(1 - 2*sbl**2)
    L.append(f"  have c0 := cos_two F (cl' := {lit(cl)}) (ch' := {lit(ch)}) sb (by norm_num) (by norm_num) (by norm_num)")
    L.append(f"  rw [show (2 : ℝ) * ({raw(b)}) = {raw(a0)} by norm_num] at c0")
    sl, sh = base_sin(a0); sl, sh = rd(sl), ru(sh)
    L.append(f"  have ba := base_sin F (a := {raw(a0)}) (by norm_num) (by norm_num)")
    L.append(f"  have s0 : {lit(sl)} ≤ F.sin ({raw(a0)}) ∧ F.sin ({raw(a0)}) ≤ {lit(sh)} := ⟨le_trans (by norm_num) ba.1, le_trans ba.2 (by norm_num)⟩")
    a = a0; prev = "s0 c0"
    for i in range(1, k+1):
        assert sl >= 0 and cl >= 0
        sl2, sh2, cl2, ch2 = rd(2*sl*cl), ru(2*sh*ch), rd(1 - 2*sh**2), ru(1 - 2*sl**2)
        L.append(f"  have d{i} := dbl F (sl' := {lit(sl2)}) (sh' := {lit(sh2)}) (cl' := {lit(cl2)}) (ch' := {lit(ch2)}) {prev}")
        L.append(f"    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)")
        L.append(f"  rw [show (2 : ℝ) * ({raw(a)}) = {raw(2*a)} by norm_num] at d{i}")
        a = 2*a; sl, sh, cl, ch = sl2, sh2, cl2, ch2; prev = f"d{i}.1 d{i}.2"
    assert a == X
    return (sl, sh, cl, ch), L, f"d{k}"
# reference values (Arb, b3c_arb.log)
c_arb = -0.1857718401393168364675068
k_arb = 0.0850387603917657602445070238
(slA, shA, clA, chA), LA, dA = chain(xlo, 5, "xlo")
(slB, shB, clB, chB), LB, dB = chain(xhi, 5, "xhi")
cot_hi = ru(chA/slA); cot_lo = rd(clB/shB)
cw_lo = rd(2 - th0/2 - cot_hi/slo); cw_hi = ru(2 - th0/2 - cot_lo/shi)
print("xlo", xlo, "xhi", xhi)
print("cW enclosure", float(cw_lo), float(cw_hi), "width", float(cw_hi-cw_lo), "vs arb", float(cw_lo)-c_arb, float(cw_hi)-c_arb)
json.dump({"xlo": str(xlo), "xhi": str(xhi), "cot_hi": str(cot_hi), "cot_lo": str(cot_lo), "cw_lo": str(cw_lo), "cw_hi": str(cw_hi),
  "sinA": str(slA), "cosA": str(chA), "cosB": str(clB), "sinB": str(shB)}, open("Riemann/Deep/B3Cert/Lean/audit/trig.json", "w"))
open("Riemann/Deep/B3Cert/Lean/audit/chainA.txt", "w").write("\n".join(LA) + "\n")
open("Riemann/Deep/B3Cert/Lean/audit/chainB.txt", "w").write("\n".join(LB) + "\n")
print("chain lines", len(LA), len(LB))

