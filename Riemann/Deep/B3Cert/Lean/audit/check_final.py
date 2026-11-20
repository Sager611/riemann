
from fractions import Fraction as Fr
import json
from decimal import Decimal, getcontext
getcontext().prec = 60
t = json.load(open("Riemann/Deep/B3Cert/Lean/audit/trig.json")); e = json.load(open("Riemann/Deep/B3Cert/Lean/audit/exp.json"))
cwlo, cwhi = Fr(t["cw_lo"]), Fr(t["cw_hi"]); klo, khi = Fr(e["klo"]), Fr(e["khi"])
CWL, CWH = Fr(-1857718405, 10**10), Fr(-18577184, 10**8)
KL, KH = Fr(850387602, 10**10), Fr(850387611, 10**10)
assert CWL <= cwlo and cwhi <= CWH and KL <= klo and khi <= KH
def sq(q): return Decimal(q.numerator)/Decimal(q.denominator)
def hc(c, k): return (3 + sq(k) - (sq((1-k)*(9-k-8*c))).sqrt())/4
print("h(CWL,KL) =", hc(CWL, KL), " h(CWH,KH) =", hc(CWH, KH))
for HL in [Fr(335, 10**7), Fr(33545, 10**9)]:
    B = 3 + KL - 4*HL; print("HL", float(HL), B >= 0 and (1-KL)*(9-KL-8*CWL) < B*B)
for HH in [Fr(33547, 10**9), Fr(33546, 10**9)]:
    B = 3 + KH - 4*HH; print("HH", float(HH), B*B <= (1-KH)*(9-KH-8*CWH))
ks_hi = 1 - Fr(2)/(2 - CWL); ks_lo = 1 - Fr(2)/(2 - CWH)
print("kappa-k* lower", float(KL - ks_hi), " upper", float(KH - ks_lo))
for DL in [Fr(47345, 10**9), Fr(473455, 10**10)]: print("DL", float(DL), KL - ks_hi > DL)
for DH in [Fr(47347, 10**9), Fr(47348, 10**9)]: print("DH", float(DH), KH - ks_lo <= DH)
print("kstar enclosure", float(ks_lo), float(ks_hi))

