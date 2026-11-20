
# b3lean check_det.py: E1Lean's xs₁ / r_P (Detector.lean) vs params_exp010.json set nu = 999/10000 (the set b3c_arb.py reads).
import json, re
from fractions import Fraction as Fr
src = open("Riemann/Deep/SantiAudit/Lean/Detector.lean").read()
blk = src[src.index("def xs₁"):]; blk = blk[:blk.index("]")]
lean_xs = [Fr(int(p), int(q)) for p, q in re.findall(r"\((-?\d+) / (\d+) : ℝ\)", blk)]
rP = re.search(r"Pfam \((\d+) / (\d+)\)", src); lean_r = Fr(int(rP.group(1)), int(rP.group(2)))
D = json.load(open("Riemann/Deep/SantiAudit/params_exp010.json"))
sets = [p for p in D["parameter_sets"] if p["nu"] == "999/10000"]
ps = sets[0]
q = lambda s: Fr(*map(int, s.split("/"))) if "/" in str(s) else Fr(int(s))
js = [q(x) for x in ps["x_j"]]
print("sets with nu=999/10000:", len(sets), "| lean xs:", len(lean_xs), "json x_j:", len(js), "| equal in order:", lean_xs == js)
print("r_P lean:", lean_r, "json:", q(ps["r_P"]), "equal:", lean_r == q(ps["r_P"]))
arb = open("Riemann/Deep/B3Cert/code/b3c_arb.py").read()
print("b3c_arb selects nu == '999/10000':", "p['nu'] == '999/10000'" in arb, "| uses r_P:", "ps['r_P']" in arb, "| index start=1, T[2*j+1]:", "enumerate(xs, start=1)" in arb and "T[2*j+1] - u" in arb)
print("all parameter_sets nu:", [p["nu"] for p in D["parameter_sets"]])

