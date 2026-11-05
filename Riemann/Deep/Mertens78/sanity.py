# sanity.py -- numerical sanity checks (not part of the proof): Phi(t) >= |1/zeta|, Rademacher transcription
import json
from flint import arb, acb, ctx
from cert import Cert
ctx.prec = 80
def rad(s, eta):
    s = acb(s); eta = arb(eta)
    return 3*abs(1+s)/abs(1-s)*(abs(1+s)/(2*arb.pi()))**((1+eta-s.real)/2)*(1+eta).zeta()
worst = 1e9
for sig in ["0.876", "0.9", "0.95", "1.0", "1.05", "1.2"]:
    for t in [0.5, 3, 14.1347, 30, 100, 1000, 12345.6, 1e5]:
        for eta in ["0.005", "0.05", "0.5"]:
            s = acb(arb(sig), arb(t))
            z = abs(s.zeta()); R = rad(s, eta)
            worst = min(worst, float((R/z).lower()))
print("Rademacher: min ratio bound/|zeta| over grid =", worst)
for name, d in [("H093", {"sigma1":"0.93","sigma0":"1.02","eta":"0.1","harnack":True}),
                ("ML093g5", dict(json.load(open("c_ml093_g5.json")), eta="0.1")),
                ("S089g3000", dict(json.load(open("c_s089_g3000.json")), eta="0.005"))]:
    c = Cert(d); mn = 1e300; rows = []
    for t in [0.0, 0.3, 1, 5, 14.1347, 50, 200, 1000, 5000, 2e4, 1e5]:
        lp = c.logPhi(arb(t))
        z = abs(acb(c.s1, arb(t)).zeta())
        inv = 1/z
        ratio = float((lp - inv.log()).lower())
        mn = min(mn, ratio); rows.append((t, round(float(lp.mid()), 2), round(float(inv.mid()), 4)))
    print(name, "min log(Phi*|zeta|) =", round(mn, 3), rows[:4], rows[-3:])

