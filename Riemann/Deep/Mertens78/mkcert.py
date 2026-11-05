# mkcert.py -- produce a rounded LP certificate JSON (masses as exact decimal strings)
import json, sys
from lp import lp
cfg = json.loads(sys.argv[1])
masses, st, dt = lp(cfg["s1"], cfg["lines"], cfg["G"], margin=cfg.get("margin", 2e-3), mabs=cfg.get("mabs", 1e-3), Ymax=cfg.get("Ymax", 3.0), dy=cfg.get("dy", 0.04), nfine=cfg.get("nfine", 601))
assert st == 0
out = []
for s0, y, l in masses:
    if abs(l) < 1e-8: continue
    out.append(["%.6g" % s0, "%.6g" % y, "%.9g" % l])
d = {"sigma1": "%.6g" % cfg["s1"], "eta": cfg["eta"], "masses": out, "T1": "1000", "G_lp": cfg["G"], "lines": cfg["lines"]}
json.dump(d, open(cfg["out"], "w"), indent=0)
print("masses", len(out), "time", round(dt, 1), "Lam", sum(float(m[2]) for m in out))

