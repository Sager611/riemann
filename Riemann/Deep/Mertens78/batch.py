# batch.py -- certify a list of configurations; appends to results.json
import json, sys, os, time
from flint import ctx
from run import run
cfgs = json.loads(open(sys.argv[1]).read())
res = json.load(open("results.json")) if os.path.exists("results.json") else []
for cf in cfgs:
    t0 = time.time()
    if cf.get("harnack"):
        d = dict(cf)
    else:
        d = json.load(open(cf["cert"])); d["eta"] = cf["eta"]; d["T1"] = "1000"
        if "ks" in cf: d["ks"] = cf["ks"]
    out, c = run(d)
    chk = c.check_C() if not cf.get("skipcheck") else "checked separately by chk.py (see chk logs)"
    assert c.a > 0 and c.a < 1   # hypotheses of Lemma B5 / Theorem B6
    out["C1C2"] = chk; out["name"] = cf["name"]; out["secs"] = round(time.time() - t0, 1)
    out["cert"] = cf.get("cert")
    print(cf["name"], "theta", out["theta"], "log10C(x>=xa)", round(out["best"]["log10C"], 3), "xa 1e%d" % out["best"]["log10xa"],
          "C'all", round(out["allx"]["Cprime_log10"], 3), "ramare", round(out["ramare_log10_allx"], 1),
          "X", {k: (round(v, 1) if v else None) for k, v in out["X"].items()}, chk, out["secs"], "s", flush=True)
    res = [r for r in res if r.get("name") != cf["name"]] + [out]
    json.dump(res, open("results.json", "w"), indent=1)

