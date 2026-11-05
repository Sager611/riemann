# run.py -- rigorous end-to-end bound for one parameter set
import json, sys, math
from flint import arb, ctx
from cert import Cert, J0_direct, final, B
ctx.prec = 120
def run(d, J0cache={}):
    c = Cert(d)
    a = c.a; s1 = c.s1
    th_star = (s1 + a)/(1 + a)
    theta = math.ceil(float(th_star.upper())*1e4)/1e4
    if theta <= float(th_star.upper()): theta += 1e-4
    theta_s = "%.4f" % theta
    logC = float(c.tail(d["T1"]).mid()); af = float(a.mid())
    logD = logC + math.log(1/af + 1/(1-af)) + af*math.log(2) - math.log(math.pi)
    logch = (math.log(4*af) + logD)/(1+af)
    l10 = logch/math.log(10); e10 = math.floor(l10); chs = "%.6fe%d" % (10**(l10 - e10), e10)
    key = (d["sigma1"], d["T1"])
    if key not in J0cache: J0cache[key] = J0_direct(d["sigma1"], d["T1"])
    J0, nJ, mx = J0cache[key]
    res = []
    for k in d.get("ks", [1.0, 1.5, 2.0, 1e4]):
        l10 = math.ceil((math.log10(k) + logch/math.log(10))/(1 - theta))
        r = final(c, d["T1"], "1e%d" % l10, chs, theta_s, J0=J0)
        lC = float(r["C"].upper().log().upper())/math.log(10)
        res.append(dict(k=k, log10xa=l10, C=10**min(lC, 300), log10C=lC, Cprime_log10=max(lC, (1-theta)*l10)))
    out = dict(sigma1=d["sigma1"], eta=d["eta"], cfg=d.get("sigma0", "LP"), Lam=float(c.Lam.mid()), a=af,
               theta_star=float(th_star.mid()), theta=theta, ch=chs, J0=float(J0.upper()), max_inv_zeta_T1=mx,
               logCinf=float(c.tail(d["T1"]).upper()), runs=res)
    best = min(res, key=lambda q: q["log10C"])            # asymptotic-ish constant (x >= xa)
    allx = min(res, key=lambda q: q["Cprime_log10"])      # valid for all x >= 1
    out["best"] = best; out["allx"] = allx
    def cross(log10C):
        Cl = log10C*math.log(10); L = 10.0
        while Cl + theta*L > math.log(0.013) + L - math.log(L): L += 0.01
        return L/math.log(10)
    out["ramare_log10_allx"] = cross(allx["Cprime_log10"])
    out["X"] = {str(t): (max(best["log10C"]/(t - theta), best["log10xa"]) if t > theta else None) for t in [0.95, 0.94, 0.93, 0.92, 0.91, 0.90]}
    return out, c
if __name__ == "__main__":
    d = json.loads(sys.argv[1])
    out, c = run(d)
    print("C1/C2:", c.check_C())
    print(json.dumps(out, indent=1))

