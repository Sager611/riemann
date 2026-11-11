# bprime.py -- routine rerun of the Mertens78 pipeline with the Lean-checked half-plane B' = 0.87338196661320752.
# Each Mertens78 certificate is translated by -(7/8 - B'): all lines sigma_j and sigma_1 move left by D, so the
# heights sigma_j - B' are identical to the old sigma_j - 7/8 and (C1), (C2) are unchanged (DERIVATION.md section 4).
# Only kappa = (1 + eta - B')/2, A = log(3 zeta(1+eta)) + kappa log(1/(1+B')), C_inf (zeta at the new lines) and J0 change.
# Run: uv run --no-project --with python-flint python3 bprime.py
import json, math, os, sys, time
from decimal import Decimal
from flint import arb, ctx
M78 = '/Users/adriansager/playground/riemann/Riemann/Deep/Mertens78'
sys.path.insert(0, M78)
import cert as C78
ctx.prec = 120
BPs = "0.87338196661320752"; BP = arb(BPs); DD = Decimal("0.875") - Decimal(BPs)
C78.B = BP                                    # heights in Cert (Harnack mass, x1, F) now relative to B'
def sh(v): return str(Decimal(v) - DD)
def translate(cf):
    if cf.get('harnack'):
        d = dict(cf); d['sigma1'] = sh(cf['sigma1']); d['sigma0'] = sh(cf['sigma0'])
    else:
        d = json.load(open(os.path.join(M78, cf['cert']))); d['eta'] = cf['eta']; d['T1'] = '1000'
        d['sigma1'] = sh(d['sigma1']); d['masses'] = [[sh(a), b, c] for a, b, c in d['masses']]
    d['name'] = cf['name'] + "_Bp"
    return d
def make(d):
    c = C78.Cert(d)
    c.kappa = (1 + c.eta - BP)/2
    c.A = (3*(1 + c.eta).zeta()).log() + c.kappa*(1/(1 + BP)).log()
    c.a = (c.Lam - 1)*c.kappa
    c.x1 = c.s1 - BP
    return c
def run(d, J0cache):
    c = make(d); a = c.a; s1 = c.s1
    assert a > 0 and a < 1
    # contour shift (M78 Lemma B2) with B': Harnack at sigma_0 = 1.02, eta = 0.05 gives a' < 2 for sigma_1 >= 0.888
    ap = ((arb("1.02") - BP)/(s1 - BP) - 1)*(1 + arb("0.05") - BP)/2
    assert ap < 2
    th_star = (s1 + a)/(1 + a)
    theta = math.ceil(float(th_star.upper())*1e4)/1e4
    if theta <= float(th_star.upper()): theta += 1e-4
    theta_s = "%.4f" % theta
    logC = float(c.tail(d["T1"]).mid()); af = float(a.mid())
    logD = logC + math.log(1/af + 1/(1 - af)) + af*math.log(2) - math.log(math.pi)
    logch = (math.log(4*af) + logD)/(1 + af)
    l10 = logch/math.log(10); e10 = math.floor(l10); chs = "%.6fe%d" % (10**(l10 - e10), e10)
    if d['sigma1'] not in J0cache: J0cache[d['sigma1']] = C78.J0_direct(d['sigma1'], d['T1'])
    J0, nJ, mx = J0cache[d['sigma1']]
    res = []
    for k in [1.0, 1.25, 1.5, 2.0, 1e4]:
        l10 = math.ceil((math.log10(k) + logch/math.log(10))/(1 - theta))
        r = C78.final(c, d["T1"], "1e%d" % l10, chs, theta_s, J0=J0)
        lC = float(r["C"].upper().log().upper())/math.log(10)
        res.append(dict(log10xa=l10, log10C=lC, Cprime_log10=max(lC, (1 - theta)*l10)))
    best = min(res, key=lambda q: q["log10C"]); allx = min(res, key=lambda q: q["Cprime_log10"])
    return dict(name=d['name'], sigma1=d['sigma1'], eta=d['eta'], Lam=float(c.Lam.mid()), a=af, theta_star=float(th_star.mid()),
                theta=theta, ch=chs, J0=float(J0.upper()), max_inv_zeta_T1=mx, logCinf=float(c.tail(d["T1"]).upper()),
                best=best, allx=allx), c
if __name__ == '__main__':
    CFG = json.load(open(os.path.join(M78, 'all.json')))
    old = {r['name']: r for r in json.load(open(os.path.join(M78, 'results.json')))}
    J0c = {}; out = []; t0 = time.time()
    for cf in CFG:
        d = translate(cf)
        r, c = run(d, J0c)
        if cf['name'] in ('H093', 'ML093g5_e10', 'ML092g8_e10'):
            r['C1C2_recheck'] = str(c.check_C())     # identical heights; rechecked for the small certificates
        o = old[cf['name']]
        r['old_theta'] = o['theta']; r['old_log10Cprime'] = o['allx']['Cprime_log10']
        out.append(r); d2 = dict(d); json.dump(d2, open(os.path.join('logs', 'cert_' + d['name'] + '.json'), 'w'))
        print("%-18s sigma1 %s  a %.5f  theta %.4f (was %.4f)  log10 C(x>=1e%d) %.3f  log10 C'(all x) %.3f (was %.3f at theta %.4f)  [%.0fs]" %
              (r['name'], r['sigma1'], r['a'], r['theta'], o['theta'], r['best']['log10xa'], r['best']['log10C'],
               r['allx']['Cprime_log10'], o['allx']['Cprime_log10'], o['theta'], time.time() - t0), flush=True)
    json.dump(out, open(os.path.join('logs', 'bprime_results.json'), 'w'), indent=1)
