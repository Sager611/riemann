
# cellcheck.py -- pbpaudit: cert.py's own cell_bound / tail_bound re-executed at a chosen precision (128 or 256 bits).
# The module source is exec'd with ctx.prec replaced, so every module constant is recomputed at that precision.
# Memory: one kernel per cell, built only for the p-grid points j/200 that this cell can request (float pre-pass, +-1 neighbours),
# which gives the same p_up choices as the full grid.
import sys, math, json, time, types
prec = int(sys.argv[1]); cells = json.load(open(sys.argv[2]))
CERT = sys.argv[4] if len(sys.argv) > 4 else 'Riemann/Deep/Ingham3/cert.py'     # optional: another copy of cert.py (e.g. HP141's theta141 patch)
src = open(CERT).read()
assert src.count('ctx.prec = 128') == 1
src = src.replace('ctx.prec = 128', 'ctx.prec = %d' % prec)
print('cert source:', CERT, ' B =', [l for l in src.split('\n') if l.startswith('B = ')][0][:60], flush=True)
mod = types.ModuleType('certP'); mod.__dict__['__name__'] = 'certP'; mod.__file__ = 'cert_copy.py'
exec(compile(src, 'cert.py@%d' % prec, 'exec'), mod.__dict__)
C = mod; arb = C.arb
t00 = time.time()
P = C.build_kln()
LOGB = float(C.B.mid())
def needed_js(k, La, Lb, tail=False):
    js = set()
    def add(p):
        j = math.ceil(200*p - 1e-9)
        for jj in (j - 1, j, j + 1):
            if 1 <= jj <= 300: js.add(jj)
    if tail:
        Dlow = Lb/k - math.log(k); add(8*(1 - 0.515)/3 + 4/Dlow); add(1 + 1/Dlow); return js
    db = k*math.log1p(math.exp(-Lb/k)); D = -math.log(db); lam0 = db*3e12; l0 = math.log(3e12)
    add(1 + 1/D); pden = l0 if lam0 >= 1 else D; add(2/pden)
    s = 0.515
    for j in range(146):
        add(8*(1 - s)/3 + (5 - 2*s)/pden); s += 0.0025
    return js
out = []
for c in cells:
    t = time.time(); k = c['k']
    tail = c.get('tail', False)
    js = needed_js(k, c['La'], c['Lb'], tail)
    pgrid = [arb(j)/200 for j in sorted(js)]
    da = k*math.log1p(math.exp(-c['La']/k))
    dmx = [d for d in ('1e-9', '0.025', '0.15') if da < float(d)*(1 - 1e-12)][0]
    K = C.Kernel(2, arb(dmx), pgrid)
    S = C.Setup(P, {2: [K]}, pgrid, arb('0.515'), arb('0.0025'))
    if tail:
        R = C.tail_bound(S, k, c['La'], 2)
    else:
        R = C.cell_bound(S, k, c['La'], c['Lb'], 2)
    out.append(dict(c, prec=prec, R=R.str(12, radius=False), Rf=float(R.mid()), np=len(pgrid), dmax=dmx))
    print('%-10s k=%d [%.6f, %.6f] prec=%d  R=%s  (%d p-points, dmax %s, %.1fs)' % (c.get('tag', ''), k, c['La'], c['Lb'], prec, R.str(10, radius=False), len(pgrid), dmx, time.time() - t), flush=True)
json.dump(out, open(sys.argv[3], 'w'), indent=1)
print('total %.1fs' % (time.time() - t00))
