# Build a certified re-optimised table: q0 = (2 theta - 1)^2 rounded up, grid t_i = i*h, last step chosen minimal.
import sys, time, json, math
from fractions import Fraction as F
from dbn import psi_matrix, bernstein_coeffs, lp_best, kexact, E, CG_exact, P_, B_
from arbsc import arb_scalars, ub_int, KEYS, DEN
t00 = time.time()
TH = F(5458637291332547, 6250000000000000)
QDEN = 10 ** 10; CDEN = 10 ** 8; LDEN = 10 ** 6; DDEN = 10 ** 7
b = F(1, 25)
h = F(1, int(sys.argv[1])) if len(sys.argv) > 1 else F(1, 4000)
mode = sys.argv[2] if len(sys.argv) > 2 else 'new'
q0exact = (2 * TH - 1) ** 2 if mode == 'new' else F(9, 16)
out = sys.argv[3] if len(sys.argv) > 3 else 'rows.jsonl'
MU, MV = 16, 1
def ceilq(x, den): return F(-((-x.numerator * den) // x.denominator), den)
def floorq(x, den): return F((x.numerator * den) // x.denominator, den)

def row_data(tm, tp, q):
    s = arb_scalars(tm, tp)
    U = {k: ub_int(s[k]) for k in KEYS}
    Cc, Gc, ep0 = CG_exact({k: F(U[k], DEN) for k in KEYS})
    C = floorq(Cc, CDEN); G = ceilq(Gc, CDEN)
    lamf, muf, kf = lp_best(float(q), float(C), float(G), MU, MV)
    mu = F(round(muf * LDEN), LDEN)
    b0 = bernstein_coeffs(psi_matrix(q, F(0), F(0)), MU, MV, F(1))
    b1 = bernstein_coeffs(psi_matrix(q, F(1), F(0)), MU, MV, F(1))
    b2 = bernstein_coeffs(psi_matrix(q, F(0), F(1)), MU, MV, F(1))
    hi = None
    for c0, x1, x2 in zip(b0, b1, b2):
        cl = x1 - c0; base = c0 + mu * (x2 - c0)
        if cl < 0:
            lim = base / (-cl); hi = lim if hi is None or lim < hi else hi
    lam = floorq(hi, LDEN)
    if lam == hi: lam -= F(1, LDEN)
    beta = [c0 + lam * (x1 - c0) + mu * (x2 - c0) for c0, x1, x2 in zip(b0, b1, b2)]
    assert min(beta) > 0, "Bernstein fail"
    k = kexact(lam, mu, C, G, q)
    assert k > 0
    return dict(U=U, C=C, G=G, lam=lam, mu=mu, k=k, eps0=ep0, minbeta=min(beta), kf=kf)

rows = []
q = ceilq(q0exact, QDEN); t = F(0); i = 0
fo = open(out, 'w')
def dump(r):
    fo.write(json.dumps({kk: (str(v) if isinstance(v, F) else v) for kk, v in r.items()}) + "\n"); fo.flush()
while True:
    rd = row_data(t, t + h, q)
    k = rd['k']
    qn = ceilq(((2 + k * q) / E(k * h, 20) - 2) / k, QDEN)
    if qn > b:
        r = dict(i=i, t=t, d=h, q=q, qnext=qn, **{kk: rd[kk] for kk in ('C', 'G', 'lam', 'mu')}, U=rd['U'], k=float(k), eps0=float(rd['eps0']), minbeta=float(rd['minbeta']))
        rows.append(r); dump(r)
        if i % 25 == 0: print("row %d t=%s q=%.10f k=%.6f eps0=%.5f minbeta=%.3e  %.1fs" % (i, t, float(q), float(k), float(rd['eps0']), float(rd['minbeta']), time.time() - t00), flush=True)
        t += h; q = qn; i += 1
        continue
    # final step: minimal d (multiple of 1/DDEN) with exact gate to b
    df = math.log((2 + float(k) * float(q)) / (2 + float(k) * float(b))) / float(k)
    d = ceilq(F(df), DDEN)
    while True:
        rd = row_data(t, t + d, q); k = rd['k']
        if E(k * d, 20) * (2 + k * b) >= 2 + k * q: break
        d += F(1, DDEN)
    # try to shorten
    while d > F(1, DDEN):
        d2 = d - F(1, DDEN); rd2 = row_data(t, t + d2, q); k2 = rd2['k']
        if E(k2 * d2, 20) * (2 + k2 * b) >= 2 + k2 * q: d, rd, k = d2, rd2, k2
        else: break
    r = dict(i=i, t=t, d=d, q=q, qnext=b, **{kk: rd[kk] for kk in ('C', 'G', 'lam', 'mu')}, U=rd['U'], k=float(k), eps0=float(rd['eps0']), minbeta=float(rd['minbeta']))
    rows.append(r); dump(r)
    tn = t + d
    print("FINAL: rows=%d t_n=%s=%.8f  Lambda <= %s = %.8f   (%.1fs)" % (len(rows), tn, float(tn), tn + b / 2, float(tn + b / 2), time.time() - t00), flush=True)
    break
fo.close()
