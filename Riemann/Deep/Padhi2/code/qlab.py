"""padhi2 step B (full scan): Padhi's certified q_arb on every published interval with his kernels, with a labelled copy
of scan_arb.logN that records which bound attains the minimum in N+ (N1 zero count, N2 KLN, N3u Bellotti uniform,
N3e29 Bellotti T<=e^29, N3r Bellotti Table 2 rows, N4 Lemma D cells). The copy is checked against scan_arb.logN.
If no N3 label ever wins, removing N3 (or raising its C, e.g. Bellotti's corrected Table 2) leaves every q unchanged.
usage: qlab.py <variant> [L1min L1max]"""
import sys, os, glob, time, math, bisect, json, collections
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(__file__))
import common as CM
from flint import arb
S = CM.S
var = sys.argv[1]; lo = float(sys.argv[2]) if len(sys.argv) > 2 else 0; hi = float(sys.argv[3]) if len(sys.argv) > 3 else 1e9
CM.set_variant(var)
ORIG_LOGN = S.logN
COUNT = collections.Counter(); N3WIN = []; CHECK = [0]
def logN_lab(s, lT):
    sF = Fr(s); P = S._lt(lT); T, lTa = P['T'], P['lTa']
    cands = [(P['all'], 'N1')]
    i = bisect.bisect_right(S.KLN_S, sF) - 1
    if i >= 0 and T > S.H0K:
        if i not in P['kln']:
            s0, A, B = S.KLN[i]; s0a = arb(s0.numerator)/s0.denominator
            P['kln'][i] = (A*lTa**(2*s0a)*lTa**(5 - 4*s0a)*T**(arb(8)/3*(1 - s0a)) + B*lTa**2).log()
        cands.append((P['kln'][i], 'N2'))
    one_s = arb(1) - arb(s)
    if sF >= Fr('0.9927') and lT > math.log(3e12):
        cands.append((S.LOG_C1 + arb('1.448')*one_s*lTa, 'N3u'))
        if lT <= 29: cands.append((S.LOG_C2 + arb('14.160')*one_s*lTa, 'N3e29'))
    for j, (a, b, al, C, B) in enumerate(S.ROWS):
        lo_ = math.log(a) if j == 0 else a
        if lo_ < lT <= b and sF >= al:
            cands.append((S.ROWLOG[j] + B*one_s*lTa, 'N3r'))
    de = 1 - sF
    if de <= Fr(3, 100) and 28 < lT <= 40:
        k = math.ceil(de*10**4)
        lTF = Fr(lT)
        lk = Fr(math.ceil(2*lTF), 2) if lT <= 31 else Fr(math.ceil(lTF))
        key = '%d,%s' % (k, (int(lk) if lk.denominator == 1 else float(lk)))
        if key in S.CELL: cands.append((S.CELL[key], 'N4'))
    best = min(cands, key=lambda c: float(c[0].upper()))
    bu = float(best[0].upper())
    labs = sorted(set(l for v, l in cands if float(v.upper()) == bu))
    if len(labs) > 1: COUNT['tie:' + '+'.join(labs)] += 1
    COUNT[best[1]] += 1
    if best[1].startswith('N3') and len(N3WIN) < 50: N3WIN.append((s, lT, labs))
    if CHECK[0] < 3000:                                   # the copy must agree with Padhi's logN
        CHECK[0] += 1; ref = ORIG_LOGN(s, lT)
        assert float(ref.upper()) == bu, (s, lT, ref, best)
    return best[0]
S.logN = logN_lab
ivs = []
for f in sorted(glob.glob(os.path.join(S.HERE, 'scan_arb_*.txt'))):
    for line in open(f):
        p = line.split(); ivs.append((float(p[0]), float(p[1]), float(p[2]), int(p[3][2:]), float(p[4][2:])))
stride = int(sys.argv[4]) if len(sys.argv) > 4 else 1
ivs.sort(); sel = [v for v in ivs if lo <= v[0] < hi][::stride]
t0 = time.time(); mx = (0, None); nd = 0; maxd = 0.0
for L1, L2, qpub, cc, fr in sel:
    before = dict(COUNT)
    qa, parts = S.q_arb(L1, L2, 3.0e12, cc, fr)
    qn = S.fup_safe(qa); d = qn - qpub; nd += (abs(d) > 1e-12); maxd = max(maxd, abs(d))
    if qn > mx[0]: mx = (qn, (L1, L2))
    dc = {k: COUNT[k] - before.get(k, 0) for k in COUNT if COUNT[k] - before.get(k, 0)}
    print('%g %g c=%d f=%.2f q_pub=%.9f q=%.9f diff=%.1e %s %.0fs' % (L1, L2, cc, fr, qpub, qn, d, json.dumps(dc), time.time() - t0), flush=True)
print('variant %s: %d intervals [%g, %g), max q %.10f on %s; %d differ from published (max |diff| %.2e); %.0fs'
      % (var, len(sel), sel[0][0], sel[-1][1], mx[0], mx[1], nd, maxd, time.time() - t0))
print('argmin counts:', json.dumps(dict(COUNT)))
print('N3 wins (first 50):', N3WIN)
