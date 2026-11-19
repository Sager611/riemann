"""padhi2 step B: recompute Padhi's certified q (scan_arb.q_arb, his kernels (c, f) as published in scan_arb_{0,1}.txt)
under N3 variants: published | both | all | none (see common.set_variant). Any admissible kernel gives a valid bound, so
keeping his kernels is legitimate. usage: qvar.py <variant> <L1min> <L1max> [stride]; prints L1 L2 q_published q_new."""
import sys, os, glob, time, math
sys.path.insert(0, os.path.dirname(__file__))
import common as CM
S = CM.S
var, lo, hi = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]); stride = int(sys.argv[4]) if len(sys.argv) > 4 else 1
CM.set_variant(var)
ivs = []
for f in sorted(glob.glob(os.path.join(S.HERE, 'scan_arb_*.txt'))):
    for line in open(f):
        p = line.split(); ivs.append((float(p[0]), float(p[1]), float(p[2]), int(p[3][2:]), float(p[4][2:])))
ivs.sort(); sel = [v for v in ivs if lo <= v[0] < hi][::stride]
t0 = time.time(); mx = (0, None); nd = 0
for L1, L2, qpub, cc, fr in sel:
    qa, parts = S.q_arb(L1, L2, S.H if hasattr(S, 'H') else 3.0e12, cc, fr)
    qn = S.fup_safe(qa); d = qn - qpub; nd += (abs(d) > 1e-12)
    if qn > mx[0]: mx = (qn, (L1, L2))
    print('%g %g c=%d f=%.2f q_pub=%.9f q_%s=%.9f diff=%.2e  %.0fs' % (L1, L2, cc, fr, qpub, var, qn, d, time.time() - t0), flush=True)
print('variant %s: %d intervals, max q %.9f on %s, %d differ from published by > 1e-12, %.0fs' % (var, len(sel), mx[0], mx[1], nd, time.time() - t0))
