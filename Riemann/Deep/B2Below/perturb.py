# Perturbations of the ceiling: kept margins (freedoms) and stronger lemma constants (sensitivities, not freedoms).
import sys
from ceiling import curve, mk_curve, best
Z = dict(LEFT=0.70, STEP=1 / 400000, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0, gate1=0)
def top(A, Bc): return (A - 0.25) / (Bc - 0.5) - 1e-9     # start where thr = 1/2
cases = [("ceiling: continuum, zero pads", dict(Z)),
         ("gate-1 margin 3e-5 kept", dict(Z, gate1=3e-5)),
         ("k pad 1/5000 kept", dict(Z, padk=1 / 5000)),
         ("STEP 1/4000, zero pads", dict(Z, STEP=1 / 4000)),
         ("TOP 0.7750 (later start)", dict(Z, TOP=0.775)),
         ("[lemma] Poisson radius 1-5u (u=(1-H)/5)", dict(Z, udiv=5)),
         ("[lemma] F0 constant 5/6 -> 5/6+0.001", dict(Z, A=5 / 6 + 0.001, TOP=top(5 / 6 + 0.001, 1.25))),
         ("[lemma] F0 constant 5/4 -> 5/4-0.001", dict(Z, Bc=1.25 - 0.001, TOP=top(5 / 6, 1.25 - 0.001)))]
for nm, kw in cases:
    Gf, lo, hi = curve(**kw)
    B, bo, h, d, x = best(mk_curve(Gf, lo, hi, (0, 1, 0, 1)))
    print("%-42s G(0.7524)=%.7f  B*=%.10f  (b=%.5f, contact %.4f)" % (nm, float(Gf(0.7524)), B, bo, d)); sys.stdout.flush()
