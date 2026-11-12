# t4_bisect.py (explore39): flip-certified ranges c_flip per parity (Galerkin crossing = upper end, Barta-LP crossing = lower end).
import sys, numpy as np, math, json, time
sys.path.insert(0, '.scratch/explore39'); from xform import Form, barta_lp
kind = sys.argv[1]; sec = sys.argv[2]; lo, hi = float(sys.argv[3]), float(sys.argv[4]); n = 20
def gal(c):
    F = Form(c, kind); Qm, p = F.Q(n, sec); Fm, _ = F.flip(n, sec); e, V = np.linalg.eigh(Qm - Fm); return float(e[0]), F, V[:, 0]
def bar(c):
    g, F, v = gal(c); b, _, _ = barta_lp(F, n, sec, flip=True, v0=v); return b
for name, fn in (('galerkin', lambda c: gal(c)[0]), ('barta', bar)):
    a, b = lo, hi
    fa, fb = fn(a), fn(b)
    if fa is None or fb is None or fa*fb > 0:
        print(json.dumps({'kind': kind, 'sec': sec, 'which': name, 'no_bracket': [a, fa, b, fb]})); continue
    for _ in range(14):
        m = (a + b)/2; fm = fn(m)
        if fm is None: fm = -1.0
        if fm > 0: a = m
        else: b = m
    print(json.dumps({'kind': kind, 'sec': sec, 'which': name, 'c_flip': [round(a, 5), round(b, 5)]}), flush=True)
