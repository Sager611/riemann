#!/usr/bin/env python3
# scan.py -- float model scan over s* and row sets (not a certificate)
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import model as M
B = 0.87334585610671857680; aB = B - 0.5
R4 = ('0.85', '0.86', '0.87', '0.88'); RALL = tuple(r['name'] for r in M.rows['main'])
for use, lab in ((R4, 'rows4'), (RALL, 'rowsAll')):
    for s0 in (0.80, 0.83, 0.845, 0.85, 0.855, 0.86):
        sig = list(np.linspace(s0, B, int(round((B - s0) / 0.00015)) + 1))
        for mt, lf in ((False, False), (True, True)):
            c = M.cmax(lambda n: M.Qall(n, s0 - 0.5) + M.Qdens(n, sig, aB, mt, lf, use), it=30)
            print('%s s*=%.3f J=%d mintot=%d logform=%d: c = %.4f' % (lab, s0, len(sig) - 1, mt, lf, c), flush=True)
