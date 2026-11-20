"""gapjinf: independent validity check of a saved form (logs/form_*.json) on explore53's polynomial basis (validate53.py,
copied): min over 600 gaps g in (0, pi] of the smallest eigenvalue of the form relative to the H^J Gram matrix.
Also h of the form on its moments.  usage: valform.py form.json [N]"""
import sys, json
import numpy as np
from validate53 import minform
F = json.load(open(sys.argv[1])); N = int(sys.argv[2]) if len(sys.argv) > 2 else 30
J = F['J']; a = F['a']; phi = {tuple(int(t) for t in k.split(',')): v for k, v in F['phi'].items()}
m = F['moments']
h = -a[0] * m[0] / sum(a[j] * m[j] for j in range(1, J + 1))
ev, g = minform(a, phi, J, N=N)
print('%s: J=%d h=%.7f (LP %.7f); polynomial-basis check N=%d: min rel eig %.3e at g/pi=%.4f' % (sys.argv[1], J, h, F['h'], N, ev, g / np.pi))
