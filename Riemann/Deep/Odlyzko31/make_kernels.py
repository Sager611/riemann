#!/usr/bin/env python3
"""odlyzko31: choose kernels (float optimisation only; the certificate is certify.py).
For each strip half-width w and target (n, r1): maximise [r1 logA + (n-r1) logB - E]/n over
F = S(x/b) sech(w x), S in {G, H}, b > 0 (wp = w is optimal in every tested case, see STATE/logs).
Writes kernels.json."""
import math, json, sys
import numpy as np
from kern import consts

W = {'classical': ('1/2', 0.5), 'ct': ('15/32', 15 / 32), 'ab': ('3733819666132076/10000000000000000', 0.3733819666132076)}
NS = list(range(2, 31)) + list(range(32, 52, 2)) + [52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100] \
    + list(range(110, 210, 10)) + list(range(220, 420, 20)) + [480, 500, 600, 660, 720, 840, 960, 1000, 1200, 1332, 1500, 2000,
    2400, 3000, 4000, 4800, 5000, 6000, 7000, 8000, 9000, 10000, 100000, 1000000]
MIXED_MAX = 12

def val(n, r1, s, b, w):
    lA, lB, E = consts(s, b, w)
    return (r1 * lA + (n - r1) * lB - E) / n

def best(n, r1, w):
    bs = np.exp(np.linspace(math.log(0.1), math.log(300.0), 160))
    cand = []
    for s in ('G', 'H'):
        vs = [val(n, r1, s, b, w) for b in bs]
        i = int(np.argmax(vs))
        lo, hi = math.log(bs[max(i - 1, 0)]), math.log(bs[min(i + 1, len(bs) - 1)])
        g = (math.sqrt(5) - 1) / 2
        a, c = lo, hi
        x1, x2 = c - g * (c - a), a + g * (c - a)
        f1, f2 = val(n, r1, s, math.exp(x1), w), val(n, r1, s, math.exp(x2), w)
        for _ in range(60):
            if f1 > f2:
                c, x2, f2 = x2, x1, f1; x1 = c - g * (c - a); f1 = val(n, r1, s, math.exp(x1), w)
            else:
                a, x1, f1 = x1, x2, f2; x2 = a + g * (c - a); f2 = val(n, r1, s, math.exp(x2), w)
        bopt = math.exp((a + c) / 2)
        bstr = f'{bopt:.5g}'
        cand.append((val(n, r1, s, float(bstr), w), s, bstr))
    return max(cand)

def crossover(r1mode, w, target):
    # smallest n (r1 = 0 formal or r1 = n) whose optimal float bound exceeds log(target) + 1e-6
    def ok(n):
        v, s, b = best(n, n if r1mode == 'tr' else 0, w)
        return v > math.log(target) + 1e-6
    lo, hi = 2, 2
    while not ok(hi):
        lo, hi = hi, hi * 2
        if hi > 10 ** 9: return None
    while hi - lo > 1:
        m = (lo + hi) // 2
        if ok(m): hi = m
        else: lo = m
    return hi

if __name__ == '__main__':
    out = []
    for wk, (wstr, w) in W.items():
        for n in NS:
            r1s = sorted({0, n % 2, n} | (set(range(n % 2, n + 1, 2)) if n <= MIXED_MAX else set()))  # r1 = 0 formal (all signatures)
            for r1 in r1s:
                v, s, b = best(n, r1, w)
                out.append({'w': wk, 'wq': wstr, 'n': n, 'r1': r1, 'shape': s, 'b': b, 'float_rd': math.exp(v), 'kind': 'table'})
        print(wk, 'done', flush=True)
    cl_c, cl_r = 4 * math.pi * math.exp(0.57721566490153286), 4 * math.pi * math.exp(1.57721566490153286)
    for wk in ('ct', 'ab'):
        wstr, w = W[wk]
        for mode, target in (('tc', cl_c), ('tr', cl_r)):
            n0 = crossover(mode, w, target)
            r1 = 0 if mode == 'tc' else n0
            v, s, b = best(n0, r1, w)
            out.append({'w': wk, 'wq': wstr, 'n': n0, 'r1': r1, 'shape': s, 'b': b, 'float_rd': math.exp(v), 'kind': 'cross_' + mode})
            print(wk, mode, 'crossover n0 =', n0, 'rd>', math.exp(v), flush=True)
    # threshold degrees quoted in REPORT.md: (strip, mode, target)
    for wk, mode, target in (('ct', 'tc', 22.0), ('ct', 'tr', 60.0), ('ab', 'tr', 81.0)):
        wstr, w = W[wk]
        n0 = crossover(mode, w, target)
        r1 = 0 if mode == 'tc' else n0
        v, s, b = best(n0, r1, w)
        out.append({'w': wk, 'wq': wstr, 'n': n0, 'r1': r1, 'shape': s, 'b': b, 'float_rd': math.exp(v), 'kind': 'thr_' + mode, 'target': target})
        print(wk, mode, 'rd >', target, 'for n >=', n0, flush=True)
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else 'kernels.json', 'w'), indent=0)
    print('entries', len(out))
