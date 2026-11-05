#!/usr/bin/env python3
"""Adversarial dominance test of the atom-coordinate C kernel (evalbox, built from the deposited
source): on boxes around local maxima of |a_n| and random boxes, the reported V and P bounds must
dominate max |a_n| over the box.  The max is searched by bounded local optimization (float64) from
several starts and re-evaluated in 256-bit Arb (log-series recursion, independent of the kernel's
product/Taylor-model formula).  usage: kz_vtest.py evalbox n Nmax boxes_per_N seed"""
import sys, math, json, random, subprocess, os
import numpy as np
from scipy.optimize import minimize
from flint import arb, acb, ctx
ctx.prec = 256
TWO_E = 2 / math.e

def an_f(x, n, N):
    lam = np.asarray(x[:N]); th = np.concatenate([[0.0], np.asarray(x[N:])])
    h = [0j] + [-2 * np.sum(lam * np.exp(1j * m * th)) for m in range(1, n + 1)]
    f = [1 + 0j]
    for k in range(1, n + 1):
        f.append(sum(j * h[j] * f[k - j] for j in range(1, k + 1)) / k)
    return abs(f[n]) * math.exp(-lam.sum())

def an_arb(x, n, N):
    mass = [arb(v) for v in x[:N]]; ang = [arb(0)] + [arb(v) for v in x[N:]]
    tot = sum(mass, arb(0))
    g = [acb(0)] + [-2 * sum((m * acb((k * t).cos(), (k * t).sin()) for m, t in zip(mass, ang)), acb(0)) for k in range(1, n + 1)]
    a = [acb(1)]
    for k in range(1, n + 1):
        a.append(sum((j * g[j] * a[k - j] for j in range(1, k + 1)), acb(0)) / k)
    return abs(a[n]) * (-tot).exp()

def normalize(x, N):
    """rotate heaviest atom to angle 0, sort others by angle in (0,2pi), conjugate if th2+thN > 2pi."""
    lam = list(x[:N]); th = [0.0] + list(x[N:])
    h = max(range(N), key=lambda j: lam[j])
    th = [(t - th[h]) % (2 * math.pi) for t in th]
    rest = sorted([j for j in range(N) if j != h], key=lambda j: th[j])
    lam2 = [lam[h]] + [lam[j] for j in rest]; th2 = [th[j] for j in rest]
    if N >= 3 and th2[0] + th2[-1] > 2 * math.pi:
        pairs = sorted(zip([(2 * math.pi - t) % (2 * math.pi) for t in th2], lam2[1:]))
        th2 = [a for a, b in pairs]; lam2 = [lam2[0]] + [b for a, b in pairs]
    if N == 2 and th2[0] > math.pi:
        th2 = [2 * math.pi - th2[0]]
    return np.array(lam2 + th2)

def local_max(n, N, starts, rng):
    out = []
    for _ in range(starts):
        x0 = np.concatenate([rng.uniform(0.05, 1.5, N), rng.uniform(0, 2 * math.pi, N - 1)])
        r = minimize(lambda x: -an_f(np.concatenate([np.abs(x[:N]), x[N:]]), n, N), x0, method='BFGS', options=dict(maxiter=400))
        x = np.concatenate([np.abs(r.x[:N]), np.mod(r.x[N:], 2 * math.pi)])
        out.append((-r.fun, normalize(x, N)))
    out.sort(key=lambda t: -t[0])
    return out

if __name__ == '__main__':
    exe, n, Nmax, per, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    thi = {6: 6.99, 7: 8.15}[n]
    rng = np.random.default_rng(seed)
    report = dict(n=n, boxes=0, finite_v=0, v_closed=0, violations=[], min_rel_slack_v=None, min_rel_slack_vclosed=None,
                  max_value_in_vclosed=0.0)
    for N in range(1, Nmax + 1):
        lm = local_max(n, N, 30 if N > 1 else 10, rng)
        boxes = []
        for b in range(per):
            if b % 2 == 0:
                c = lm[b % min(6, len(lm))][1].copy()
            else:
                c = normalize(np.concatenate([rng.uniform(0.02, 2.5, N), rng.uniform(0, 2 * math.pi, N - 1)]), N)
            rad = 10 ** rng.uniform(-6, -1.2, size=2 * N - 1)
            lo = np.concatenate([np.maximum(0, c[:N] - rad[:N]), c[N:] - rad[N:]])
            hi = np.concatenate([c[:N] + rad[:N], c[N:] + rad[N:]])
            hi = np.minimum(hi, 7.99); lo = np.maximum(lo, np.concatenate([np.zeros(N), -7.99 * np.ones(N - 1)]))
            boxes.append((lo, hi))
        lines = [f'{n} {N} {thi} ' + ' '.join(format(v, '.17g') for v in list(lo) + list(hi)) for lo, hi in boxes]
        res = subprocess.run([exe], input=chr(10).join(lines) + chr(10), text=True, capture_output=True, env=dict(os.environ, KRZ_P='1'), check=True)
        for (lo, hi), line in zip(boxes, res.stdout.splitlines(), strict=True):
            st, vb, pb, *_ = line.split(); vb = float(vb); pb = float(pb); report['boxes'] += 1
            best = 0.0; bestx = None
            for s in range(6):
                x0 = lo + rng.uniform(0, 1, len(lo)) * (hi - lo) if s else (lo + hi) / 2
                r = minimize(lambda x: -an_f(x, n, N), x0, method='L-BFGS-B', bounds=list(zip(lo, hi)), options=dict(maxiter=200))
                v = -r.fun
                if v > best: best, bestx = v, np.clip(r.x, lo, hi)
            val = an_arb(list(bestx), n, N)
            for name, bound in (('V', vb), ('P', pb)):
                if math.isfinite(bound):
                    if not (val < arb(bound)):
                        report['violations'].append(dict(N=N, test=name, bound=bound, value=float(val), lo=list(lo), hi=list(hi)))
            if math.isfinite(vb):
                report['finite_v'] += 1
                slack = (vb - float(val)) / vb
                report['min_rel_slack_v'] = slack if report['min_rel_slack_v'] is None else min(report['min_rel_slack_v'], slack)
                if vb < TWO_E:
                    report['v_closed'] += 1
                    report['max_value_in_vclosed'] = max(report['max_value_in_vclosed'], float(val))
                    report['min_rel_slack_vclosed'] = slack if report['min_rel_slack_vclosed'] is None else min(report['min_rel_slack_vclosed'], slack)
    report['violations'] = report['violations'][:10]
    print(json.dumps(report))
