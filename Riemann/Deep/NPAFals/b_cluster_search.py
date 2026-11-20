#!/usr/bin/env python3
"""npafals (B): can J off-line pairs fake a word for free? For a configuration X of J off-line pairs (depths s_j <= s0, angles th_j),
p_k = sum_j 2cosh(k s_j) e^{i k th_j}. A free fake of word W in a set S is a sequence with eps*v_X(W) -> const while eps*v_X(D) -> 0
for every diagonal word D in S, i.e. Q(X) = |v_X(W)| / max_D v_X(D) unbounded. We maximise Q over (s_j, th_j) by random-restart
L-BFGS for growing depth caps s0: Q bounded in s0 => no free fake from J pairs; Q ~ e^{c s0} => free fake.
Diagonal sets: 'P' (pair words |p_k|^2/k, k<=N; positive control: top triples are freely faked) and 'L4' (pair words and quartic
diagonals |p_d p_e|^2/z, d+e<=N). usage: b_cluster_search.py N J s0list words...   word = a,b[,c]:n  (e.g. 2,2,2:6)"""
import sys, math, itertools
import numpy as np
from scipy.optimize import minimize

def zl(lam):
    z = 1
    for j in set(lam):
        a = lam.count(j); z *= j ** a * math.factorial(a)
    return z

def diag_words(N, S):
    D = [((k,), zl((k,))) for k in range(1, N + 1)]
    if S == 'L4':
        D += [((d, e), zl((d, e))) for d in range(1, N) for e in range(d, N + 1 - d)]
    return D

def make_Q(N, J, s0, lam, n, S, T=40.0):
    D = diag_words(N, S); zW = math.sqrt(zl(lam) * zl((n,))); k = np.arange(N + 1)
    def pk(x):
        s = s0 / (1 + np.exp(-x[:J])); th = x[J:]; th = np.concatenate([[0.0], th[1:]])
        return (2 * np.cosh(np.outer(s, k)) * np.exp(1j * np.outer(th, k))).sum(axis=0)
    def logQ(x):
        p = pk(x)
        num = np.log(abs(np.prod([p[q] for q in lam]) * np.conj(p[n])) / zW + 1e-300)
        ld = np.array([2 * np.log(abs(np.prod([p[q] for q in d])) + 1e-300) - math.log(z) for d, z in D])
        m = ld.max(); den = m + np.log(np.exp(T * (ld - m)).sum()) / T      # smooth max (log-sum-exp, temperature T)
        return num - den, num - m
    return logQ, pk

def search(N, J, s0, lam, n, S, starts=30, rng=None):
    logQ, pk = make_Q(N, J, s0, lam, n, S); best = (-np.inf, None)
    for _ in range(starts):
        x0 = np.concatenate([rng.normal(2, 2, J), rng.uniform(0, 2 * np.pi, J)])
        r = minimize(lambda x: -logQ(x)[0], x0, method='L-BFGS-B', options={'maxiter': 200})
        v = logQ(r.x)[1]
        if v > best[0]: best = (v, r.x)
    s = s0 / (1 + np.exp(-best[1][:J]))
    return best[0], s, best[1][J:]

if __name__ == '__main__':
    N = int(sys.argv[1]); J = int(sys.argv[2]); s0s = [float(t) for t in sys.argv[3].split(',')]
    rng = np.random.default_rng(7)
    for w in sys.argv[4:]:
        S = 'L4'
        if w.startswith('P:'): S = 'P'; w = w[2:]
        l, n = w.split(':'); lam = tuple(int(t) for t in l.split(',')); n = int(n)
        out = []
        for s0 in s0s:
            v, s, th = search(N, J, s0, lam, n, S, rng=rng); out.append(f"s0={s0:g}: logQ={v:+.3f}")
        print(f"N={N} J={J} set={S} word={lam}|({n}) " + "  ".join(out) + f"  [last shape s={np.round(s,2).tolist()}]", flush=True)
