#!/usr/bin/env python3
"""fam.py - symmetric search with colour-rotating symmetries.
Generators are (vertex permutation, d) meaning colour(g(e)) = colour(e) + d (mod 3).
usage: fam.py <spec> <mode: sat|ls|both> <seconds>
spec examples:  A:m:f      Z_m x Z_3 orbits, sigma:(x,i)->(x,i+1) rotating colours, f extra vertices
                B:m:a:f    sigma:(x,i)->(a x, i+1) rotating colours (a^3=1 mod m)
                P:m:r:f    plain block circulant r orbits (no colour symmetry)
                E:m:r:a:f  block circulant r orbits + multiplier x->a x on all orbits rotating colours (ord(a)|3? any with a^3=1)
f fixed vertices: in A/B they are permuted by sigma in one 3-cycle if f==3, fixed if f==1 (then edges to them rotate)."""
import sys, os, subprocess, itertools, time, json, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
KK = int(os.environ.get('XI6_K', '6'))
sys.path.insert(0, HERE)
from sym import write_colouring, check_fast

def wuf_classes(n, gens):
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    idx = {p: i for i, p in enumerate(pairs)}
    P = list(range(len(pairs))); O = [0] * len(pairs)
    def find(x):
        o = 0; path = []
        while P[x] != x:
            path.append(x); x = P[x]
        # path compression
        r = x
        # recompute offsets
        acc = 0
        for y in reversed(path):
            acc = (acc + O[y]) % 3
        # second pass: set each to root with its total offset
        tot = 0; offs = []
        for y in path:
            pass
        # simple non-compressing walk for offsets
        return r
    def findo(x):
        o = 0
        while P[x] != x:
            o = (o + O[x]) % 3; x = P[x]
        return x, o
    for (g, d) in gens:
        for i, (a, b) in enumerate(pairs):
            x, y = g[a], g[b]
            if x > y: x, y = y, x
            j = idx[(x, y)]
            r1, o1 = findo(i); r2, o2 = findo(j)
            if r1 == r2:
                if (o2 - o1 - d) % 3 != 0: return None
            else:
                P[r2] = r1; O[r2] = (o1 + d - o2) % 3
    roots = {}; C = [[-1] * n for _ in range(n)]
    for i, (a, b) in enumerate(pairs):
        r, o = findo(i)
        if r not in roots: roots[r] = len(roots)
        C[a][b] = C[b][a] = 3 * roots[r] + o
    # vertex orbits
    VP = list(range(n))
    def vf(x):
        while VP[x] != x: x = VP[x]
        return x
    for (g, d) in gens:
        for a in range(n):
            ra, rb = vf(a), vf(g[a])
            if ra != rb: VP[max(ra, rb)] = min(ra, rb)
    vroots = sorted({vf(a) for a in range(n)})
    return C, len(roots), vroots

def build(spec):
    t = spec.split(':'); kind = t[0]
    if kind in ('A', 'B'):
        m = int(t[1]); a = 1 if kind == 'A' else int(t[2]); f = int(t[-1])
        n = 3 * m + f
        V = lambda x, i: (i % 3) * m + (x % m)
        tr = list(range(n)); sg = list(range(n))
        for i in range(3):
            for x in range(m):
                tr[V(x, i)] = V(x + 1, i); sg[V(x, i)] = V(a * x, i + 1)
        if f == 3:
            sg[3 * m], sg[3 * m + 1], sg[3 * m + 2] = 3 * m + 1, 3 * m + 2, 3 * m
        gens = [(tr, 0), (sg, 1)]
    elif kind == 'P':
        m, r, f = int(t[1]), int(t[2]), int(t[3]); n = r * m + f
        tr = list(range(n))
        for i in range(r):
            for x in range(m): tr[i * m + x] = i * m + (x + 1) % m
        gens = [(tr, 0)]
    elif kind == 'E':
        m, r, a, f = int(t[1]), int(t[2]), int(t[3]), int(t[4]); n = r * m + f
        tr = list(range(n)); mu = list(range(n))
        for i in range(r):
            for x in range(m): tr[i * m + x] = i * m + (x + 1) % m; mu[i * m + x] = i * m + (a * x) % m
        gens = [(tr, 0), (mu, 1)]
    elif kind == 'F':  # E but multiplier colour-preserving (d=0)
        m, r, a, f = int(t[1]), int(t[2]), int(t[3]), int(t[4]); n = r * m + f
        tr = list(range(n)); mu = list(range(n))
        for i in range(r):
            for x in range(m): tr[i * m + x] = i * m + (x + 1) % m; mu[i * m + x] = i * m + (a * x) % m
        gens = [(tr, 0), (mu, 0)]
    else:
        raise SystemExit('bad spec')
    return n, gens

def main():
    spec, mode, secs = sys.argv[1], sys.argv[2], float(sys.argv[3])
    n, gens = build(spec)
    res = wuf_classes(n, gens)
    if res is None:
        print(json.dumps(dict(spec=spec, n=n, res='INCONSISTENT'))); return
    C, m, vroots = res
    plain = all(d == 0 for (_, d) in gens)
    if plain:
        return main_plain(spec, mode, secs, n, C, m, vroots)
    if 3 * m > 128:
        print(json.dumps(dict(spec=spec, n=n, m=m, res='TOO_MANY_CLASSES'))); return
    inp = f"{n} {3*m} {KK} {len(vroots)}\n" + " ".join(map(str, vroots)) + "\n" + "\n".join(" ".join(map(str, row)) for row in C) + "\n"
    out = subprocess.run([os.path.join(HERE, 'hyp')], input=inp, capture_output=True, text=True, check=True)
    lines = out.stdout.strip().split('\n')
    E = [list(map(int, l.split())) for l in lines[1:]]
    E2 = []
    for T in E:
        ks = {}
        for it in T: ks.setdefault(it // 3, set()).add(it % 3)
        if any(len(s) == 3 for s in ks.values()): continue
        E2.append(T)
    tag = spec.replace(':', '_')
    os.makedirs(os.path.join(HERE, 'inst'), exist_ok=True)
    hp = os.path.join(HERE, 'inst', 'f%d_' % KK + tag)
    with open(hp + '.hyp', 'w') as fh:
        fh.write(f"{3*m} {len(E2)}\n" + "\n".join(" ".join(map(str, T)) for T in E2) + "\n")
    pickle.dump((n, C, m), open(hp + '.pkl', 'wb'))
    info = dict(spec=spec, n=n, m=m, roots=len(vroots), hyper=len(E2), minsize=min((len(T) for T in E2), default=0))
    print('BUILT', json.dumps(info), flush=True)
    col = None; res = None
    if mode in ('sat', 'both'):
        from pysat.solvers import Solver
        v = lambda k, c: 3 * k + c + 1
        cl = []
        for k in range(m):
            cl.append([v(k, 0), v(k, 1), v(k, 2)])
            for x, y in ((0, 1), (0, 2), (1, 2)): cl.append([-v(k, x), -v(k, y)])
        for T in E2:
            for c in range(3): cl.append(sorted({v(it // 3, (c - it % 3) % 3) for it in T}))
        cl.append([v(0, 0)])
        t0 = time.time()
        with Solver(name=os.environ.get('XI6_SOLVER','cadical153'), bootstrap_with=cl) as s:
            ok = s.solve()
            if ok:
                mod = set(x for x in s.get_model() if x > 0)
                col = [next(c for c in range(3) if v(k, c) in mod) for k in range(m)]
                res = 'SAT'
            else: res = 'UNSAT'
        info['t_sat'] = round(time.time() - t0, 1)
    if col is None and mode in ('ls', 'both') and res != 'UNSAT':
        o = subprocess.run([os.path.join(HERE, 'ls2'), hp + '.hyp', str(secs), '1', '1', str(m)], capture_output=True, text=True)
        first = o.stdout.split('\n')[0]
        info['ls'] = first
        if o.stdout.startswith('SOLVED'):
            col = list(map(int, o.stdout.split('\n')[1].split()[1:])); res = 'SAT(ls)'
        else: res = res or 'LSFAIL'
    info['res'] = res
    if col is not None:
        M = [[-1 if a == b else (col[C[a][b] // 3] + C[a][b] % 3) % 3 for b in range(n)] for a in range(n)]
        ok, _ = check_fast(M, KK); info['check'] = ok
        os.makedirs(os.path.join(HERE, 'sol'), exist_ok=True)
        write_colouring(os.path.join(HERE, 'sol', 'f%d_' % KK + tag + '.txt'), M, note='fam ' + spec + ' k=%d' % KK)
    print(json.dumps(info), flush=True)
    with open(os.path.join(HERE, 'results.jsonl'), 'a') as fh: fh.write(json.dumps(info) + '\n')


def main_plain(spec, mode, secs, n, C, m, vroots):
    if m > 128:
        print(json.dumps(dict(spec=spec, n=n, m=m, res='TOO_MANY_CLASSES'))); return
    Cp = [[-1 if a == b else C[a][b] // 3 for b in range(n)] for a in range(n)]
    inp = f"{n} {m} {KK} {len(vroots)}\n" + " ".join(map(str, vroots)) + "\n" + "\n".join(" ".join(map(str, row)) for row in Cp) + "\n"
    out = subprocess.run([os.path.join(HERE, 'hyp')], input=inp, capture_output=True, text=True, check=True)
    lines = out.stdout.strip().split('\n')
    E = [list(map(int, l.split())) for l in lines[1:]]
    tag = spec.replace(':', '_')
    hp = os.path.join(HERE, 'inst', 'f%d_' % KK + tag)
    os.makedirs(os.path.join(HERE, 'inst'), exist_ok=True)
    with open(hp + '.hyp', 'w') as fh:
        fh.write(f"{m} {len(E)}\n" + "\n".join(" ".join(map(str, T)) for T in E) + "\n")
    pickle.dump((n, Cp, m), open(hp + '.pkl', 'wb'))
    info = dict(spec=spec, n=n, m=m, roots=len(vroots), hyper=len(E), minsize=min((len(T) for T in E), default=0), plain=True)
    print('BUILT', json.dumps(info), flush=True)
    col = None; res = None
    if any(len(T) < 3 for T in E): res = 'UNSAT(trivial)'
    if res is None and mode in ('sat', 'both'):
        from pysat.solvers import Solver
        v = lambda k, c: 3 * k + c + 1
        cl = []
        for k in range(m):
            cl.append([v(k, 0), v(k, 1), v(k, 2)])
            for x, y in ((0, 1), (0, 2), (1, 2)): cl.append([-v(k, x), -v(k, y)])
        for T in E:
            for c in range(3): cl.append([v(k, c) for k in T])
        cl.append([v(0, 0)])
        if m > 1: cl.append([-v(1, 2)])
        t0 = time.time()
        with Solver(name=os.environ.get('XI6_SOLVER','cadical153'), bootstrap_with=cl) as s:
            if s.solve():
                mod = set(x for x in s.get_model() if x > 0)
                col = [next(c for c in range(3) if v(k, c) in mod) for k in range(m)]; res = 'SAT'
            else: res = 'UNSAT'
        info['t_sat'] = round(time.time() - t0, 1)
    if col is None and mode in ('ls', 'both') and res is None:
        o = subprocess.run([os.path.join(HERE, 'ls2'), hp + '.hyp', str(secs), '1', '0', str(m)], capture_output=True, text=True)
        info['ls'] = o.stdout.split('\n')[0]
        if o.stdout.startswith('SOLVED'):
            col = list(map(int, o.stdout.split('\n')[1].split()[1:])); res = 'SAT(ls)'
        else: res = 'LSFAIL'
    info['res'] = res
    if col is not None:
        M = [[-1 if a == b else col[Cp[a][b]] for b in range(n)] for a in range(n)]
        ok, _ = check_fast(M, KK); info['check'] = ok
        write_colouring(os.path.join(HERE, 'sol', 'f%d_' % KK + tag + '.txt'), M, note='fam ' + spec + ' k=%d' % KK)
    print(json.dumps(info), flush=True)
    with open(os.path.join(HERE, 'results.jsonl'), 'a') as fh: fh.write(json.dumps(info) + '\n')

if __name__ == '__main__':
    main()

