#!/usr/bin/env python3
"""batch.py - run named symmetric constructions through hyp + SAT with a time limit each.
usage: batch.py TIMELIMIT name1 name2 ...   (names defined in CONS below; 'list' prints them)"""
import sys, os, json, time, threading
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sym import *
from grp import *
HERE = os.path.dirname(os.path.abspath(__file__))
KK = int(os.environ.get('XI6_K', '6'))
os.makedirs(os.path.join(HERE, 'sol'), exist_ok=True)

def cay(Tf): return lambda: (lambda T: (len(T), left_regular(T)))(Tf())
def bc(mod, r, fixed=0, mults=()): return lambda: block_circ(mod, r, fixed, mults)
def cycsd(nn, a, d): return lambda: cyclic_sd(nn, a, d)
CONS = {}
# Cayley groups
for nm, f in [
    ('Z9xZ3', lambda: abelian([9, 3])), ('Z3^3', lambda: abelian([3, 3, 3])), ('Heis3', lambda: heis(3)), ('Z9sdZ3', lambda: cyclic_sd(9, 4, 3)),
    ('Z14xZ2', lambda: abelian([14, 2])), ('D14', lambda: dihedral(14)), ('Dic7', lambda: dicyclic(7)), ('Z7sdZ4', lambda: cyclic_sd(7, 6, 4)),
    ('D15', lambda: dihedral(15)), ('Z5xS3', lambda: product(abelian([5]), S3())), ('Z3xD5', lambda: product(abelian([3]), dihedral(5))),
    ('Z16xZ2', lambda: abelian([16, 2])), ('Z8xZ4', lambda: abelian([8, 4])), ('Z8xZ2^2', lambda: abelian([8, 2, 2])), ('Z4^2xZ2', lambda: abelian([4, 4, 2])),
    ('Z4xZ2^3', lambda: abelian([4, 2, 2, 2])), ('Z2^5', lambda: abelian([2, 2, 2, 2, 2])), ('D16', lambda: dihedral(16)), ('Dic8', lambda: dicyclic(8)),
    ('SD32', lambda: cyclic_sd(16, 7, 2)), ('M32', lambda: cyclic_sd(16, 9, 2)), ('Z8sdZ4_3', lambda: cyclic_sd(8, 3, 4)), ('Z8sdZ4_5', lambda: cyclic_sd(8, 5, 4)),
    ('Z8sdZ4_7', lambda: cyclic_sd(8, 7, 4)), ('Z4sdZ8', lambda: cyclic_sd(4, 3, 8)), ('D8xZ2', lambda: product(dihedral(8), abelian([2]))),
    ('Dic4xZ2', lambda: product(dicyclic(4), abelian([2]))), ('D4xZ4', lambda: product(dihedral(4), abelian([4]))), ('D4xZ2^2', lambda: product(dihedral(4), abelian([2, 2]))),
    ('Q8xZ4', lambda: product(dicyclic(2), abelian([4]))), ('Q8xZ2^2', lambda: product(dicyclic(2), abelian([2, 2]))), ('Dih(Z4xZ4)', lambda: gen_dihedral([4, 4])),
    ('Dih(Z8xZ2)', lambda: gen_dihedral([8, 2])), ('D4xD4/?', lambda: product(dihedral(4), dihedral(2))),
    ('D17', lambda: dihedral(17)),
    ('Z6xZ6', lambda: abelian([6, 6])), ('Z12xZ3', lambda: abelian([12, 3])), ('Z18xZ2', lambda: abelian([18, 2])), ('D18', lambda: dihedral(18)), ('Dic9', lambda: dicyclic(9)),
    ('Z3xA4', lambda: product(abelian([3]), A4())), ('S3xS3', lambda: product(S3(), S3())), ('Z6xS3', lambda: product(abelian([6]), S3())), ('Z3xDic3', lambda: product(abelian([3]), dicyclic(3))),
    ('Dih(Z3xZ6)', lambda: gen_dihedral([3, 6])), ('Z3^2sdZ4', lambda: semidirect([3, 3], [[0, 2], [1, 0]], 4)), ('Z2xZ2xZ9', lambda: abelian([2, 2, 9])), ('Z9sdZ4', lambda: cyclic_sd(9, 8, 4)),
    ('D19', lambda: dihedral(19)), ('Z13sdZ3', lambda: cyclic_sd(13, 3, 3)),
    ('Z20xZ2', lambda: abelian([20, 2])), ('Z10xZ2^2', lambda: abelian([10, 2, 2])), ('D20', lambda: dihedral(20)), ('Dic10', lambda: dicyclic(10)),
    ('Z5sdZ8', lambda: cyclic_sd(5, 2, 8)), ('F20xZ2', lambda: product(cyclic_sd(5, 2, 4), abelian([2]))), ('Z4xD5', lambda: product(abelian([4]), dihedral(5))),
    ('Z5xD4', lambda: product(abelian([5]), dihedral(4))), ('Z5xQ8', lambda: product(abelian([5]), dicyclic(2))), ('Z2^2xD5', lambda: product(abelian([2, 2]), dihedral(5))),
    ('Dic5xZ2', lambda: product(dicyclic(5), abelian([2]))), ('Z5sdZ8_4', lambda: cyclic_sd(5, 4, 8)),
    ('D10', lambda: dihedral(10)), ('Dic5', lambda: dicyclic(5)), ('F20', lambda: cyclic_sd(5, 2, 4)), ('Z10xZ2', lambda: abelian([10, 2])), ('Z7sdZ3', lambda: cyclic_sd(7, 2, 3)), ('D11', lambda: dihedral(11)), ('D12', lambda: dihedral(12)), ('Dic6', lambda: dicyclic(6)), ('Z12xZ2', lambda: abelian([12, 2])), ('Z6xZ2^2', lambda: abelian([6, 2, 2])), ('SL23', lambda: None), ('Z3xQ8', lambda: product(abelian([3]), dicyclic(2))), ('Z3xD4', lambda: product(abelian([3]), dihedral(4))), ('S3xZ4', lambda: product(S3(), abelian([4]))), ('S3xZ2^2', lambda: product(S3(), abelian([2, 2]))), ('Z2xA4', lambda: product(abelian([2]), A4())), ('Z3sdZ8', lambda: cyclic_sd(3, 2, 8)), ('Dic3xZ2', lambda: product(dicyclic(3), abelian([2]))), ('Z11xZ2', lambda: abelian([22])),
    ('Z5^2', lambda: abelian([5, 5])), ('Z7^2', lambda: abelian([7, 7])), ('Z2^6', lambda: abelian([2]*6)), ('Z4^3', lambda: abelian([4, 4, 4])), ('Z8^2', lambda: abelian([8, 8])), ('Z3^4', lambda: abelian([3]*4)), ('Z6xZ6', lambda: abelian([6, 6])), ('Z5sdZ5?', lambda: abelian([25])), ('D13', lambda: dihedral(13)), ('Z3^2xZ3?', lambda: abelian([3, 9])),
    ]:
    CONS['cay_' + nm] = cay(f)
for nn in range(20, 80):
    CONS['cyc_%d' % nn] = (lambda nn: (lambda: (nn, cyc_gen(nn))))(nn)
for m in range(9, 22):
    for r in (2, 3, 4):
        if 26 <= r * m <= 44: CONS[f'bc_{m}x{r}'] = bc(m, r)
    for r in (2, 3):
        for fx in (1, 2):
            if 26 <= r * m + fx <= 44: CONS[f'bc_{m}x{r}+{fx}'] = bc(m, r, fx)

def run(name, tl):
    t0 = time.time()
    n, gens = CONS[name]()
    C, m, roots = orbits_from_gens(n, gens)
    E, info = hyperedges(n, C, m, roots, KK)
    t1 = time.time()
    from pysat.solvers import Solver
    v = lambda c, col: 3 * c + col + 1
    cl = []
    for c in range(m):
        cl.append([v(c, 0), v(c, 1), v(c, 2)])
        for a, b in ((0, 1), (0, 2), (1, 2)): cl.append([-v(c, a), -v(c, b)])
    triv = any(len(T) < 3 for T in E)
    res = 'UNSAT(trivial)' if triv else None
    col = None
    if not triv:
        for T in E:
            for cc in range(3): cl.append([v(c, cc) for c in T])
        cl.append([v(0, 0)])
        if m > 1: cl.append([-v(1, 2)])
        s = Solver(name='cadical153', bootstrap_with=cl)
        timer = threading.Timer(tl, lambda: s.interrupt()); timer.start()
        r = s.solve_limited(expect_interrupt=True)
        timer.cancel()
        if r is True:
            model = set(x for x in s.get_model() if x > 0)
            col = [next(cc for cc in range(3) if v(c, cc) in model) for c in range(m)]
            res = 'SAT'
        elif r is False: res = 'UNSAT'
        else: res = 'TIMEOUT'
        s.delete()
    rec = dict(k=KK, name=name, n=n, vars=m, roots=len(roots), hyper=len(E), res=res, t_gen=round(t1 - t0, 1), t_sat=round(time.time() - t1, 1))
    if col is not None:
        M = full_colouring(n, C, col)
        ok, miss = check_fast(M, KK)
        rec['check'] = ok
        write_colouring(os.path.join(HERE, 'sol', ('k%d_' % KK if KK != 6 else '') + name.replace('/', '_') + '.txt'), M, note=name)
    print(json.dumps(rec), flush=True)
    with open(os.path.join(HERE, 'results.jsonl'), 'a') as f: f.write(json.dumps(rec) + '\n')


def dump(name, path):
    import pickle
    n, gens = CONS[name]()
    C, m, roots = orbits_from_gens(n, gens)
    inp = f"{n} {m} 6 {len(roots)}\n" + " ".join(map(str, roots)) + "\n" + "\n".join(" ".join(map(str, row)) for row in C) + "\n"
    out = subprocess.run([os.path.join(HERE, "hyp")], input=inp, capture_output=True, text=True, check=True)
    open(path + '.hyp', 'w').write(out.stdout)
    pickle.dump((n, C, m), open(path + '.pkl', 'wb'))
    print(out.stderr.strip())

def load_sol(path, solline, outfn, note=''):
    import pickle
    n, C, m = pickle.load(open(path + '.pkl', 'rb'))
    col = list(map(int, solline.split()[1:]))
    M = full_colouring(n, C, col)
    ok = check_fast(M, 6)
    write_colouring(outfn, M, note=note)
    return ok


if __name__ == '__main__':
    if sys.argv[1] == 'dump':
        dump(sys.argv[2], sys.argv[3]); sys.exit()
    if sys.argv[1] == 'load':
        lines = open(sys.argv[3]).read().split('\n')
        sol = next(l for l in lines if l.startswith('SOL'))
        print(load_sol(sys.argv[2], sol, sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else '')); sys.exit()

    tl = float(sys.argv[1]); names = sys.argv[2:]
    if names == ['list']:
        print(' '.join(CONS)); sys.exit()
    for nm in names:
        if nm.endswith('*'):
            for k in CONS:
                if k.startswith(nm[:-1]): run(k, tl)
        else: run(nm, tl)

