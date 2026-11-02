#!/usr/bin/env python3
"""verify.py -- rigorous re-check of the lower bound
        t_Z([0,1]) >= 0.4222354990      (integer Chebyshev constant of [0,1])
from CERTIFICATE_LOWER_new.json, via Theorem "cert" of arXiv 2610.00386:
  t_Z >= exp( min(inf_[0,1] A, min_i G_i) / M ),   M = exact total mass, needs B < M,
  A(x) = U_mu(x) - sum_j c_j log|Q_j(x)|,
  G_i  = (1/d_i)[ int log|Q_i| dmu - sum_{j!=i} c_j log|Res(Q_i,Q_j)| ].
Checks: (1) exact (Fractions of the shipped doubles): edges increasing in [0,1], masses/multipliers >= 0,
M and B, library irreducible/primitive/distinct with leading coefficient >= 1, exact x->1-x symmetry
(then inf A over [0,1] = inf over [0,1/2]); (2) gcert.c: G_i in arb balls, exact integer resultants;
(3) lowcert.c: inf A >= V by bisection in arb balls; (4) final exp by mpmath interval arithmetic.
Optional: --upper FILE re-checks the published upper certificate (arXiv 2610.00386 ancillary JSON),
          --lower FILE re-checks the published lower certificate (same theorem).
Needs: cc, FLINT >= 3 (with arb), GMP, MPFR (e.g. Homebrew), python3 with mpmath and sympy
(python-flint used if present).  Runtime about 1 minute.
"""
import json, os, subprocess, sys, tempfile, math, time, importlib
from fractions import Fraction as Fr
import mpmath as mp


def exp_lower(x):
    """rigorous lower bound of exp(x), x rational: exp is increasing, round x down first."""
    mp.iv.dps = 30
    return mp.mpf(mp.iv.exp(mp.iv.mpf("%de-40" % math.floor(x * 10 ** 40))).a)


def exp_upper(x):
    """rigorous upper bound of exp(x), x rational."""
    mp.iv.dps = 30
    return mp.mpf(mp.iv.exp(mp.iv.mpf("%de-40" % math.ceil(x * 10 ** 40))).b)

HERE = os.path.dirname(os.path.abspath(__file__))


def compile_c(name, tmp):
    exe = os.path.join(tmp, name)
    flags = ["-I/opt/homebrew/include", "-I/usr/local/include", "-L/opt/homebrew/lib", "-L/usr/local/lib"]
    subprocess.run(["cc", "-O2", "-o", exe, os.path.join(HERE, name + ".c")] + flags +
                   ["-lflint", "-lgmp", "-lmpfr", "-lm"], check=True)
    return exe


def irreducible_primitive(q):
    try:
        flint = importlib.import_module("flint")
        cont, fac = flint.fmpz_poly(list(reversed(q))).factor()
        return abs(int(cont)) == 1 and len(fac) == 1 and fac[0][1] == 1
    except ImportError:
        import sympy as sp
        x = sp.Symbol("x")
        res = sp.factor_list(sp.Poly(q, x).as_expr())
        c, fl = res[0], res[1]
        return abs(c) == 1 and sum(m for _, m in fl) == 1


def reflect(q):
    asc = list(reversed(q)); n = len(asc) - 1; out = [0] * (n + 1)
    for k, a in enumerate(asc):
        for i in range(k + 1):
            out[i] += a * math.comb(k, i) * (-1) ** i
    d = list(reversed(out))
    return tuple(-v for v in d) if d[0] < 0 else tuple(d)


def exact_checks(edges, masses, lib, mults):
    E = [Fr(e) for e in edges]; Ms = [Fr(m) for m in masses]; Cs = [Fr(c) for c in mults]
    assert all(E[i] < E[i + 1] for i in range(len(E) - 1)) and E[0] >= 0 and E[-1] <= 1, "edges"
    assert min(Ms) >= 0 and min(Cs) >= 0, "negative mass or multiplier"
    seen = set()
    for q in lib:
        assert q[0] >= 1, "leading coefficient"
        assert tuple(q) not in seen, "duplicate"
        seen.add(tuple(q))
        assert irreducible_primitive(q), ("reducible or imprimitive", q)
    M = sum(Ms); B = sum(c * (len(q) - 1) for c, q in zip(Cs, lib))
    assert B < M, "budget B < M fails"
    K = len(Ms)
    idx = {tuple(q): i for i, q in enumerate(lib)}
    sym = all(E[k] + E[K - k] == 1 for k in range(K + 1)) and all(Ms[k] == Ms[K - 1 - k] for k in range(K))
    sym = sym and all(reflect(q) in idx and Cs[idx[reflect(q)]] == Cs[i] for i, q in enumerate(lib))
    return M, B, sym


def write_lists(f, edges, masses):
    f.write("%d\n" % len(masses))
    f.write(" ".join("%.17g" % e for e in edges) + "\n")
    f.write(" ".join("%.17g" % m for m in masses) + "\n")


def check_lower(path, tmp, low, gc, label):
    C = json.load(open(path))
    edges, masses, lib, mults, V = C["cell_edges"], C["masses"], [list(map(int, q)) for q in C["library"]], C["multipliers"], C["V"]
    t0 = time.time()
    M, B, sym = exact_checks(edges, masses, lib, mults)
    print("[%s] exact checks OK: K=%d J=%d  M-1=%.3e  B-M=%.3e  symmetric=%s (%.0fs)" % (
        label, len(masses), len(lib), float(M - 1), float(B - M), sym, time.time() - t0), flush=True)
    gi = os.path.join(tmp, "g.txt")
    with open(gi, "w") as f:
        write_lists(f, edges, masses)
        f.write("%d\n" % len(lib))
        for q, c in zip(lib, mults):
            f.write("%.17g %d %s\n" % (c, len(q) - 1, " ".join(str(a) for a in reversed(q))))
    out = subprocess.run([gc], stdin=open(gi), capture_output=True, text=True, check=True).stdout.split("\n")
    gmin = float([l for l in out if l.startswith("MIN")][0].split()[1])
    print("[%s] min_i G_i >= %.15f (arb)" % (label, gmin), flush=True)
    li = os.path.join(tmp, "a.txt")
    X1 = 0.5 if sym else 1.0
    with open(li, "w") as f:
        write_lists(f, edges, masses)
        act = [(q, c) for q, c in zip(lib, mults) if c > 0]
        f.write("%d\n" % len(act))
        for q, c in act:
            f.write("%.17g %d %s\n" % (c, len(q) - 1, " ".join(str(a) for a in reversed(q))))
        f.write("0 %.17g %.17g 45\n" % (X1, V))
    t0 = time.time()
    r = subprocess.run([low], stdin=open(li), capture_output=True, text=True)
    print("[%s] lowcert on [0,%g]: %s (%.0fs)" % (label, X1, r.stdout.strip(), time.time() - t0), flush=True)
    assert r.returncode == 0 and r.stdout.startswith("OK"), "inf A >= V NOT proved"
    m = min(Fr(V), Fr(gmin))
    lo = exp_lower(m / M)
    print("[%s] PROVED t_Z([0,1]) >= %s" % (label, mp.nstr(lo, 15, min_fixed=-1, max_fixed=1)), flush=True)
    return lo


def check_upper(path, tmp, up):
    C = json.load(open(path))
    fi = os.path.join(tmp, "u.txt")
    with open(fi, "w") as f:
        f.write("%d\n" % len(C["factors"]))
        for it in C["factors"]:
            asc = list(reversed(it["coeffs_desc_w"]))
            f.write("%d %d %s\n" % (it["e"], len(asc) - 1, " ".join(str(int(a)) for a in asc)))
        f.write("%d\n0 0.25 %.17g 50\n" % (C["N_w_degree"], C["proved_threshold_T"]))
    r = subprocess.run([up], stdin=open(fi), capture_output=True, text=True)
    print("[upper] upcert: %s" % r.stdout.strip(), flush=True)
    assert r.returncode == 0 and r.stdout.startswith("OK")
    T = Fr(C["proved_threshold_T"])
    print("[upper] PROVED t_Z([0,1]) <= %s" % mp.nstr(exp_upper(T / 2), 15), flush=True)


def main():
    args = sys.argv[1:]
    with tempfile.TemporaryDirectory() as tmp:
        low = compile_c("lowcert", tmp); gc = compile_c("gcert", tmp)
        check_lower(os.path.join(HERE, "CERTIFICATE_LOWER_new.json"), tmp, low, gc, "new lower")
        if "--lower" in args:
            check_lower(args[args.index("--lower") + 1], tmp, low, gc, "published lower")
        if "--upper" in args:
            check_upper(args[args.index("--upper") + 1], tmp, compile_c("upcert", tmp))
    return 0


if __name__ == "__main__":
    sys.exit(main())
