#!/usr/bin/env python3
"""Verify the certified upper bound V*(K_M) <= T for all M > u (3D snaky-corridor moving sofa).

Usage:  python3 verify.py [CERT_DIR]      (default: directory of this script)
CERT_DIR contains meta.json, codes.bin.gz (the branch-and-bound tree) and sofaub.c (the interval-arithmetic core).

Steps (all must pass):
 1. exact checks (fractions): c = 1 + u is exact, K^2 = 8u^2/(u^2-2)^2 (tilt), max(u, T) = T, sha256 of the tree;
 2. compile sofaub.c (-O2 -ffp-contract=off -fno-fast-math) and replay the tree: every leaf box must be infeasible
    (misses a window or violates the tilt) or carry a rigorous upper bound <= T; the tree must cover the root grid;
 3. independent cross-checks on a deterministic sample of leaves:
    a. exact rational re-integration (exact_check.py) of each sampled leaf's dilated-slice integral: the C bound
       must be >= the exact value (up to 1e-12 for the rounding of the edge vectors);
    b. Monte Carlo containment: random tubes inside sampled boxes have exact volume <= the leaf bound, and random
       tubes inside 'infeasible' boxes are indeed not admissible;
    c. root coverage: every admissible tube lies in some root box (exact rational arithmetic).
"""
import gzip, hashlib, json, math, os, random, subprocess, sys, tempfile, time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from exact_check import Dfam  # noqa: E402


def die(msg):
    print("FAIL:", msg)
    sys.exit(1)


def frac_sqrt_bounds(x, bits=80):
    """rational lo <= sqrt(x) <= hi for a nonnegative Fraction x"""
    if x == 0:
        return F(0), F(0)
    s = F(math.isqrt(int(x * 4**bits)), 2**bits)
    lo = s if s * s <= x else s - F(1, 2**bits)
    hi = lo + F(1, 2**bits)
    while lo * lo > x: lo -= F(1, 2**bits)
    while hi * hi < x: hi += F(1, 2**bits)
    return lo, hi


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else HERE
    meta = json.load(open(os.path.join(d, "meta.json")))
    u, T, ya = meta["u"], meta["target"], meta["ya"]
    grid = meta["grid"]
    print(f"certificate: u={u!r} target T={T!r} ya={ya!r} grid={grid}")
    # ---- 1. exact checks
    U, TT = F(u), F(T)
    C = 1 + U
    if not (float(C) - 1.0 == u and F(float(C)) == C): die("1+u not exact in double")
    K2 = 8 * U * U / (U * U - 2) ** 2
    mu0 = (U * U - 2) / (U * U + 2)
    if not (U * U > 2 and K2 == 1 / mu0**2 - 1): die("tilt constant mismatch")
    if not TT >= U: die("need T >= u so that max(u, V_ii) <= T")
    raw = gzip.open(os.path.join(d, "codes.bin.gz"), "rb").read()
    h = hashlib.sha256(raw).hexdigest()
    if h != meta["codes_sha256"]: die("tree sha256 mismatch")
    print(f"[1] exact: c = {C} exact; mu0 = {float(mu0):.12f}; K^2 = {float(K2):.12f}; T >= u; tree {len(raw):,} nodes, sha256 {h[:16]}... OK")
    # ---- 2. compile and replay
    tmp = tempfile.mkdtemp(prefix="sofaub_")
    exe = os.path.join(tmp, "sofaub")
    cc = subprocess.run(["cc", "-O2", "-ffp-contract=off", "-fno-fast-math", "-o", exe, os.path.join(d, "sofaub.c"), "-lm"],
                        capture_output=True, text=True)
    if cc.returncode != 0: die("compile failed: " + cc.stderr[:500])
    codes = os.path.join(tmp, "codes.bin")
    open(codes, "wb").write(raw)
    args = [str(u), repr(T), repr(ya)] + [str(g) for g in grid] + [codes]
    t0 = time.time()
    nleaf_hint = int(meta["leaves"])
    every = max(1, nleaf_hint // int(meta.get("n_sample", 300)))
    dump = os.path.join(tmp, "dump.txt")
    env = dict(os.environ, DUMP_EVERY=str(every), DUMP_FILE=dump)
    rv = subprocess.run([exe, "dump"] + args, capture_output=True, text=True, env=env)   # full replay + sample dump
    out = rv.stdout
    summ = [l for l in out.splitlines() if l.startswith("roots=")]
    if rv.returncode != 0 or "RESULT PASS" not in out or not summ: die("C replay failed:\n" + out[-2000:] + rv.stderr[-2000:])
    stats = dict(kv.split("=") for kv in summ[0].split())
    maxb = float(stats["max_leaf_bound"])
    if not maxb <= T: die("max leaf bound above target")
    if int(stats["nodes"]) != len(raw): die("node count mismatch")
    print(f"[2] replay: {int(stats['leaves']):,} leaves ({int(stats['infeasible']):,} infeasible, {int(stats['fallback']):,} near-vertical, "
          f"{int(stats['integral']):,} slice-integral, {int(stats['degenerate'])} degenerate); max leaf bound {maxb:.12f} <= T; {time.time()-t0:.0f}s OK")
    # ---- 3a. exact rational re-integration of sampled leaves
    rows = [list(map(float, l.split())) for l in open(dump)]
    if len(rows) < 50: die("too few sampled leaves")
    nint = nchk = 0; worst = 1e9
    for r in rows:
        st, v = int(r[1]), r[2]
        if st != 0: continue
        box = r[3:13]; kc = r[13:15]; xc = r[15:17]; gg = r[17:25]; PA0, PA1, PB0, PB1 = r[27:31]
        g1 = ((gg[0] + gg[1]) / 2, (gg[2] + gg[3]) / 2); g2 = ((gg[4] + gg[5]) / 2, (gg[6] + gg[7]) / 2)
        fam = Dfam(kc, xc, g1, g2, PA0, PA1, PB0, PB1, ya)
        ex = fam.integral(float(C))
        nint += 1
        slack = v - float(ex)
        worst = min(worst, slack)
        if slack < -1e-12: die(f"leaf {int(r[0])}: C bound {v!r} below exact rational integral {float(ex)!r}")
        nchk += 1
        if nint >= 120: break
    print(f"[3a] exact rational re-integration of {nint} sampled leaves: C bound - exact >= {worst:.3g} OK")
    # ---- 3b. Monte Carlo containment on sampled leaves (independent float volume of the true tube)
    from mc_check import admissible, tube_volume
    rng = random.Random(12345)
    nfeas = ninf = 0; maxratio = -1e9
    for r in rows[:300]:
        st, v = int(r[1]), r[2]
        box = r[3:13]
        for _ in range(4):
            p = [box[2*i] + rng.random() * (box[2*i+1] - box[2*i]) for i in range(5)]
            ok = admissible(p, u, ya)
            if st == 2:
                if ok: die(f"leaf {int(r[0])}: box declared infeasible but contains admissible tube {p}")
                ninf += 1
            elif ok:
                vol = tube_volume(p, u, ya)
                if vol > v + 1e-9: die(f"leaf {int(r[0])}: tube volume {vol} exceeds leaf bound {v}")
                maxratio = max(maxratio, vol - v)
                nfeas += 1
    print(f"[3b] Monte Carlo: {nfeas} admissible sampled tubes all below their leaf bound (max vol - bound = {maxratio:.3g}); "
          f"{ninf} samples in infeasible boxes all non-admissible OK")
    # ---- 3c. root coverage (exact)
    from mc_check import check_roots
    meta["_exe"] = exe
    check_roots(u, ya, grid, meta)
    print("[3c] root grid covers all admissible tubes (exact rational check) OK")
    print(f"CERTIFIED: V_ii(u) <= {maxb:.12f} <= {T}, hence V*(K_M) <= {T} for every M > u = {u!r}")


if __name__ == "__main__":
    main()
