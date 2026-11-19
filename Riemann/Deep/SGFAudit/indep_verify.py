# Independent re-verification (sgfaudit) of the lattice duals CertB, CertB1 at 256 bits.
# Object checked: the integer lists certB/certB1 in Riemann/Deep/ExploreDPP/Lean/TheoremB.lean (eta_k = L[k-1]/10^12),
# i.e. exactly the Lean propositions CertB / CertB1:  (alpha <= theta ? alpha : 0) <= Pfin(h, L, alpha) on [0, 1/(2h)].
# Method differs from verify_dual.py: every cosine/sine evaluated separately from an exact rational argument
# (arb.sin_cos_pi_fmpq), centred third-order Taylor cells with M3 = 2h sum|eta_k|(2 pi k h)^3, adaptive bisection,
# theta a mandatory breakpoint, cells closed. PASS iff every accepted cell has lower bound > 0.
import re, sys, time, random, resource
from fractions import Fraction as Fr
from flint import arb, acb, acb_poly, fmpq, ctx
ctx.prec = 256
src = open('Riemann/Deep/ExploreDPP/Lean/TheoremB.lean').read()
def lean_list(name):
    m = re.search(r'def ' + name + r' : List ℤ := \[(.*?)\]', src, re.S)
    return [int(x) for x in m.group(1).replace('\n', ' ').split(',') if x.strip()]
CASES = {'B': ('certB', Fr(31, 280), Fr(4, 7), Fr(31, 40), 7), 'B1': ('certB1', Fr(119, 600), Fr(1), Fr(119, 200), 3)}
def A(q): return arb(fmpq(q.numerator, q.denominator))
pi = arb.pi()
def setup(case):
    name, h, theta, lam, k0 = CASES[case]
    L = lean_list(name); D = 10**12
    eta = [Fr(x, D) for x in L]
    K = len(eta)
    assert k0 * h == lam
    assert all(eta[k - 1] == -1 for k in range(1, k0)), 'hard core'
    assert all(e >= -1 for e in eta), 'weights'
    assert theta <= Fr(1, 2) / h
    ha = A(h)
    M2 = 2 * ha * sum((abs(A(e)) * (2 * pi * k * ha) ** 2 for k, e in enumerate(eta, 1)), arb(0))
    M3 = 2 * ha * sum((abs(A(e)) * (2 * pi * k * ha) ** 3 for k, e in enumerate(eta, 1)), arb(0))
    return name, h, theta, lam, k0, eta, K, M2, M3
def derivs(h, eta, al):
    # P, P', P'' at rational al; cos(2 pi k h al) = cos(pi * x), x = 2 k h al exact rational
    ha = A(h); s0 = arb(0); s1 = arb(0); s2 = arb(0)
    for k, e in enumerate(eta, 1):
        x = 2 * k * h * al
        s, c = arb.sin_cos_pi_fmpq(fmpq(x.numerator, x.denominator))
        w = 2 * pi * k * ha; ea = A(e)
        s0 += ea * c; s1 -= ea * w * s; s2 -= ea * w * w * c
    return (1 - ha) + 2 * ha * s0, 2 * ha * s1, 2 * ha * s2
def run(case, w0=Fr(1, 500), rmin=Fr(1, 10**9)):
    name, h, theta, lam, k0, eta, K, M2, M3 = setup(case)
    amax = Fr(1, 2) / h
    segs = []
    for (a, b, left) in [(Fr(0), theta, True), (theta, amax, False)]:
        n = -(-(b - a) // w0)
        segs += [(a + (b - a) * i / n, a + (b - a) * (i + 1) / n, left) for i in range(int(n))]
    stack = segs[::-1]; nev = 0; worst = None; wcell = None; gmin = None; gcell = None; ncell = 0; covered = Fr(0)
    t0 = time.time()
    while stack:
        a, b, left = stack.pop()
        c = (a + b) / 2; r = (b - a) / 2
        P, dP, d2P = derivs(h, eta, c); nev += 1
        g = P - A(c) if left else P
        dg = dP - 1 if left else dP
        ra = A(r)
        q = arb(0) if d2P.lower() >= 0 else d2P * ra * ra / 2
        lo = g - abs(dg) * ra + q - M3 * ra ** 3 / 6
        if lo.lower() > 0:
            ncell += 1; covered += (b - a)
            if worst is None or lo.lower() < worst.lower(): worst = lo; wcell = (a, b)
            if gmin is None or g.mid() < gmin.mid(): gmin = g; gcell = c
        else:
            assert r > rmin, f'FAIL near {float(c)}: g={g.str(10)} lo={lo.str(10)}'
            stack += [(c, b, left), (a, c, left)]
    assert covered == amax, 'coverage'
    print(f'{case} ({name}): h={h} theta={theta} lam={lam} k0={k0} K={K} 1/(2h)={amax} M2={M2.str(12)} M3={M3.str(12)}')
    print(f'  cells={ncell} evals={nev} covered=[0,{amax}] exactly; min cell lower bound={worst.lower()} on [{float(wcell[0]):.7f},{float(wcell[1]):.7f}]')
    print(f'  smallest centre value of g = {gmin.str(15)} at alpha={float(gcell):.7f}  ({time.time()-t0:.1f}s, maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB)')
    print('  PASS' if worst.lower() > 0 else '  FAIL')
    return h, eta, M2, theta
def theirs_cell(case, a0, dl):
    # their reported extreme cell, recomputed with their endpoint-Taylor formula but our P, P' (independent evaluation)
    name, h, theta, lam, k0, eta, K, M2, M3 = setup(case)
    P, dP, _ = derivs(h, eta, a0)
    left = a0 + dl <= theta
    g0 = P - A(a0) if left else P; dg = dP - 1 if left else dP
    lo = g0 + dg * A(dl); m = (g0 if g0.lower() < lo.lower() else lo) - M2 * A(dl) ** 2 / 2
    print(f'  {case} extreme cell [{float(a0)}, {float(a0+dl)}]: their-formula bound = {m.str(40)}')
def horner_vs_direct(case, nrand=25, seed=1):
    name, h, theta, lam, k0, eta, K, M2, M3 = setup(case)
    Q = acb_poly([acb(0)] + [acb(A(e)) for e in eta]); ha = A(h)
    rng = random.Random(seed); amax = Fr(1, 2) / h; worst = arb(0)
    pts = [Fr(0), theta, amax] + [Fr(rng.randrange(0, int(amax * 40000) + 1), 40000) for _ in range(nrand)]
    for al in pts:
        if al > amax: continue
        z = acb(0, 2 * pi * A(al) * ha).exp(); P1 = (1 - ha) + 2 * ha * Q(z).real
        P2 = derivs(h, eta, al)[0]
        assert P1.overlaps(P2), f'mismatch at {al}'
        d = abs(P1.mid() - P2.mid()); worst = d if d > worst else worst
    print(f'  {case}: {len(pts)} points (0, theta, 1/(2h), {nrand} random): Horner (acb) and direct sin_cos_pi agree, max |diff of midpoints| = {worst.str(5)}')
if __name__ == '__main__':
    for case in sys.argv[1:]:
        run(case)
        horner_vs_direct(case)
        if case == 'B': theirs_cell('B', Fr(25101, 20000), Fr(1, 20000))
        if case == 'B1': theirs_cell('B1', Fr(78155, 40000), Fr(1, 40000))
