#!/usr/bin/env python3
"""Rigorous certificate (python-flint / Arb ball arithmetic) for the polylogarithmic effective isogeny theorem and
Serre-constant bound (DERIVATION.md).  One command, run from this directory:

  uv run --no-project --with python-flint --with mpmath --with numpy python3 cert.py > logs/cert_output.txt

Reuses the Arb kernel code of ../ChebExp/cert.py (imported, not copied): Kernel (Phi, ||k||, S_j, sup of the
zero ratio by branch and bound), triv_ratio, ram_factor, pow_sum, g_gamma.  What is new here (DERIVATION section 3):
  * the excluded-primes term n_L ||k|| log m / x (Lemma 3.1: a prime p not dividing m);
  * chains in rho = log R, R = rad(2 N_E N_E'), instead of Lambda = log Delta_L: for rho in [ra, rb] the field L
    has n_L <= nmax and Lambda <= Lmax(rb) = (nmax-1) rb + nmax log nmax (Hensel-Serre bound, Lemma 2.3),
    log m = log R <= rb, and x = (a ra + b)^(1/delta);  LHS > 0 is checked at n = 2 and n = nmax (affine in n);
  * a monotone tail for all rho >= rN.
Every chain must start at rho_0 <= log 2 (R >= 2) and end with a tail.
"""
import json, math, sys, time, os, importlib.util
sys.dont_write_bytecode = True                         # leave no __pycache__ in ../ChebExp
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('chebexp_cert', os.path.join(HERE, '..', 'ChebExp', 'cert.py'))
CE = importlib.util.module_from_spec(_spec)            # ChebExp certificate module (Arb), imported under its own name
_spec.loader.exec_module(CE)
from flint import arb, acb, ctx                        # noqa: E402

ctx.prec = 120
up, lo = CE.up, CE.lo
LOG2, LOG3 = CE.LOG2, CE.LOG3
DELTAS = CE.DELTAS

def lmax(nmax, rho):
    """Hensel-Serre: log Delta_L <= (n-1) log rad(d_L) + n log n <= (nmax-1) rho + nmax log nmax (n <= nmax)."""
    n = arb(nmax)
    return (n - 1) * rho + n * n.log()

def lhs_interval_R(delta, nmax, a, b, ra, rb, K, A1, A2):
    """lower bound of LHS_m(x) (DERIVATION (3.3)) for every L with n_L in [2, nmax], log rad(d_L) <= rb,
       log m <= rb, at x = (a ra + b)^(1/delta)."""
    X = a * arb(ra) + b
    lx = X.log() / delta
    assert lo(lx) > up(K.U)                            # x > e^U
    x = lx.exp()
    Mx = up(arb.max(A1, up(((2 * delta - 1) * lx).exp() * A2)))
    xd = up((-delta * lx).exp())
    eps = CE.triv_ratio(K, lx)
    Rf = CE.ram_factor(K, lx)
    pw, mode = CE.pow_sum(K, lx)
    g = CE.g_gamma(K.s0)
    cB = up(1 / K.s0 + 1 / (K.s0 - 1))
    rb_ = arb(rb)
    Lb = lmax(nmax, rb_)
    res = []
    for n in (arb(2), arb(nmax)):
        Bb = Lb / 2 + cB - n * g
        val = K.P1 - xd * Mx * Bb - ((n - 1) * K.P0 + n * K.P0 * eps + 2 * Lb * Rf + n * pw + n * K.kn * rb_) / x
        res.append(lo(val))
    return min(res[0], res[1]), mode, float(lx.mid()), float(Mx.mid()), float(K.P1.mid())

def lhs_tail_R(delta, nmax, a, b, rN, K, A1, A2):
    """lower bound valid for all rho >= rN, x = (a rho + b)^(1/delta) (DERIVATION section 5.3)."""
    rN_ = arb(rN)
    XN = a * rN_ + b
    lx = XN.log() / delta
    assert lo(lx) > up(K.U)
    assert lo(lx) > 2                                  # J <= log x/log 2 + 1 <= 2 log x, and lx x^{-3/4} decreasing
    x = lx.exp()
    Mx = up(arb.max(A1, up(((2 * delta - 1) * lx).exp() * A2)))   # decreasing in x
    eps = CE.triv_ratio(K, lx)                         # decreasing in x
    g = CE.g_gamma(K.s0)
    assert g > 0                                       # -n g <= -2 g
    cB = up(1 / K.s0 + 1 / (K.s0 - 1))
    n = arb(nmax)
    nl = n * n.log()
    al = (n - 1) / 2
    be = nl / 2 + cB - 2 * g
    # zero term: x^{-delta} M B <= M (al rho + be)/(a rho + b); sup over rho >= rN of a Moebius function
    zr = arb.max(al / a, (al * rN_ + be) / XN)
    # 1/x terms: Lmax(rho)/x(rho) and rho/x(rho) decrease for rho >= rN (sufficient conditions asserted)
    assert lo((1 / delta - 1) * a * rN_) > up(b)       # rho/x(rho) decreasing: a rho + b < (a/delta) rho
    assert up(b * (n - 1)) < lo((a / delta) * nl + (1 / delta - 1) * a * (n - 1) * rN_)  # Lmax/x decreasing
    Rfc = arb.min(K.kn * (K.U + LOG2), K.P0 + K.kn * LOG2) / LOG2   # R_f upper bound for every x
    THETA_C = CE.THETA_C
    powx = THETA_C * ((-lx / 2).exp() * (K.S2 + K.Ph2 / 2) + (-2 * lx / 3).exp() * (K.S3 + K.Ph3 / 3)
                      + 2 * lx * (-3 * lx / 4).exp() * (K.kn + K.P0 / 4))
    val = K.P1 - Mx * zr - (n * K.P0 * (1 + eps) + 2 * lmax(nmax, rN_) * Rfc + n * K.kn * rN_) / x - n * up(powx)
    return lo(val), float(lx.mid()), float(Mx.mid()), float(K.P1.mid())

def run_chain(ch):
    name, dexpr, nmax = ch['name'], ch['delta_expr'], ch['nmax']
    a, b = arb(ch['a_str']), arb(ch['b_str'])
    assert float(a.mid()) == float(ch['a']) and float(b.mid()) == float(ch['b'])
    delta = DELTAS[dexpr]
    print(f'\n=== chain {name}: delta = {dexpr}, n_L <= {nmax}, x = ({ch["a_str"]} log R + {ch["b_str"]})^(1/delta), '
          f'log R >= {ch["rho0"]} ===', flush=True)
    assert arb(ch['rho0']) <= LOG2, 'chain must start at or below log 2'
    prev_b = ch['rho0']; ok = True; worst = None; nint = 0
    for k, iv in enumerate(ch['intervals']):
        assert iv['rho_a'] <= prev_b, ('gap', name, k)
        K = CE.Kernel(iv['h'], iv['ss'], iv['s0']).consts()
        A1, A1lo, n1 = K.sup_line(1 - delta, 1e-3)
        A2, A2lo, n2 = K.sup_line(delta, 0.5)
        if iv['kind'] == 'int':
            v, mode, lx, Mx, P1 = lhs_interval_R(delta, nmax, a, b, iv['rho_a'], iv['rho_b'], K, A1, A2)
            prev_b = iv['rho_b']; nint += 1
            tag = f'[{iv["rho_a"]:.6f}, {iv["rho_b"]:.6f}] log x={lx:.3f} pow:{mode}'
        else:
            v, lx, Mx, P1 = lhs_tail_R(delta, nmax, a, b, iv['rho_a'], K, A1, A2)
            prev_b = float('inf')
            tag = f'tail [{iv["rho_a"]:.6f}, oo) log x>={lx:.3f}'
        rel = v / K.P1
        good = v > 0
        print(f'  {tag} h={iv["h"]:.5f} s*={iv["ss"]:.5f} s0={iv["s0"]:.5f} P1={P1:.5g} A1={float(A1.mid()):.5g} '
              f'({n1} boxes) M={Mx:.5g} LHS/P1 >= {float(rel.mid()):.5f} {"OK" if good else "FAIL"}', flush=True)
        ok = ok and good
        worst = rel if worst is None or rel < worst else worst
    ok = ok and prev_b == float('inf')
    print(f'  chain {name}: {"CERTIFIED" if ok else "NOT CERTIFIED"} ({nint} intervals + tail, min LHS/P1 = {float(worst.mid()):.5f})')
    return ok, float(worst.mid())

def ceil4(v):
    """an upper bound for the arb v, printed with 4 decimals and rounded up (float error absorbed by +1e-6)."""
    return '%.4f' % (math.ceil(float(up(v).mid()) * 1e4 + 1e-6) / 1e4)

def table(ch):
    """bounds at sample radicals (upper bounds, rounded up)."""
    a, b = arb(ch['a_str']), arb(ch['b_str'])
    ex = 1 / DELTAS[ch['delta_expr']]
    rows = []
    for lab, R in [('R = 2', arb(2)), ('R = 22 (N_E = 11)', arb(22)), ('R <= 2*10^6 (N_E <= 10^6)', arb(2) * arb(10) ** 6),
                   ('R <= 2*10^100 (N_E <= 10^100)', arb(2) * arb(10) ** 100)]:
        X = a * R.log() + b
        log10p = ex * X.log() / arb(10).log()
        log10C = (ex / 2 * X.log() + LOG2) / arb(10).log()
        rows.append((lab, ceil4(X), ceil4(log10p), ceil4(log10C)))
    return rows

# ---------------------------------------------------------------- algebra checks (DERIVATION section 2, 5)
def _mats(q):
    return [(a, b, c, d) for a in range(q) for b in range(q) for c in range(q) for d in range(q)]

def _mul(x, y, q):
    a, b, c, d = x; e, f, g, h = y
    return ((a * e + b * g) % q, (a * f + b * h) % q, (c * e + d * g) % q, (c * f + d * h) % q)

def _order(x, q):
    one = (1, 0, 0, 1); y = x; k = 1
    while y != one:
        y = _mul(y, x, q); k += 1
    return k

def algebra_checks():
    """brute-force checks of the finite group facts used in DERIVATION sections 2 and 5."""
    out = []
    GL2 = [m for m in _mats(2) if (m[0] * m[3] - m[1] * m[2]) % 2 == 1]
    assert len(GL2) == 6 and len(GL2) ** 2 == 36                       # case (a): [Q(E[2],E'[2]):Q] <= 36
    assert all(((m[0] + m[3]) % 2 == 1) == (_order(m, 2) == 3) for m in GL2)   # trace 1 <=> order 3 (section 2.5)
    out.append('GL_2(F_2): 6 elements, 36 = 6^2; trace 1 exactly on the 2 elements of order 3')
    for k in (1, 2, 3):                                                  # |GL_2(Z/2^k)| = 6*16^(k-1)
        q = 2 ** k
        n = sum(1 for m in _mats(q) if (m[0] * m[3] - m[1] * m[2]) % 2 == 1)
        assert n == 6 * 16 ** (k - 1)
    out.append('|GL_2(Z/2^k)| = 6*16^(k-1) for k = 1, 2, 3 (enumerated)')
    cap = 2 ** 8 - 1                                                     # rank M <= 8 => |G| <= 255
    divs = [2 ** a * 3 ** b for a in range(9) for b in range(3) if 2 ** a * 3 ** b <= cap]
    assert max(divs) == 192 and max(d for d in divs if d % 3) == 128
    out.append('largest 2^a 3^b (b <= 2) <= 255 is 192; largest power of 2 <= 255 is 128')
    M0 = [m for m in _mats(2) if (m[0] + m[3]) % 2 == 0]
    assert len(M0) == 8
    inv = {g: next(h for h in GL2 if _mul(g, h, 2) == (1, 0, 0, 1)) for g in GL2}
    assert all(_mul(_mul(g, a, 2), inv[g], 2) in M0 for g in GL2 for a in M0)   # GL_2 acts on M_2^0
    assert len(M0) * len(GL2) == 48
    out.append('M_2^0(F_2): 8 elements, stable under conjugation; |M_2^0(F_2) x| GL_2(F_2)| = 48')
    C3 = [m for m in GL2 if _order(m, 2) in (1, 3)]
    assert len(C3) == 3
    fixed = [(v0, v1) for v0 in range(2) for v1 in range(2) if (v0, v1) != (0, 0)
             and all(((m[0] * v0 + m[1] * v1) % 2, (m[2] * v0 + m[3] * v1) % 2) == (v0, v1) for m in C3)]
    assert fixed == []                                                   # C_3 image is irreducible (no rational 2-torsion)
    out.append('C_3 < GL_2(F_2) fixes no nonzero vector (irreducible; not absolutely irreducible: eigenvalues in F_4)')
    for ell in (5, 7, 11, 13, 41, 43, 47):
        eps = next(e for e in range(2, ell) if pow(e, (ell - 1) // 2, ell) == ell - 1)   # a non-square
        Cns = [(a, eps * c % ell, c, a) for a in range(ell) for c in range(ell) if (a, c) != (0, 0)]
        J = (1, 0, 0, ell - 1)
        Cplus_minus = [_mul(J, x, ell) for x in Cns]
        assert all((m[0] + m[3]) % ell == 0 for m in Cplus_minus)       # C_ns^+ minus C_ns: trace 0
        gens = {g for g in range(1, ell) if all(pow(g, (ell - 1) // r, ell) != 1 for r in range(2, ell) if (ell - 1) % r == 0 and all(r % s for s in range(2, r)))}
        for m in Cns + Cplus_minus:
            t = (m[0] + m[3]) % ell; d = (m[0] * m[3] - m[1] * m[2]) % ell
            if (1 - t + d) % ell == 0:                                   # eigenvalue 1; other eigenvalue = d
                assert d not in gens, ('ordinary inertia element in C_ns^+', ell, m)
        for m in Cns:                                                    # C_ns: eigenvalue in F_ell => scalar
            t = (m[0] + m[3]) % ell; d = (m[0] * m[3] - m[1] * m[2]) % ell
            if any((x * x - t * x + d) % ell == 0 for x in range(ell)):
                assert m[1] == 0 and m[2] == 0
    out.append('C_ns^+(l), l in {5,7,11,13,41,43,47}: no element with eigenvalues {1, generator of F_l^x}; '
               'C_ns^+ minus C_ns has trace 0; elements of C_ns with an F_l-eigenvalue are scalar (section 5)')
    return out

def main():
    t0 = time.time()
    print('Algebra checks (DERIVATION sections 2, 5):')
    for line in algebra_checks():
        print('  OK', line)
    path = sys.argv[1] if len(sys.argv) > 1 else 'chains.json'
    chains = json.load(open(path))
    allok = True; summary = []
    for ch in chains:
        ok, w = run_chain(ch)
        allok = allok and ok
        summary.append((ch['name'], ch['nmax'], ch['a_str'], ch['b_str'], ok, w, ch['delta_expr']))
        print('  table: R | a log R + b <= | log10 (a log R + b)^(1/delta) <= | log10 (2 (a log R + b)^(1/(2 delta))) <=')
        for r in table(ch):
            print('   ', ' | '.join(r))
    print('\nSUMMARY')
    for s in summary:
        print('  %-6s n_L <= %-4d p <= (%s log R + %s)^(1/delta), delta = %-22s %s  min margin %.5f' % (s[0], s[1], s[2], s[3], s[6], 'CERTIFIED' if s[4] else 'FAILED', s[5]))
    print('ALL CERTIFIED' if allok else 'SOME FAILED', f'({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
