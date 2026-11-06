#!/usr/bin/env python3
"""lincomb3 certificate (DERIVATION.md section 6).

F(s) = sum_{j=1}^N b_j e^{i theta_j} L(s, chi_j): chi_1, ..., chi_N distinct primitive characters mod q >= 3 of the same parity,
b_j nonzero reals, e^{2 i theta_j} = i^a conj(tau(chi_j))/sqrt(q).
Proposition N (DERIVATION section 4, Lamzouri's Pro:BombieriHejhal for N functions) and section 5:
    liminf (zeros of F on the line with odd multiplicity) / N_F  >=  1 - sum_j min(1 - delta_j, 1/N)  >=  sum_j delta_j - (N - 1),
delta_j = liminf proportion of c-isolated zeros of L(s, chi_j) in dyadic blocks, for one fixed c > 0.
delta_j >= p(c) = IsoZero cert.py hinge_p(c) (all constants chi-free, IsoZero DERIVATION 8.2-8.4).
Run: cd Riemann/Deep/LinComb3 && uv run --no-project --with python-flint --with mpmath python3 cert.py   (about 15 s; logs/cert.log)
Every printed lower bound is a float <= the lower endpoint of an Arb ball; every upper bound >= the upper endpoint.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'IsoZero'))
import cert as iso                                   # Riemann/Deep/IsoZero/cert.py (read-only use)
from flint import arb

def main():
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    out = open(os.path.join(HERE, 'logs', 'cert.log'), 'w')
    def P(s):
        print(s); out.write(s + '\n'); out.flush()
    P('== lincomb3 cert (python-flint, prec %d)' % iso.ctx.prec)
    CMT = iso.CMT
    P('C_MT in [%.12f, %.12f]; isolated-zero ceiling of the method 2 - C_MT <= %.9f' % (iso.flo(CMT), iso.fhi(CMT), iso.fhi(2 - CMT)))
    # 1. delta_j = p(c): certified lower bounds (IsoZero hinge bound), several c
    ps = {}
    for c in ('0.001', '0.0001', '0.00001'):
        p, rho, ar, Bv, lam, r0 = iso.hinge_p(c)
        ps[c] = p
        P('c = %-8s: lambda_C = %.3e, rho = %.4f, alpha_rho <= %.6f, B(rho) <= %.4e -> delta >= p = %.9f' % (c, lam, rho, ar, Bv, p))
    c_best = max(ps, key=lambda k: ps[k])
    p = arb(ps[c_best])                               # exact float, <= true p(c)
    P('using c = %s, p >= %.9f' % (c_best, ps[c_best]))
    # 2. the N-function bounds
    P('N | threshold (N-1)/N | naive N p - (N-1) >= | capacity 1 - N min(1-p, 1/N) >= | ceiling N(2-C_MT)-(N-1) <=')
    res = {}
    for N in range(2, 7):
        naive = N * p - (N - 1)
        cap = 1 - N * min(1 - ps[c_best], 1.0 / N) if (1 - ps[c_best]) >= 1.0 / N else None
        capA = (1 - N * (1 - p)) if (1 - ps[c_best]) < 1.0 / N else arb(0)
        ceil = N * (2 - CMT) - (N - 1)
        res[N] = (iso.flo(naive), iso.flo(capA), iso.fhi(ceil))
        P('%d | %.6f | %+.9f | %.9f | %+.9f' % (N, (N - 1) / N, iso.flo(naive), max(iso.flo(capA), 0.0), iso.fhi(ceil)))
    # 3. headline N = 3 and checks
    d3 = res[3][0]
    assert d3 > 0.0174737, d3
    assert res[2][0] > 0.3449825 and res[4][0] < 0 and res[4][2] < 0      # N = 4: even the method's ceiling is negative
    P('HEADLINE N = 3: proportion >= 3p - 2 >= %.9f  (> 0.0174737)' % d3)
    P('N = 2 (reproduces IsoZero cert8 at the same c): 2p - 1 >= %.9f' % res[2][0])
    P('N = 4: naive and ceiling negative; capacity bound = 0 (1 - p >= 1/4): no positive proportion from this method')
    # 4. with Lamzouri's own inputs (p_A = 2 - C_MT - 0.0806 alpha_0) N = 3 fails
    ah, al, it = iso.alpha_tilde(arb(0), arb(0))
    pA = 2 - CMT - iso.D_PAPER * arb(ah)
    P("paper's inputs: p_A in [%.9f, %.9f]; 3 p_A - 2 <= %.6f < 0 (so N = 3 needs p > 2/3: our isolated-zero bound)" % (iso.flo(pA), iso.fhi(pA), iso.fhi(3 * pA - 2)))
    assert iso.fhi(3 * pA - 2) < 0
    # 5. margins
    P('margin over the N = 3 threshold: p - 2/3 >= %.9f' % iso.flo(p - arb(2) / 3))
    P('gap to the method ceiling (N = 3): (4 - 3 C_MT) - (3p - 2) <= %.3e' % iso.fhi((4 - 3 * CMT) - (3 * p - 2)))
    # 6. stage 2 (DERIVATION sections 8-9): delta^R_j >= s0 - eta for every eta > 0, s0 from arXiv 2610.08965 Thm thm:line (exact rationals)
    from fractions import Fraction as Fr
    C0p = Fr(6639737, 5000000)                      # the paper's rational cap C'_0 >= C(f) for its second window (TeX l.1199)
    x2, b2 = Fr(8722, 10**6), Fr(4389, 10**6)        # certificate for the critical line (x, b) = (8722, 4389) 10^-6 (l.1231)
    s0 = (2 - C0p - b2) / (1 - x2)
    assert s0 == Fr(1669159, 2478195), s0           # the paper's s_0
    P("stage 2: C'_0 = %s, (x, b) = (%s, %s); s0 = (2 - C'_0 - b)/(1 - x) = %s = %.9f (exact; equals 2610.08965 thm:line)" % (C0p, x2, b2, s0, float(s0)))
    P('N | 1 - N min(1 - s0, 1/N) (exact) | value')
    st2 = {}
    for N in range(2, 7):
        v = 1 - N * min(1 - s0, Fr(1, N))
        st2[N] = v
        P('%d | %s | %.9f' % (N, v, float(v)))
    assert st2[3] == Fr(51087, 2478195) and st2[2] == Fr(860123, 2478195) and st2[4] == 0
    assert float(st2[3]) > 0.0206146 and float(st2[2]) > 0.3470764
    P('HEADLINE (stage 2) N = 3: >= 3 s0 - 2 = 51087/2478195 = %.9f; N = 2: >= 2 s0 - 1 = 860123/2478195 = %.9f' % (float(st2[3]), float(st2[2])))
    P('stage 2 margin over the N = 3 threshold: s0 - 2/3 = %s = %.6f; 1 - s0 = %.6f > 1/4, so N = 4 gives 0' % (s0 - Fr(2, 3), float(s0 - Fr(2, 3)), float(1 - s0)))
    P('GRH benchmark (Bombieri-Hejhal Thm C): 34.5%%; stage 2 N = 2 exceeds it by %.6f' % (float(st2[2]) - 0.345))
    # 6b. missing lemma (DERIVATION 7.4): positivity needs theta < 1/(N(1 - delta))
    P('missing-lemma thresholds theta < 1/(N(1-delta)): ' + '; '.join('N=%d: %.4f (delta=s0), %.4f (delta=p)' % (N, float(1 / (N * (1 - s0))), 1 / (N * (1 - ps[c_best]))) for N in (4, 5, 10)))
    P('separation constant 2/sqrt(2 pi) = %.4f (DERIVATION 3.2)' % (2 / (2 * 3.141592653589793) ** 0.5))
    # 7. comparison (not certified: quoted) Dousselin, Forum Math. 37 (2025), arXiv 2311.10285 Thm "application"
    dous = {2: 7.38e-9, 3: 4.91e-9, 4: 3.68e-9, 5: 2.94e-9}
    P('before (Dousselin, quoted): N=3 kappa_F >= %.2e; ratio stage 1/old = %.3e, stage 2/old = %.3e' % (dous[3], d3 / dous[3], float(st2[3]) / dous[3]))
    P('Dousselin (quoted): N=2 %.2e, N=4 %.2e, N=5 %.2e, N=10 1.46e-09; large N: 2.161e-06/(N log N)' % (dous[2], dous[4], dous[5]))
    out.close()

if __name__ == '__main__':
    main()
