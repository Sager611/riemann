#!/usr/bin/env python3
"""Side result (REPORT.md "Side result", DERIVATION §8): Lamzouri's Thm 1.4 for F = b1 e^{i th1} L(s,chi1) + b2 e^{i th2} L(s,chi2).
Pro:BombieriHejhal + Remark 1.5 (2610.26841 l.251-252, 868-964): proportion of zeros of F on the line with odd multiplicity >= delta1 + delta2 - 1,
where delta_j is a liminf lower bound for the proportion of c-isolated zeros of L(s,chi_j) in dyadic blocks (T, 2T], c > 0 fixed.
delta_j = p(c) of cert.py (hinge bound): every constant in p(c) (C_MT, lambda_C, theta-bar, rho_max, alpha_rho, B(rho)) is independent of chi.
Run: uv run --no-project --with python-flint python3 cert8.py   (about 10 s; output logs/cert8.log)
"""
import os, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cert
from flint import arb

def main():
    os.makedirs("logs", exist_ok=True)
    out = open("logs/cert8.log", "w")
    def P(s):
        print(s); out.write(s + "\n"); out.flush()
    P("== cert8: Lamzouri Thm 1.4 with the hinge bound (python-flint, prec %d)" % cert.ctx.prec)
    # reproduce the paper's delta_0 = 2(0.6552 - 1/2) from its own inputs (c -> 0, D = 0.0806)
    ah, al, it = cert.alpha_tilde(arb(0), arb(0))
    pA = 2 - cert.CMT - cert.D_PAPER * arb(ah)
    P("paper's inputs (c -> 0): p >= %.9f, 2p - 1 >= %.9f (paper: delta_0 = 2(0.6552 - 1/2) = 0.3104)" % (cert.flo(pA), cert.flo(2 * pA - 1)))
    for c in ('0.001', '0.0001'):
        p, rho, ar, Bv, lam, r0 = cert.hinge_p(c)
        pl = arb(p).lower()                         # certified lower bound (float, rounded down)
        d = 2 * pl - 1
        P("c = %-7s: lambda_C = %.3e, rho = %.4f (rho_max %.4f), alpha_rho <= %.6f, B(rho) <= %.4e -> p >= %.9f, 2p - 1 >= %.9f"
          % (c, lam, rho, r0, ar, Bv, p, cert.flo(d)))
    P("GRH-conditional benchmark (Bombieri-Hejhal Thm C): 2(2 - C_MT) - 1 = %.9f" % cert.flo(2 * (2 - cert.CMT) - 1))
    out.close()

if __name__ == "__main__":
    main()
