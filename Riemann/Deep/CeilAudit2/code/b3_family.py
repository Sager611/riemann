# B3 float test: the 13/15 floor's count-free generic-row end sums the FULL family (R = 1) in absolute value. Freedom: keep the signs
# (a family first moment). Finite model = BeyondB2/finite_model.py rows (Dirichlet characters mod prime q, exact L-values by FFT), with the
# generic-row summand T = L(1/2 - th, chi) * c(chi), c = 1/L(0.75, eta chibar) (g = 0: no prime product). Diagonal n = m of
# sum_chi chi(n) chibar(m): (#chi) / L(1.25 - th, eta) (eta = chi_7 mod q, or 1). Measure r = |sum T|/sum|T| and r_off after removing it.
import os, sys, math, numpy as np
src = open(os.path.join(os.path.dirname(__file__), '..', '..', 'BeyondB2', 'finite_model.py')).read()
g = {}; exec(compile(src[:src.index("for q in [")], 'fm_defs', 'exec'), g)
hurwitz, proot = g['hurwitz'], g['proot']
def family(q, th, eta):
    gq = proot(q); n = q - 1; ak = np.empty(n, dtype=np.int64); a = 1
    for k in range(n): ak[k] = a; a = a * gq % q
    x = ak / q
    def Lall(s): return q ** (-s) * n * np.fft.ifft(hurwitz(s, x))
    Lw = Lall(0.5 - th); Ls = Lall(0.75); Ld = Lall(1.25 - th)
    j = np.arange(n); ij = (eta - j) % n
    c = 1.0 / Ls[ij]; T = Lw * c
    msk = (j != 0) & (ij != 0)
    diag = msk.sum() / Ld[eta % n] if eta else msk.sum() / Ld[0]
    return T[msk], diag
print("q        th      eta   rows     r_full   mean/diag      r_off   r_off*sqrt(q)")
for q in (12289, 40961, 65537, 163841):
    for th in (0.0, 0.04736):
        for eta in (7, 0):
            T, diag = family(q, th, eta)
            S = T.sum(); A = np.abs(T).sum()
            print("%-8d %.5f %-5s %-8d %.5f  %.5f%+.5fi  %.3e  %.3f" % (q, th, eta if eta else 1, len(T), abs(S) / A, (S / diag).real, (S / diag).imag,
                  abs(S - diag) / A, abs(S - diag) / A * math.sqrt(q))); sys.stdout.flush()
