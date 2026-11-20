# P2 (alpkato): does a PUBLISHED bound imply LemmaK? kBound(C,e,n,D) = exp((2en)^C log 2D). Work with L = log D.
# (i) Bayer-Fluckiger-Emery-Houriet 2014 (arXiv:1406.5953, Thm 1.1), K-degree m = 2e, degree d = n, ell = max(d+1, 2m+2):
#     log card_ell K_m(O_F)_tors <= exp(F),  F = 71 m^4 d^3 log(2m+1) + 293 m^5 d^5 log d + 528 m^5 d^4 log D.
#     small-torsion loss: #A <= card_ell(A) * ell^rank(A), rank <= log2 h_L + n, h_L <= sqrt(D)(1+log sqrt D)^(n-1).
#     claim: F + log(1+X) <= (2en)^C log(2D), X = (log2 h + n) log ell; then exp(F) + X <= exp(F)(1+X) <= kBound.
# (ii) Soule 2003 Prop 4 i): loglog k(m) <= 2^20 m^4 log(m) n^(4n) sqrt(D) log(D)^(n-1): where it stops implying kBound(24).
import math
def F(e, n, L):
    m = 2*e; d = n
    return 71*m**4*d**3*math.log(2*m+1) + 293*m**5*d**5*math.log(d) + 528*m**5*d**4*L
def X(e, n, L):
    ell = max(n+1, 4*e+2)
    log2h = (0.5*L + (n-1)*math.log(1+0.5*L))/math.log(2)
    return (log2h + n)*math.log(ell)
grid = [(e, n, L) for e in range(1, 41) for n in range(2, 81, 2) for L in [math.log(3), 2, 5, 10, 30, 100, 1e3, 1e4, 1e6]]
def lhs(e, n, L): return F(e, n, L) + math.log1p(X(e, n, L))
def rhs(C, e, n, L): return (2*e*n)**C*(math.log(2) + L)
worst = min((rhs(16, *g)/lhs(*g), g) for g in grid)
print("BEH => LemmaK(C=16): min RHS/LHS on grid = %.3g at (e,n,logD) = %s" % worst)
for C in range(8, 17):
    print("  C = %2d holds on grid: %s" % (C, all(rhs(C, *g) >= lhs(*g) for g in grid)))
print("Analytic: u = 2en >= 4, m = 2e, d = n: m^4 d^3 <= u^4, m^5 d^5 = u^5, m^5 d^4 <= u^5, log(2m+1) <= u, log d <= u, log D <= log 2D, log 2D > 1:")
print("  F <= 71 u^5 + 293 u^6 + 528 u^5 log 2D <= 892 u^6 log 2D;  log(1+X) <= u log 2D;  so F + log(1+X) <= 893 u^6 log 2D <= u^11 log 2D (4^5 = 1024 > 893).")
print("  grid check of F <= 892 u^6 log 2D:", all(F(e, n, L) <= 892*(2*e*n)**6*(math.log(2)+L) for (e, n, L) in grid),
      "; of log(1+X) <= u log 2D:", all(math.log1p(X(e, n, L)) <= (2*e*n)*(math.log(2)+L) for (e, n, L) in grid))
print("Soule Prop.4(i) vs kBound(24), e = 1 (K_2): smallest log D where Soule's loglog bound exceeds kBound(24)'s loglog")
for n in [2, 4, 6, 8, 12]:
    m = 2
    def sl(L): return math.log(2**20) + 4*math.log(m) + math.log(math.log(m)) + 4*n*math.log(n) + 0.5*L + (n-1)*math.log(max(L, 1e-9))
    def kl(L): return 24*math.log(2*n) + math.log(math.log(2) + L)
    lo = next((i*0.05 for i in range(1, 40000) if sl(i*0.05) > kl(i*0.05)), None)
    print("  n = %2d: Soule's bound is weaker than kBound(24) for log D > %s" % (n, lo))
print("  (n = 2: Soule Prop.4(iii) instead has the kBound shape: loglog k(m) <= 1120 m^4 log(m) log D.)")
