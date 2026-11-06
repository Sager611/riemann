RECORD: J^{d,n} hyperbolic whenever n^2 log(n+2) >= K d^3 with K = 24 (every d >= 1, using the all-n range d <= 9.0e24; d distinct negative zeros when d > 9.0e24), and with K = 2.76e25 using no all-n range (d distinct negative zeros for every d).

# Explicit constant in Savva's two-thirds hyperbolicity wedge (node N6)

## Theorem

Jensen polynomials J^{d,n}(X) = sum_j binom(d,j) gamma(n+j) X^j of Riemann's Xi (xi(1/2+w) = sum gamma(n) w^{2n}/n!).

1. Let D_all >= 9.0e24 be a degree up to which J^{d,n} is hyperbolic for all n (Chasse's sector lemma with RH to height
   3.0001753328e12; our TREE N6 range d0 = 1.614e25). Then **n^2 log(n+2) >= 24 d^3 implies J^{d,n} hyperbolic**, for all
   d >= 1, n >= 0; for d > D_all it has d distinct negative zeros.
2. Self-contained (no all-n range): **n^2 log(n+2) >= 2.76e25 d^3 implies d distinct negative zeros**, for all d >= 1.

Before: Savva (Zenodo 22293642, 3 Sep 2026, unrefereed, AI-assisted) states this wedge with an absolute K that is "effective in
the usual proof-theoretic sense" but prints no number: C_B6 and N_analytic stay symbolic (TeX l. 95-101, App l. 826-831). His
Lean code (not built here) fixes the cutoff n >= ceil(e^(10^15)) + 1, which makes his constant at least e^(3*10^15). Our previous
explicit region was Holland's wedge n^3 log^2(n+2) >= 7000 d^5 (Riemann/Deep/HollandK); the new region contains it for every
d > D_all and is larger by a factor of about 3.9 (d/log n)^{1/6} in the n-threshold.

## Method

Savva's architecture is kept: six matched coefficients (four second-difference log quotients plus two normalizations), the
comparison 3F2 as a multiplicative finite free convolution of two Jacobi polynomials (MMP real roots and log-mesh; MSS max-root
product in both orientations), his four-term recurrence and global-maximum argument for the critical radius, six-node complex
Hermite-Genocchi, and multiplier stability. Replaced: his analytic input (Thm 7.1, sectorial asymptotic with existential R, C)
by our explicit complex saddle bound on |z - n| <= n/10 (HollandK Part A, extended to Q_6), and every unquantified step by Arb
enclosures: the six-match branch by a Brouwer/Krawczyk certificate in each of 1000 cells of eps = 1/L; Gershgorin localization
2.0 and 1.414 instead of 8; C_loc = 4.0006 instead of 12 + 8 sqrt 6; K_r = 2.3 instead of 4096; Cauchy radius 1.33 r
instead of 2 rho, using all six matches (tail e/(lambda^5(lambda-1))). The matching system becomes a signed power-sum problem
whose eps -> 0 limit is the Gauss-Radau rule for 2x dx with node 1 (Savva's y* = (alpha,t,w,delta) = (3, 2, 16/3, 1/3)).
Derivation: DERIVATION.md. Claim-by-claim audit: AUDIT.md.

## Audit (AUDIT.md)

No mathematical error in Savva's argument; we re-derive and use Sections 9-16. Findings: (a) the expanded recurrence coefficient
P_{1,m} (App l. 688) has (m - d) where (m - d + 1) is needed, and the sentence "exact algebra proves equality" with (13.6)
(App l. 692-694) is false as printed; the decomposed form (13.6), which the inequalities use, is correct (exact check 44/44);
(b) the elimination identities (l. 717-722) use 0-based component indices; (c) the sector theorem, the branch contraction and
C_B6 are existential in the paper. Savva does not use the unproved one-sentence complex GORTTW claim that Holland used: he
proves his own sector theorem (and claims a Lean proof); our HollandK bound is an explicit substitute already from n >= 1e12.

## Constants (Savva TeX line -> explicit value)

| Savva | TeX line | here |
|---|---|---|
| K (Thm 1.1) | 81-89 | 24 (with D_all >= 9.0e24); 2.76e25 self-contained |
| N_analytic (sector, contour, branch) | 471-530, 746-760; App 826-831 | N* = Z(40) >= 6.9625e36 (L* = 80); Z(12) = 9.986e11 (L* = 24). Savva's Lean: e^(10^15) |
| R, C of Thm 7.1 | 471-484 | replaced: on abs(z - n) <= n/10, abs(G - G(n)) <= 0.0534, relative saddle error <= 7.9e-18 (L* = 80) |
| derivative tower 2, -2, 4, -12, 48 | 568-578 | Q_2..Q_6 exact; Q_6 = 48 + 26t - 28t^2 + .., t = 2/L |
| contraction ||T_n(y*) - y*|| <= R/2, ||DT_n|| <= 1/2 | 750-753 | Brouwer/Krawczyk per cell; a >= 0.3235, t in [1.9557, 2.0006], w in [5.180, 5.337], delta in [0.3332, 0.3420] |
| C_J = 8 (Jacobi Gershgorin) | 823-830; App 625-657 | 2.00001 (A,B), 1.41441 (C,D) |
| C_loc = 12 + 8 sqrt 6 | 861-866 | 4.0006 |
| K_pre = 256, K_0 = 262144 | 869-870, 1145-1146 | not needed (d/n <= 2.3e-12) |
| P2 >= n/8, abs(P3) <= 2, C_1 < 96, C_0 = 48 | 933-940 | P2 >= 1.9781 n, abs(P3) <= 1, abs(P1) <= 4.0007 n sqrt(Bd), P0 <= 2.0306 n^2 d |
| K_r = 4096 | 943, 955-962 | 2.3 (budget <= 0.97277) |
| Omega radius 2 rho; abs(z - j) <= 3 rho | 985-990, 1039-1045 | radius 1.33 r; abs(z - j) <= d + 1.33 r |
| d + 2 rho <= eta n | 991-998 | (d + 1.33 r)/n <= 6.6e-6 |
| C_B6 | 1004-1011 | n^5 L sup abs(E_F^(6)) <= 2.0504 (leading value -2 = 48 - 10 - 40) |
| C_* in sup abs(E_F) <= C_* d^3/(n^2 log(n+2)) | 1048-1050 | Cfin <= 18.752 (in n^2 L units) |
| sup abs(c_F - 1) <= 1 < 16 | 1058-1061, 1091-1096 | exp(sup abs(E_F)) - 1 < lambda^5(lambda-1) = 1.37332, i.e. sup abs(E_F) <= 0.86329 |
| K_geometry = 5.2e44, K_finite = N0^2(N0+2)+1 | 1124-1134, 1153 | K_N(9e24) = 5.641; self-contained N*^2 log(N*+2) <= 2.756e25 |

## Certificate (one command each)

```
cd Riemann/Deep/SavvaK && uv run --no-project --with python-flint --with sympy python3 certify6.py 80 2.3 1000 1.33
cd Riemann/Deep/SavvaK && uv run --no-project --with python-flint --with sympy python3 certify6.py 24 4 400 1.33 -self
```
First (30 s, certify.log): Kres = 23.033, KN(9e24) = 5.641, "K <= 24: PASS". Second (2 min, certify_self.log): Kres = 830.0,
N*^2 log(N*+2) <= 2.756e25, "KSELF <= 2.76e25: PASS". Both print the documented constants, rounded outward, and write per-cell
witnesses (cells_main.csv: 1000 cells, max Cfin 18.7514, max radius budget 0.97277; cells_self.csv: 2210 cells after bisection).

## Checks (sanity, not part of the proof)

* chk_alg.py, chk_rec.py (audit): factorization coefficientwise; limit system, det = -1/144, inverse row sums (87, 15, 304/3,
  10/3); leading E6 = -2; Jacobi roots within 1.71 sqrt(Vd) of V (Savva's 8 holds); recurrence (13.6) exact 44/44, printed
  expanded P_1 0/44.
* chk_jacobi6.py (chk.log): the closed forms for diag(k) - V and off(k)^2 used in the Gershgorin bound match the Jacobi
  recurrence to 1e-55; the matrix eigenvalues equal the roots of 2F1(-d,U;V;Y/U) (relative 3e-35); Gershgorin radius >= true
  deviation in all 5 cases.
* sanity6.py a (sanity.log): model polynomials at n = 1e8, 1e9, 1e10 (d = 12, 20, 30): all roots simple, positive, within
  2.50 sqrt(Bd) of B (certified 4.0006); max_k abs(T_k)^{1/k} <= 0.86 sqrt(Bd) at every critical point (certified 2.3).
* sanity6.py b: at n = 1e40, d = 7.27e26 on the wedge boundary n^2 log(n+2) = 24 d^3, solving the six matches from the saddle main
  term at 420 digits: branch (0.3246, 1.9598, 5.196, 0.3411) inside the certified ranges; E_F(0..5) = 0 to 1e-376; true
  max abs(E_F) on the boundary of Omega = 0.7998 vs certified 0.8254 (need <= 0.8633); n^5 L E_F(6) = -1.939 (certified <= 2.0504).
* Small Jensen polynomials (HollandK sanity2.py a, re-run today): all 319 J^{d,n}, d = 2..30, n in {0,..,300}, hyperbolic by Arb
  root isolation; this includes the K = 24 wedge-boundary pairs at small n (e.g. n = 300, d = 27).

## Review

Reviewer subagent, round 1 on the argument pasted inline: ERROR, six points. (1) self-contained constant misread as using
L* = 80 (our brief's ambiguity; the L* = 24 run is separate); (2) the all-n range gives hyperbolicity, not distinct zeros:
accepted, statement 1 weakened accordingly; (3) Krawczyk on an approximate map: the argument is now Brouwer existence with a
uniform error bound (|rho| <= 5kL/n on the box, box conditions now checked); (4), (6) Cfin and the radius budget are per-cell
correlated maxima (hull substitution gives 20.6 and 0.979); (5) 1 + lamL is 1.06043, N* exact. Round 2: GAP; points 1, 2, 3, 6
resolved; remaining: (4) the per-cell Cfin witnesses were not in the pasted text (now written to cells_main.csv / cells_self.csv),
(5) from rounded constants K_res <= 23.035 (not 23.034; certificate value 23.0332), new: state a > 0 and 0 <= a eps < 1 on the
cell containing eps = 0 and the continuity of the exact matching map (both now in DERIVATION.md Sec. 3; the code checks
a > 0 and a eps < 1). No mathematical error found in either round.

Orchestrator's verdict-only reviewer (6 Nov 2026, on the chain pasted inline plus DERIVATION.md): CORRECT, conditional on the cited
analytic input and the reproduced Arb enclosures; Sec. 3's uniform perturbation bound plus Brouwer gives exact matching (uniqueness
not needed); Sec. 8 proves simplicity through critical-point signs; the assembly closes (K_res < 23.04 < 24, K_N < 5.65; 2.76e25
exceeds both thresholds); one typo fixed (lambda^5 (lambda - 1) = 1.37332, not 1.37263; theta_E = 0.86329 was computed from the
correct value). SIGNIFICANT: an explicit two-thirds region, beyond the previous explicit bound; exponent and architecture are
Savva's. The orchestrator reran both certificates on a copy: logs and cells_main.csv identical.

## Novelty

* arXiv listings 2026-10/11, math.NT, math.CA, math.CV (.scratch/lst.py 'Jensen|hyperbolic|Tur[aá]n|Laguerre'): no explicit
  Jensen region (2610.38870 is a GRH criterion; 2610.35783 is about Stirling-type polynomials).
* Zenodo '"Jensen polynomials"': Savva 22293642 (the audited paper; no number in the text; Lean cutoff e^(10^15)); Chung
  23987495 and Chen 22951677 are RH claims; Schieber 22936342 is conditional; Dumas 22897068 transfers published bounds;
  Fagiri 23287787 formalizes hyperbolicity of J^{2,0} only.
* HAL '"Jensen polynomials"': Planat hal-05675741 (n >= C0 d^4 with a numerically fitted C0, not rigorous).
* arXiv search (approach arxiv): GORTTW 1910.01227, GORZ 1902.07321, 2007.13582 (Hermite-combination criterion), Farmer
  2008.07206 - none gives an explicit region unbounded in d.
* Our own: HollandK (K = 7000 for the three-fifths wedge, 2.43e114 self-contained) is weaker; nothing else in DISCOVERIES.
* Verdict: an explicit two-thirds wedge constant (24, resp. 2.76e25) is new as far as found; Savva's Lean development implies an
  explicit but astronomical one (>= e^(3e15)).

## Caveats

* Cited, not re-proved: MSS max-root theorem, MMP Props 2.7(iii) and 2.17, Jacobi zeros = Jacobi-matrix eigenvalues, Chasse +
  the verified RH height for D_all, Polya's integral; the HollandK saddle derivation (three reviewer rounds there).
* Statement 1 claims only hyperbolicity for d <= D_all (the all-n range does not give distinct zeros).
* We did not build Savva's Lean development; its claims are reported, not verified.

RH relevance: N6 (RH iff all J^{d,n} hyperbolic; this is the largest explicit hyperbolic region: all n >= (24 d^3/log n)^{1/2})

## Files

certify6.py (certificate), certify.log, certify_self.log, cells_main.csv, cells_self.csv, DERIVATION.md, AUDIT.md, chk_alg.py,
chk_rec.py, chk_jacobi6.py (closed forms of the transported Jacobi matrix; eigenvalues = roots), chk.log, sanity6.py, sanity.log.
