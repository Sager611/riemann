# Audit of Savva, "A Two-Thirds Hyperbolicity Wedge for Jensen Polynomials of Riemann's Xi-Function"

Source: J. Savva, Zenodo 22293642, version 1.1, 3 Sep 2026 (unrefereed, AI-assisted; files JENSEN_TWO_THIRDS_MAIN.pdf,
TECHNICAL_SUPPLEMENT.pdf). Line numbers "l." refer to the TeX source of the same release (GitHub aimathpapers/Jensen-2-Thirds,
tag v1.1.0, PAPER_SOURCE/JENSEN_TWO_THIRDS_MAIN.tex); "App l." to PAPER_SOURCE/c48_detailed_appendices.tex; "p." to the PDF.
Checks: chk_alg.py / chk_rec.py (exact sympy, run during the audit, outputs quoted below), certify6.py, sanity6.py.

**Verdict.** No mathematical error found in the architecture. Two typos (one false "exact algebra" sentence). The paper's
analytic constants (C_B6, N_analytic, and the branch contraction threshold) are existential in the text; the accompanying Lean
code (not built here) fixes the cutoff n >= ceil(e^(10^15)) + 1, i.e. an explicit but astronomical K >= e^(3*10^15). Savva does
NOT use the unproved complex GORTTW sentence that Holland used: he proves his own sector theorem (Thm 7.1). Our certificate
replaces his analytic input by our explicit complex saddle bound (Riemann/Deep/HollandK) and all of his constants by
certified ones, giving K = 24 (with the all-n range d <= 9.0e24) and K = 2.76e25 (self-contained).

| # | Savva's claim (section, p., TeX line) | Status | Evidence / note |
|---|---|---|---|
| 1 | Thm 1.1 (p. 4, l. 81-89): n^2 log(n+2) >= K d^3 => d distinct negative zeros, absolute K | **true**, now with explicit K | certify6.py: K = 24 (d > D_all = 9.0e24 via the wedge; d <= D_all by the all-n range), K = 2.76e25 self-contained |
| 2 | "effective ... two analytic constants can be read from the sectorial saddle proof" (l. 95-101; App l. 813-837; Sec. 17 l. 1137-1169) | true in principle; no number in the paper | paper: K_geometry = 5.2e44 (l. 1153) and K_final = max(K_geometry, C_mult C_B6, K_finite) with C_B6, N_analytic symbolic (App l. 826-831). Lean: cutoff ceil(exp(10^15))+1 (XiNaturalExplicitBranch.lean l. 25-31), global constant max(endpoint, N0^2(N0+2)+1) (XiNaturalMultiplierCertificate.lean l. 391-393) => K >= e^(3e15) |
| 3 | Def. 2.1, normalization (p. 6, l. 178-196) | true | standard; our gamma = C 4^{cz} Gamma(z+1) Mt(z)/Gamma(2z+1) agrees up to constant and linear terms, which no second difference sees |
| 4 | Prop. 4.1 factor eight (p. 8, l. 276-302) | not used | our normalization is independent (HollandK) |
| 5 | Lemma 5.1 sector saddle branch (p. 9, l. 323-381) | not used | Rouche with explicit numbers; plausible |
| 6 | Sec. 6 contour (p. 10, l. 384-457): R_4 <= C|K_s||r|^4, tails O(e^{-c|K_s|^{1/5}}), "all sufficiently large |s|" (l. 419) | not used; constants unquantified in the paper | supplement Sec. 6 and Lean modules Leading*.lean claim explicit witnesses; not built here |
| 7 | Thm 7.1 sector asymptotic, |E(M)| <= C log|M|/|M| on |arg M| <= 1/400 (p. 11, l. 471-530) | not used; R, C existential in the paper | replaced by our complex saddle bound on |z - n| <= n/10 (HollandK Part A), valid for n >= 1e12 (L* = 24) and n >= 7e36 (L* = 80) |
| 8 | Derivative tower h^(k) = (2, -2, 4, -12, 48)/(x^{k-1} L)(1 + O(1/L)) (p. 12, l. 568-578) | true | our Q_2..Q_6 (exact, sympy) have leading terms 2, -2, 4, -12, 48; Q_6 = 48 + 26t - 28t^2 + ..., t = 2/L |
| 9 | H_6 majorant < 10^4, |G_0^(6)| <= 20000/(|N|^5 log|N|) (l. 588-600) | not checked; not needed | we enclose Q_6 directly in Arb |
| 10 | Lemma 9.1, four quotients + two normalizations (p. 13, l. 615-627) | true | trivial induction |
| 11 | Cube calculus, limits H^inf (l. 689-697), S^inf = (-2,-1,-2,-2) (l. 700-702) | true | rederived as the power-sum system of DERIVATION.md Sec. 3 (Gauss-Radau rule for 2x dx with node 1) |
| 12 | Limit root y* = (3, 2, 16/3, 1/3), DF(y*), det = -1/144, ||DF^{-1}||_inf = 304/3 (p. 15, l. 709-742) | true | chk_alg.py: F(y*) = 0, det = -1/144, row sums (87, 15, 304/3, 10/3) |
| 13 | Elimination identities "3F_1 - F_2", "(4/3)F_2 - F_3" (l. 717-722) | true with 0-based component indices | with 1-based indices (as "Writing F_j for the four components" suggests) the identity fails; notation only |
| 14 | Contraction ||T_n(y*) - y*|| <= R/2, ||DT_n|| <= 1/2 "after increasing the analytic threshold" (l. 746-760) | true in principle; unquantified | replaced by a per-cell Krawczyk certificate (unique solution in a box, every n >= N*); certified ranges a >= 0.3235, t in [1.9557, 2.0006], w in [5.180, 5.337], delta in [0.3332, 0.3420] |
| 15 | Factorization p_F = Q^(1)_{A,B} [x] Q^(D)_{C,D} (p. 16, l. 784-790) | true | chk_alg.py: coefficientwise, 5 random rational cases |
| 16 | MMP Prop. 2.7(iii), 2.17 give d simple positive roots (l. 803-821) | true (cited) | same citations as Holland; hypotheses A >= B + d, C >= D + d certified |
| 17 | Jacobi localization constant 8 (l. 823-830; App l. 625-657) | true, generous | numerics: max|root - V|/sqrt(Vd) <= 1.71 (chk_alg.py), always <= 2 sqrt((U-V)/U); our Gershgorin: 2.00001 for (A,B), 1.4144 for (C,D) |
| 18 | C_loc = 12 + 8 sqrt 6 < 32 via MSS in both orientations (l. 835-871) | true | arithmetic checked; sharp version C_loc <= 4.0006; numerics <= 2.50 (sanity6.py a) |
| 19 | ODE / four-term recurrence, decomposed coefficients (13.6) (p. 18, l. 871-926) | true | chk_rec.py: recurrence holds in 44/44 exact cases with (13.6) |
| 20 | Expanded coefficients (App l. 676-690) and "exact algebra proves equality of (app:P-expanded) with (13.6)" (App l. 692-694) | **typo; the sentence is false as printed** | P_{1,m} must read ... + (m - d **+ 1**)(A + C + 2m + 1); printed form fails the recurrence 0/44; difference = -Dy(A+C+2m+1)/(AC). Not load-bearing (inequalities use (13.6)) |
| 21 | Recurrence bounds P2 >= n/8, |P3| <= 2, |P1| <= C1(n sqrt(Bd) + nd), C1 < 96, 0 <= P0 <= 48 n^2 d (p. 19, l. 933-940) | true | true sizes: P2 ~ 2n, |P1| ~ 4 n sqrt(Bd), P0 ~ 2 n^2 d (certified: >= 1.978 n, <= 4.0007, <= 2.031) |
| 22 | Global-maximum argument, Prop. 14.1 with K_r = 4096 (l. 943-969; App l. 696-720) | true | logic verified; reused with K_r = 2.3 (budget 0.9728 < 1); numerics: max_k |T_k|^{1/k}/sqrt(Bd) <= 0.86 |
| 23 | E_F^(6) formula (p. 20, l. 975-983), domain containment (l. 991-998) | true | rederived; certified R0/n <= 6.6e-6 |
| 24 | sup |E_F^(6)| <= C_B6/(n^5 log(n+2)), C_B6 "kept symbolic" (l. 1004-1011) | true; now explicit | n^5 L sup|E_F^(6)| <= 2.0504 (leading value -2 = 48 - 10 - 40); sanity: -1.939 at n = 1e40 |
| 25 | Six-node complex Hermite-Genocchi (l. 1014-1037) | true | standard (simplex mass 1/720) |
| 26 | |z - j| <= 3 rho and (eq:exponent) (l. 1039-1052) | true | we use the sharper (d + lambda r)^6 |
| 27 | Multiplier stability, Lemma 16.1 (p. 21, l. 1064-1106) | true | divided differences via real Hermite-Genocchi + Cauchy on circles centred on [0,d]; we use radius lambda r and six matches: tail e/(lambda^5(lambda-1)) |
| 28 | Assembly, P(y) = J^{d,n}(-y/S)/gamma(n), S = B R_1 (l. 1108-1134) | true | checked by coefficient comparison |
| 29 | Finite absorption K_finite = N_0^2(N_0+2)+1 (l. 1124-1134; App l. 786-811) | true | elementary |
| 30 | Lean: whole chain kernel-checked modulo typed Jacobi/MMP/MSS inputs (Sec. 18, l. 1171-1251) | not verified (no Lean build in this task) | source read only: no sorry/axiom in JensenWedge/*.lean (rg); the literature inputs are hypotheses of the headline theorem |

Uses the unproved complex GORTTW sentence? **No.** Holland's analytic input is the one-sentence complex form of GORTTW (3.2);
Savva instead proves his own sector theorem (Secs. 5-7, supplement Secs. 3-7) and claims a Lean proof of it with explicit
witnesses. In the paper the constants R, C, c of that theorem are existential. Our HollandK bound is an explicit substitute that
already holds from n >= 1e12.
