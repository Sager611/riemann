# BellottiFix NOTE: the flawed steps, corrected inequalities, proofs

Sources: arXiv 2405.12545v2 TeX (final.tex, 1574 lines; v1 has the same passages); code github.com/ChiaraBellotti/
Log-free-zero-density-estimates-for-zeta constants.py (main a43c9e5b = 62c897c0 of 25 May 2024 up to a comment) = code/constants_orig.py;
Rafiq, Zenodo 23131177. Line numbers below are final.tex lines; equation numbers are the compiled ones ((4.11) l.844, (4.13) l.859).

## 1. Flaw (i): counting representatives (l.654 rectangles; l.830 (4.9); l.843-849 (4.11); l.858-871 (4.13), c4)
Setup (l.654): rectangles alpha <= sigma <= 1, (n-1)delta <= t <= n delta, one representative rho_r per rectangle with a zero.
Fix r. Representatives above gamma_r lie in distinct rectangles n_r+1, n_r+2, ..., so the k-th nearest above has
d_k := gamma_s - gamma_r >= (k-1) delta; the same below. The paper claims ONE representative with n delta <= |gamma_r-gamma_s| <=
(n+1) delta for each n (l.843, l.858) and one with |gamma_r-gamma_s| < delta (l.830). The absolute value covers both sides:
with one zero at the centre of every rectangle, every distance n delta (n >= 1) occurs twice.
Per-pair bounds (z = 2 - rho_r - conj rho_s, d = |gamma_r-gamma_s|):
  (P1) X^{Re z}|Gamma(z)| <= 1.0067 X^{2-2alpha} min(1/d, FA/(2 delta))   [|Gamma(z)| <= 1.0067/|z| (4.8); |z| >= max(Re z, d);
       2nu <= Re z <= 2 delta; X^y/y has its minimum at y = 1/log X, so max_{[2nu, 2delta]} X^y/y <= FA X^{2delta}/(2delta),
       FA = e^{a-1}/a if a = 2 nu(T) log X < 1, else 1. The paper's "X^{2-y}/(2-y) decreasing" (l.830) fails when a < 1:
       rows 26-38, FA 1.0011-1.1362 (logs/mono.log; worst case over delta, the step fails for delta with (delta/nu) X^{-2(delta-nu)} > 1):
       a third, small gap in (4.9)]
  (P2) d > 1: X^{Re z}|Gamma(z)| <= X^{2-2alpha} g(d), g(d) = sqrt(2 pi) d^{1.5-2alpha0} e^{-pi d/2 + 1/(6d)}  [Stirling, l.359];
       g decreasing, g(1) < 1 <= 1/1.
Printed (4.9)+(4.11), per r:  sum_{s: d <= 1} X^{Re z}|Gamma(z)| <= 1.0067 X^{2-2alpha} delta^{-1} (1/2 + sum_{n=1}^{1/delta-1} 1/n)
                                                              <= 1.0067 X^{2-2alpha} delta^{-1} (1/2 + 1.12 log(1/nu_i(T1))).   FALSE in general.
CORRECTED (4.9)+(4.11):        sum_{s: d <= 1} X^{Re z}|Gamma(z)| <= 1.0067 X^{2-2alpha} delta^{-1} (1/2 + 2 H_{floor(1/delta)})
                                                              <= 1.0067 X^{2-2alpha} delta^{-1} (1/2 + 2 * 1.12 log(1/nu_i(T1))).
Proof (FA = 1; in general add 2.5 (FA - 1)/delta inside, i.e. c2 -> (5 FA - 4) c2). s = r: 1/(2 delta) by (P1). One side: at most
floor(1/delta)+1 representatives with d <= 1, the k-th with d_k >= (k-1) delta; by (P1) k = 1, 2 give <= FA/(2 delta), k = 3 gives
<= 1/(2 delta), k >= 4 give <= 1/((k-1) delta): total (H_{floor(1/delta)} + FA - 1)/delta per side. Finally
delta > nu_i(T) >= nu_i(T1) (else N(alpha,T) = 0) and H_{floor(1/nu)} <= 1.12 log(1/nu) (checked in all 38 rows; tightest row 1:
5.663 <= 5.694). In the paper's constants: c2 unchanged, c3 -> 2 c3.
Printed (4.13), per r:  sum_{s: d > 1} X^{Re z}|Gamma(z)| <= X^{2-2alpha} sqrt(2 pi) 2.4 I / delta,  I = int_1^inf e^{1/(6y)} y^{1.5-2alpha0} e^{-pi y/2} dy
                        (stated reason: one representative per band; the factor 2.4 is not explained).
CORRECTED (4.13):       sum_{s: d > 1} X^{Re z}|Gamma(z)| <= X^{2-2alpha} sqrt(2 pi) (2 I/delta + 4 e^{1/6 - pi/2}).
Proof. One side, distances beyond 1 sorted d'_1 < d'_2 < ...: the rectangles of d'_1 and d'_i are >= i-1 apart, so d'_i >= 1 + (i-2) delta
(i >= 2). g decreasing: sum <= 2 g(1) + sum_{j>=1} g(1 + j delta) <= 2 g(1) + delta^{-1} int_1^inf g. Two sides: factor 2.
Consequence: 2 + 4 delta e^{1/6-pi/2}/I <= 2 + 4 (1 - alpha0) e^{1/6-pi/2}/I <= 2.120 < 2.4 in every row (F4need in logs), so the
PRINTED c4 IS VALID: its unexplained 2.4 absorbs both sides. The earlier reading "c4 lacks a factor 2" (R160, R162, REFUTATIONS)
is wrong; the missing factor 2 sits in c3 (the dominant term, 82-93% of c2+c3+c4), which has no slack.
Reviewer (A, B): CORRECT.

## 2. Flaw (ii): |G(1/2+it)| (l.933, feeds c5 via (4.15)-(4.21))
G(s) = sum_{j,k<=W} theta_j theta_k [j,k]^{-s}, theta_j = mu(j) log(W/j)/log W (l.388).
Printed (l.933): |G(1/2+i tau)| <= sum theta_j theta_k [j,k]^{-1/2} <= (W^2)^{1/2} |G(1)| <= d5 W/log W.
Both "<=" need theta_j theta_k >= 0 (|sum a| <= sum a and termwise comparison); theta_j changes sign. FALSE as a proof.
CORRECTED: |G(1/2+i tau)| <= sum_{j,k<=W} |theta_j||theta_k| [j,k]^{-1/2} <= c^2 P W/log^2 W <= 9.7572 W/log^2 W,
  P = prod_p (1 + (sqrt p - 1)/p^2) <= 1.47252045 (code/prodP.py), c = sup_{y>=1} S(y)/sqrt y <= 2.574121012 (code/sqfree.py; rounding repaired after review, Arb upper end 2.5741210108618),
  S(y) = sum_{m<=y} mu^2(m) log(y/m)/sqrt m.
Proof. For squarefree j, k: [j,k]^{-1/2} = (jk)^{-1/2} (j,k)^{1/2} = (jk)^{-1/2} sum_{r|(j,k)} h(r), h multiplicative, h(p) = sqrt p - 1 >= 0.
So the sum = sum_r mu^2(r) h(r) Y_r^2, Y_r = sum_{j<=W, r|j} |theta_j| j^{-1/2} <= S(W/r)/(sqrt r log W) <= c sqrt W/(r log W),
and sum_r mu^2(r) h(r)/r^2 <= P. c: exact on every gap between squarefree numbers for y <= 2*10^6 (sup 2.429643); for y >= 2*10^6,
S(y) = int_1^y Q(u)(-f'(u)) du (f = log(y/u)/sqrt u) with Q(u) <= 6u/pi^2 + 3 sqrt u (Q(x) = sum_{d<=sqrt x} mu(d) floor(x/d^2)) gives
S(y) <= (6/pi^2)(4 sqrt y - 4 - log y) + 3(log y + log^2 y/4), ratio <= 2.574121012 (decreasing in y). Crude: c = 4 gives 16P <= 23.5604.
Certified (code/cert.py, arb): P <= 1.472520449, c_tail = 2.574121011 > c_small = 2.42964 (float64: A, B carry relative error
< 1.5e-10, so |error of S/sqrt y| <= (A log y + B)/sqrt y * 3e-10 < 1e-8 against a 0.14 margin), c^2 P <= 9.757067 < 9.7572;
16P <= 23.560328. Reviewer (C): method CORRECT, constants then uncertified -> closed by cert.py; second pass (F): asked for the error propagation
and 9.757067 instead of 9.757066, both now stated (the code uses 9.7572).
Effect: log W in [8.83, 14.90] over the rows and 9.7572/log W < d5 in every row: corrected c5 is 5-40% smaller;
C alone moves by -5.9% ... -0.6%. With the crude 16P, C rises 3.5-17% (logs/fix16.json).

## 3. Falsification of the corrected constants: further slips (independent 50-digit implementation, code/indep.py)
indep.py re-implements every constant from the TeX at mp.dps = 50 (mpmath quad); with the code's deviations switched on it reproduces
his code in all 38 rows to 2.9e-10, and the corrected code to 4.7e-10 (double-precision quad and dps = 10 cost < 1e-9: harmless).
Spot rows 1 (3e12, e^29) and 33 (e^500, e^1000): logs/indep2.json. Deviations found on the way (mode 'all' fixes them all):
(a) (3.10) l.538, (estimatesx) l.605: the residue term carries (u+w)^4; it must be (v+w)^4, since |M(1)| <= d2 log^4(VW), log VW =
    (v+w) log T. Code (S1, cond1) uses (u+w)^4 too. Rows with log T0 <= 35: R1 x ((v+w)/(u+w))^4 ~ 12; row 1 c1 0.9812 -> 0.9226,
    C +13%. Rafiq (his section 2) measured "up to 13%" for 9 rows yet called the term negligible.
(b) code d3 = (TeX d3)/log^2 T0 (Psisquare divides by log t where (specialcase) needs a factor log t): d4 +0.15 ... +1.1%.
(c) code b1 drops sqrt 2 in the r^3 terms of the Corollary (r = sqrt2 delta): b1 +0.0013 (C +0.01%).
(d) fracd2 evaluates d1 at T^(w+x), not T^(v+w) (Rafiq's argument-order slip), and Lemma d2overn applies d1(X) at every u <= X
    (the lower-order terms of Cully-Hugill-Trudgian give 0.746 l^3/3, not /4): d2 +1.7% (row 1); C <= +0.3%.
(e) code cond1 uses e^{1/(6 log T)} for e^{1/(6 gamma)} <= e^{1/(6*0.447 log T)}: c1 -4e-5.
(f) not used or harmless: Patel-Yang's 66.7 (J1-J4, c1) is unproved as published (REFUTATIONS) but true by our 54.6 (RECORDS);
    (tailsul) "Psi(n) = Lambda(n)/log(V/U) for n > U" is false, but |Psi(n)Theta(n)| <= d(n)^2 and the tail carries e^{-log^2 T}.
(g) Rafiq's 7 non-reproduced rows: rows 32, 34-37 are printed with 3 significant figures rounded UP (code 1.11072e6 -> 1.12e6);
    only rows 20-21 (e^50-e^70) differ beyond rounding, printed 0.27% and 0.29% above the code (conservative).
Not re-derived (used as published): Thorner-Zaman Prop. 2.11 (b1, b2; Lemma l.199 prints 4.7098, the Corollary 4.7908: harmless),
Cully-Hugill-Trudgian Thm 2, Ramare's 1.00303 and G1(D), Ramare-Zuniga Alterman 3.09, Graham's lemma, zero-free regions nu_1..nu_4,
the numeric constants 122.092, 1.0486*3.62642, 4.05206, 1.461 (|zeta(1/2+it)| <= 1.4604 for |t| <= 3: true).

