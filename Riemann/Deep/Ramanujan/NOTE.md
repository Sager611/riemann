# Audit: Zenodo 10.5281/zenodo.23440770 (Padhi, v2, 2026-10-17)

Claim: Ramanujan's inequality pi(x)^2 < (e x/log x) pi(x/e) holds for every real x >= 38,358,837,683, unconditionally
(so N_R = 38,358,837,683). Paper: main.tex of the record; code: github.com/jaideepsaipadhi/ramanujans-inequality (HEAD, 2026-10-22).
Audit 2026-10-22: paper read in full; their Arb pipeline rerun with python-flint 0.9; independent mpmath checks.

## Verdict
REFUTATION of the proof as published. Lemma 6.4 with Lemma 6.2(v) (Section 6.2) undercounts the far off-diagonal
term of the new zero-density bound (Theorem 6.1) by a factor 2. With that fixed, the paper's own certified pipeline gives
q > 1 near log x = 872-882, so Theorem 1.1 (e^58 <= x <= e^3159) and N_R = 38,358,837,683 are not proved.
The gap looks repairable: an unused slack in the same estimate, plus a rigorous replacement for an unproved bound taken
from Bellotti, gives max q = 0.896 and transfer ratio 0.895 on [57, 3159), without Bellotti's published density tables.
Parts we did not re-derive (zero-free-region table, repulsion table, Theorem 6.5 remainders) still carry the whole proof.

## 1. How e^103 < x < e^3158.442 is covered
- Reduction (Lemmas 2.1-2.2, Prop. 2.4): need |pi - li| <= eps_R(L) y/log y, eps_R ~ 1/(3 L^4); from
  |psi(t) - t| <= q eps_R(L2+1)/(1 + 2/L1) t + a(t) on intervals [L1, L2] covering [43, 3159), q < 1 for L2 > 57.
- psi bound: Buthe's smoothed explicit formula, Logan kernel (Lemma 3.1, Props 3.2, 3.5). Smoothing error 0.80319 eps/sqrt(c);
  zeros up to H0 = 3e12 on the line (Platt-Trudgian); zeros with H0 < gamma <= T1 = c/eps bounded by
  (a) zero-free region (1-beta) log gamma >= lambda*(log gamma) (Thm 4.6: BTY inequality evaluated per height, 0.21384 at 29
      instead of BTY's uniform 1/4.896); (b) N+(sigma,T) = min of N1 (zero count), N2 (Kadiri-Lumley-Ng), N3 (Bellotti 2024
      tables), N4 (new Thm 6.1: Bellotti's log-free method with an asymptotically sharp Barban-Vehov mean square d4 ~ 1.6, Thm 6.5).
- Binding range L ~ 875, T1 ~ e^31: needs N+ of order e^4-e^6 per cell for hypothetical zeros with (1-beta) log gamma ~ 0.21-0.32
  just above H0. Only N4 is strong enough; there the zero part is ~50-60% of the error budget, margin 0.9%.

## 2. Requirement eps_R (mpmath, 60 digits; matches the paper and its Arb lower enclosures)
| L | eps_R(L) | 3 L^4 eps_R |
|---|---|---|
| 58 | 3.7428e-8 | 1.27066 |
| 103 | 3.3756e-9 | 1.13977 |
| 877 | 5.7193e-13 | 1.01499 |
| 1000 | 3.3771e-13 | 1.01313 |
| 3158 | 3.3653e-15 | 1.00413 |
Published q at the critical points: 0.400 at L = 103 (T1 < H0: partial RH suffices), 0.99137 on [872,877], 0.8886 on [1000,1010],
0.40001 on [3150,3159] (zero part negligible).

## 3. The false step (Lemma 6.4, Lemma 6.2(v))
Bellotti's far term (2405.12545, (4.13)) sums |Gamma(2 - rho_r - conj rho_s)| X^(2 - beta_r - beta_s) over |gamma_r - gamma_s| > 1
assuming "a unique representative with n delta <= |gamma_r - gamma_s| <= (n+1) delta". Lemma 6.4 fixes the crowding (n(omega) per
window of height omega) but keeps one window per distance band: I' = delta n(omega) sum_i g(1 + i omega). A distance band
[1 + i omega, 1 + (i+1) omega] is two windows (above and below gamma_r), so the bound is 2 n(omega). Correct:
C = 1.0067 (1 + 2 delta L(h)) + 4 sqrt(2 pi) 2.4 I'. Their notes/LEMMA_C.md says "windows ... on each side" but sums once;
certify/arb/cells_arb.py adds 2 sqrt(2 pi) 2.4 Iprime once.
Effect (their code, all 3320 cells recomputed; old values reproduced to 2e-14):
- far share of C 9-23%; log N+ rises by 0.089-0.204 in every cell.
- scan with their kernels: [872,877] 0.99137 -> 1.0726; worst on [797,1000) 1.0872 ([877,882]); [1000,1010] 0.889 -> 0.968.
- re-optimised kernel (c 50-110, f 0.40-0.80, then refined): best q 1.0554 on [872,877], 1.0327 on [874.5,875], 1.0413 on
  [880,880.5] (0.5 is the paper's bisection floor). The cell parameters (u, tau, w, x) cannot help: C is a common factor.
So Prop. 2.4 fails near log y = 875-880 and Theorem 1.1 is not proved.

## 4. Other unjustified steps (minor)
- Lemma 6.2(iv): c5 keeps Bellotti's (4.15) |G(1/2+it)| <= d5 W/log W, "proved" by sum theta_j theta_k [j,k]^(-1/2) <= W G(1):
  invalid (theta_j = mu(j) log(W/j)/log W change sign; absolute values dropped). Rigorous: [j,k]^(-1/2) = (jk)^(-1/2) sum_{r|(j,k)} h(r),
  h(p) = sqrt(p) - 1, Y_r <= 4 sqrt(W)/(r log W), so |G| <= 16 prod_p(1 + (sqrt p - 1)/p^2) W/log^2 W = 23.6 W/log^2 W:
  up to 2.06x the used bound at wL in [15, 24]; c5 (0.012-0.031) stays < 0.07, N+ moves by <= 5%.
- N3 uses Bellotti's published bounds, which carry both flaws plus the (4.11) gap the paper itself flags; "not affected in
  practice" is not proved. Avoidable (section 5 drops N3).
- c1 drops Bellotti's residue M(1) Gamma(1-rho) X^(1-rho) (~ e^(-pi gamma/2), gamma > 3e12): negligible.
- README numbers (worst q 0.99467, transfer 0.99382) are stale; paper and certify/arb logs say 0.99137 and 0.99052.

## 5. Repair (our computation, not in the paper)
In Bellotti (4.6) the residue uses |X^w - Y^w| <= 2 X^Re(w), Y = U/L^2. Since Re w = 2 - beta_r - beta_s >= 2 nu (nu = lambda*(l)/l)
and Y < X: |X^w - Y^w| <= (1 + (Y/X)^(2 nu)) X^Re(w). In the cells (Y/X)^(2 nu) is 0.22-0.41, so C can be multiplied by 0.61-0.70.
With the far term doubled, this factor, and the 23.6 W/log^2 W bound in c5: every cell's log N+ drops by 0.127-0.370 (all still
admissible, c1^2 - c5 >= 0.2). Full rescan with their kernels and N3 removed: max q 0.8957 on [847,857], coverage OK,
Prop. 2.4 transfer ratio 0.8949 at log y = 857: TRANSFER OK. The theorem is plausibly true; the published certificate is not a proof.

## 6. Checked independently (pass)
- eps_R and g(L) < 0 (section 2). Lemma 3.1 constant: max at c = 110 is 0.80314 <= 0.80319.
- Prop. 3.2 = Buthe 2016 (arXiv 1410.7015) Props 2-3, tail 0.16 (x+1)/sinh(c) e^(0.71 sqrt(c eps)) log(3c) log(c/eps);
  |Theta| <= 2 (constant is ~ log 2 pi). Lemmas 2.3, 3.3, 3.4, 3.6, Prop. 3.5: algebra re-derived.
- (B1): max over |t| <= 450 of 1/|zeta(0.6+it)| = 13.13 (t = 14.134) <= 14.57; max |zeta(0.6+it)| = 6.22, top edge <= 0.65, circle |s-1| = 0.1 <= 10.59 (zeta(1.1)), all <= 11.89.
- I_Gamma(-0.47) = 5.92164 <= 5.922; second integral J(-0.44) = 4.71363 <= 4.72.
- Thm 6.5: residue computation of Prop. A.2 re-derived; k1 = H_c(0) = -0.60748, k0 = 0.35203 (own quadrature on Re s = 0.2,
  primes <= 2e5) vs paper -0.607063 +/- 0.0217, 0.35026 +/- 0.0711. d4 in cells (1.56-1.93) matches the main term.
- Lemma 6.3 pair LP, n(omega) in Lemma 6.4, budget (6.1): logic OK. Scan [872,877] q = 0.99137 reproduced.

## 7. Not checked
Thm 4.6 lambda* table (BTY inequality re-evaluated, 4 iterations, re-derived constants E2, E3, G, c10), Thm 5.2 d* table,
Thm 6.5 remainders (REM <= 0.016, S#), Lemma 2.2 tail bound, Hiary-Patel-Yang / Patel-Yang inputs as used.

## 8. Literature
- Johnston 2022 (arXiv 2109.02249) Thm 5.1: 38,358,837,683 <= x <= e^103 (Axler 2018 Thm 3 below e^43). Axler 2022
  (arXiv 2207.02486) Thm 1.1: x >= e^3158.442 via Fiori-Kadiri-Swidinsky. Chain 9658 -> 9032 -> 3915 -> 3604 -> 3361 -> 3158.442
  as in Axler. BTY = Bellotti-Trudgian-Yang arXiv 2603.21490: sigma >= 1 - 1/(4.896 log t), t >= 3. Under RH N_R = 38,358,837,683
  (Dudek-Platt 2015). De's preprints 2407.12052, 2408.02591: the paper's critique looks right (not re-checked).
- Established unconditional status after this audit: 38,358,837,683 <= N_R <= e^3158.442 (unchanged).
- Side finding: Bellotti 2405.12545 (4.13) (one-sided "unique representative") and (4.15) (|G(1/2+it)| <= W G(1)) are flawed;
  her published constants were not re-checked.

## Reproduce
scripts/: cellmod.py (cells: old, far term doubled, repaired), scanmod.py (fixed kernels / kernel search), scanrep.py (full rescan
+ use certify/coverage_check.py and certify/arb/transfer_arb.py), kconst.py (k1, k0), epsr.py (eps_R). They need the upstream repo,
python-flint, mpmath, numpy; set RAMAUDIT to a directory holding repo/ramanujans-inequality-HEAD and the outputs.

