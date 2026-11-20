AUDIT PASSED. (A) jensen-two-thirds-wedge (n² log(n+2) ≥ 24d³; self-contained 2.76·10^25 d³): every step of the assembly re-derived, all certificates reproduced byte-identically, and the exact matching solutions computed independently at n = 10^40 and 10^12 lie in their certified cells; no false step. (B) jensen-wedge-constant (n³ log²(n+2) ≥ 7000d⁵; self-contained 2.43·10^114): the shared explicit saddle bound and every constant where it enters Holland's chain re-derived and recomputed; no false step.

# JWAudit — line-by-line audit of the two Jensen-wedge discoveries (node N6)

Agent jwaudit, 2026-11-20 07:55–08:40 UTC. Details, with source line numbers: [FULLPROOF.md](FULLPROOF.md). Logs: [logs/](logs/).

## 1. Shared input: explicit complex saddle bound (HollandK DERIVATION §2–3, certify.py l. 49–90)

- Re-derived line by line (FULLPROOF §1): the saddle as the fixed point of T_z (invariance, contraction, β ≤ 0.0533474); the contour
  [0, α/2] ∪ vertical ∪ the horizontal ray through the complex saddle (Re(z/u²) > 0 on the ray, concavity); core remainder 16C3/κ0² with
  C3 ≤ c3e^{2α}, C3 s1 ≤ κ0/4; Gaussian, ray and amplitude tails; the pieces near 0 (exponent b·x0, my own derivation, HD states only
  the result); the θ-terms k ≥ 2; the log/Cauchy step for ηG. All bounds are sups over |z − x0| ≤ x0/10 and decrease in a0, so the
  values at a0 = a* hold for every x0 ≥ N*. By hand: rel_core + rel_am = 3.3590·10⁻²⁰ + 1.9624·10⁻²⁰ = 5.3214·10⁻²⁰ = certified δA;
  N* = 45(πe^90 − 9/4) = 1.72530·10^41.
- Q2..Q6 re-derived by a second route (in u) and by numerical differentiation of Ψ (agreement 10⁻⁵⁰, 10⁻⁹⁸; qk_check.py).
- Falsification (saddle_num.py, mpmath 70–200 digits, Mt by quadrature on the ray through the complex saddle, validated against the
  real-axis integral including k ≤ 3 theta terms): no counterexample.

| x0 (edge, 12 points + centre) | true max\|δ\| | claimed δA | true max\|G − G(x0)\| | claimed ηG |
|---|---|---|---|---|
| N*(24) = 9.99·10^11 (SavvaK -self) | 1.9·10⁻¹¹ | 1.16·10⁻⁵ | 0.05048 | 0.05522 |
| N*(80) = 6.96·10^36 (SavvaK main) | 9.1·10⁻³⁶ | 7.90·10⁻¹⁸ | 0.05202 | 0.05337 |
| N*(90) = 1.73·10^41 (HollandK) | 4.1·10⁻⁴⁰ | 5.32·10⁻²⁰ | 0.05209 | 0.05329 |
| 10^45, 10^50, 10^60 | ≤ 7.8·10⁻⁴⁴ | 5.32·10⁻²⁰ | ≤ 0.05228 | 0.05329 |

  True δ ≍ e^{−2a0} (= first-order Laplace term to 4 digits), the bound ≍ e^{−a0}: valid, loose. ηG is the tight constant (2% margin at
  z = 0.9x0) and it holds.
- Reviewer (verdict-only, formulas pasted inline, logs/reviewer.txt): **CORRECT**; no missing term or uniformity loss.

## 2. (A) SavvaK assembly (FULLPROOF §4)

Six-coefficient matching (equivalence with Δ² log γ = Δ² log m, the power-sum form, the 5kL/n error uniform on the box), Brouwer via
Krawczyk with every datum as an interval, positivity and simple roots (both factors Jacobi with positive parameters), sharp Gershgorin
(both bounds re-derived; S + d + 1 = U), MSS localization in both orientations, Savva's recurrence (correct with (m − d + 1); App l. 688's
(m − d) fails, exact check), critical-radius budget, six-node Hermite–Genocchi, the E_F^{(6)} expansion and its error terms, the
stability step (λ⁵(λ − 1) = 1.37332, θE = 0.86329), the signs at critical points, K_res = 23.0332 < 24, K_N(9e24) = 5.641, the regime
(d > D_all forces n > N*, d/n ≤ 2.3·10⁻¹²) and the -self form (K = N*² log(N*+2) = 2.7553·10^25): all check.
Cited results and hypotheses: MSS [PTRF 182 (2022), Thm 1.6] max-root product, nonnegative roots ✓ (positive lower endpoints certified
for the reversal); MMP = Martínez-Finkelshtein–Morales–Perales [IMRN 2024], statements checked in arXiv 2309.10970v3 (Szegő–Walsh (iii),
lmesh preservation), positive roots ✓. Witness files: full coverage of [0, 1/L*], no gaps, maxima as logged (csv_check.py).

## 3. (B) HollandK constants (FULLPROOF §5)

C3 s1 = 2.2901 ≤ κ0/4-bound 3.1237 ✓; δA, ηG as §1 ✓; a + b + ϑ = 0.97084 < 1 suffices (Holland l. 1797–1807: at the maximizer
M ≤ (a + b + ϑ)M) ✓; a_m, b_m, ϑ and an explicit C_V for Holland's C* re-derived (U-backward induction from U_{d+1} = 0) ✓;
K_res = (Cfin/θE)²(1 + lamL)² = 6416.59 certified (6416.74 from the rounded displays) ≤ 7000 ✓; K_N(9e24) = 784.13 ✓; K_SELF = 2.4291·10^114 ✓.
Outward rounding: fup/flo round the exact rational ball endpoints with ROUND_CEILING/FLOOR plus asserts; every quoted constant matches.

## 4. Reruns and independent spot check (FULLPROOF §3, §6)

- HollandK certify.py (L* = 90 and 80) and SavvaK certify6.py (main and -self) from copies: all four logs and both witness CSVs
  byte-identical (7–134 s, ≤ 85 MB). SavvaK chk_*.py, HollandK check_*.py: identical / all True.
- Spot check, different method (spotcheck.py): n0 = 10^40 (cell 917), h(n0 + j) by 330-digit quadrature, the exact matching equations
  solved by Newton (no Q_k, no R_k series, no error ball): (a, t, w, δ) = (0.324612, 1.959802, 5.196247, 0.341067), inside the
  certified box [0.324233, 0.324993] × [1.959145, 1.960466] × [5.192006, 5.200518] × [0.340898, 0.341235]; X0..X3 inside their
  enclosures; n^5 L E_F^{(6)} = −1.939 against C6 ≤ 1.998. Same test in the weakest regime (-self, L* = 24, n0 = 10^12, last cell):
  (0.303347, 1.869791, 4.883370, 0.362326) inside [0.302672, 0.304027] × [1.868638, 1.870967] × [4.875773, 4.891046] × [0.361976, 0.362670].

## 5. All-n input d ≤ 9.0·10^24 (FULLPROOF §7)

Follows from GORTTW (Adv. Math. 397 (2022) 108186) Theorem 1.2 with m = 0 (TeX l. 312–315; m = 0 is Chasse's Thm 1.8): RH_0(T) and
d ≤ ⌊T⌋² ⇒ J^{d,n} hyperbolic for all n; proof re-read (sector S(π, 2 arctan(1/2T)), d ≤ T² + ½ + 1/(16T²)). Platt–Trudgian
T = 3,000,175,332,800: ⌊T⌋² = 9.00105·10^24 ≥ 9.0·10^24 ✓.

## 6. Findings (no effect on any claim)

- Discovery page (A) and the brief name "Mirabelli–Marcus–Peters"; the source is Martínez-Finkelshtein–Morales–Perales (IMRN 2024).
- Discovery page (B): "Kres < 6416.68 even from rounded displays" — the rounded displays give 6416.74 (certified 6416.59).
- HD l. 85 quotes the |Q5| disc as |t| ≤ 1/(a* − 0.07), the code uses 1/(alo − 0.01); reported value only.
- GORTTW l. 752 (σ ∈ S ⊂ C^×) should admit σ = 0; l. 775's positivity sentence is superfluous for m = 0.

## 7. Not checked

Holland's Laguerre-defect, Jacobi-model and four-coefficient-matching blocks (HD §4–5; outside the brief, orchestrator reviewer CORRECT);
the cited theorems themselves (MSS, MMP, Chasse/GORTTW, Jacobi matrix, Pólya's integral); polynomial prefactors of the
super-exponentially small terms (irrelevant against exponents ≤ −10^10); Savva's Lean development; the alternative D_all = 1.6141·10^25.

## Files

FULLPROOF.md; saddle_num.py, qk_check.py, spotcheck.py, cert_one_patch.py, csv_check.py; logs/ (reruns.txt, saddle_num_4pts.log,
saddle_num_12pts.log/.out, qk_check.log, spotcheck_1e40.log, spotcheck_1e12.log, cert_one_cell917.log, cert_one_self_cell399.log,
csv_check.log, savvak_chk_rerun.log, reviewer.txt). Every run under
memcap 290 MB (peaks ≤ 85 MB).

RH relevance: N6
