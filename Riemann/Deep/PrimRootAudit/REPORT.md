AUDIT PASSED: every step of DERIVATION.md, cert.py, tail.py and sieve.py checks out against its sources, independent code reproduces the certified numbers (sieve C_ω exactly, T_ω to every printed digit; κ_ω, c_ω, tail margin, τ-inputs), and exact g*(p) for all 5,761,454 primes 3 ≤ p ≤ 10⁸ stays below every bound (worst log-ratio 0.109); items A, B′, C, D, E hold as written, conditional on H and the cited refereed inputs.

# PrimRootAudit — referee report on discovery prime-primitive-root (Riemann/Deep/PrimRootNum), 19 Nov 2026

Agent ppraudit (Orca). Read-only web (one arXiv e-print, 2608.01498, for the Rosser–Schoenfeld ranges); nobody contacted; no Lean;
PrimRootNum/ and discoveries/ untouched. Chain with D:/C:/T:/EPD: line numbers: FULLPROOF.md. One compute process at a time,
all under memcap 290 MB (peaks 8–163 MB; one unpruned sweep attempt reached 292 MB and was killed, see Files).

## Verdict per item (claim of discoveries/prime-primitive-root.md, verbatim order; H = zero-free Re s > B = 0.8733819666132076)

| Item | Holds? | Decisive evidence |
| --- | --- | --- |
| H from Bebel's Re s > B₂ + 2.2·10⁻²⁵⁴, B₂ < B | yes | cert core: B > B₂ + 10⁻¹⁷ certified; HPAudit/ExplProgAudit (every q, χ, so ζ) |
| A: g*(p) ≤ (κ_ω log p)^{1/(1−B)}, every 1 ≤ ω ≤ 10⁴; κ₁ = 3.53, κ₂ = 11.1, κ₅ = 54.1, κ₁₀ = 216, κ₁₀₀ = 2.72·10⁴, κ₁₀₀₀ = 4.32·10⁶, κ₁₀⁴ = 6.85·10⁸ | yes | own 256-bit recomputation of 14 rows incl. both extremes (logs/indep_kappa.txt): C_ω identical, T_ω equal, L_ω equal, margin (4.4) > 0 at the table κ, own least κ = table κ except κ₂ (own 8.01 < 11.1: conservative); 53-bit margin sweep with the table's T_ω, C_ω over rows 1–2750, every 25th row to 10⁴ and rows 9981–10⁴ (≈ 3,060 rows): no failure, least relative margin 4.4·10⁻⁷ (ω = 2727; 256-bit own-sieve recheck positive) |
| B′: (0.814 log p)^{3/(1−B)} ≤ (…)^{23.70}, p ≥ 5, ω(p−1) ≤ 10⁴ | yes | p ≤ 10⁷: max g* = 211 < 601.44 (own exact table, matches page records); p > 10⁷: κ re-certified at L_m = max(L_ω, log 10⁷); own recomputation for ω ≤ 8 at L_m = log 10⁷ (logs/bprime.txt): max c = 0.813767 (ω = 8, κ = 140) ≤ 0.814, ω ≥ 9 have L_ω > log 10⁷ and table c_ω ≤ 0.79; 23.70 ≥ 3/(1−B) = 23.6933 and 0.814 log p > 1 |
| C: (1.955 log p)^{3.4/(1−B)}, ω(p−1) > 10⁴ | yes | own mpmath (5.6) at c_T = 1.955: margin 4.05·10⁻⁴ (cert 4.08·10⁻⁴), negative at 1.954; true G_q ≥ certified G_lo (max ratio 0.964), true Φ_q − 1 ≤ ψ(a_q) (max 0.765); own I, J, K₃ below the certified upper sums |
| D: (1.956 log p)^{26.86}, every p ≥ 3 | yes | tail 1.955 < 1.956; c^{3.4}L_ω^{2.4} ≥ κ_ω for all 10⁴ rows (cert; every swept row rechecked); 26.86 ≥ 26.8524 with base > 1 |
| E: g* < √p − 2 for p ≥ 10^{97.547}, g* < p for p ≥ 10^{41.032} | yes | binding cases recomputed: ω = 27 at log p = 94.4782 → (1/b)log(κL) = 94.4086 < L; ω = 52 at 224.609 → 112.2502 < 112.3045 = L/2 (− 10⁻³ slack); ω > 120 by the table κ at L_ω (cert; swept rows rechecked), tail by L ≥ 10⁵ |
| Sizes (10^{49.0}, 10^{53.9}, 10^{71.3} at 10¹⁰⁰; 10^{29.2}, 10^{30.2}, 10^{44.4} at 10¹⁰) | yes | recomputed from the κ table |
| "g*(p) ≤ 211 for p ≤ 10⁷" | yes | own segmented exact computation (gstar.py, validated vs two brute-force methods for p ≤ 10⁵): 664,578 primes, max 211 at 4022911, records identical |

## What was checked
- **A. Derivation, line by line (FULLPROOF.md).** Detection (prime powers: ⟨r^k⟩ ⊆ ⟨r⟩), orthogonality on G/G^d, the explicit formula and
  zero sums imported from ExplProg (Lemma EF, P, Z with conductor p ≤ e^{L₉}, L₉ = max(log p, log 9), Zζ, K1; trivial zeros incl. the
  even-χ zero at 0); freezing K1's k's at Y_m is valid because the envelope decreases POINTWISE in Y (Lemma Z is then applied to the frozen
  envelope; no monotonicity of Z_G in k is needed, which would fail since ρ₀ < 0 near 0). Lemma S1 (Rosser sets), Lemma S2 (shifted
  sieve: inclusion-probability identity, h(p_j) = 1/(l_j − 1) ≤ 1/(p_j − 1), σ_m ≤ 0), (2.2), (4.1)–(4.4), Theorems A/B/B′; the tail:
  Buchstab–Selberg weights (T1), Selberg's minimum and Σ|ρ| ≤ c_σξ/G (T2), T3(a)(b)(c), (5.3)–(5.5) re-derived (Abel error terms
  (2Ψ(1) + 4J + K₃)/(2u³) confirmed), R_B(c), the cells and the u ≥ 60 extension, cost K_A/K_B, (5.6) and its extension to L ≥ L₁;
  Theorem D, Corollary E. Every inequality direction checked. Reviewer subagent on S1/S2/(2.2)/T1: CORRECT (remarks: (2.2) needs T ≥ 0
  — certified T > 0 for every weight used; r need not lie in the bump's support, only r ≤ x is used).
- **Rosser–Schoenfeld 1962 ranges** (used as published, only at x ≥ 286 resp. 41): arXiv 2608.01498 (TeX l.686, 806, 826, 1197–1200)
  quotes RS62: |Σ_{p≤x}1/p − log log x − M| ≤ 1/(2log²x) for x ≥ 286 (Thm 5, Cor. p. 70); Π_{p≤x}(1 − 1/p) ≥ e^{−γ}(1 − 1/(2log²x))/log x
  for x ≥ 285 and ≤ e^{−γ}(1 + 1/(2log²x))/log x for x ≥ 1. All uses are inside these ranges (x₀ = 288, Q₂ = 10⁵, z ≥ 104743). R5 (Thm 4,
  θ(x) > x(1 − 1/log x), x ≥ 41) used at z ≥ 104743.
- **B. Code audit.** sieve.py: odd/even positions, Rosser test prod·q³ < y, bulk leaves (prod·q·8 ≥ y ⇒ no grandchild), floor for μ = +1,
  ceil for μ = −1, d = 1 added exactly in cumulative(); cert.py: L_ω as the Arb log of the exact integer P_ω + 1, Y_m > 1 enforced, every κ
  and c certified by pos(), ZoverL's else-branch valid for every L_m, Theorem B′ re-certifies κ at max(L_ω, log 10⁷), small_section's proof is
  the Arb inequality 211 < (0.814 log 5)^{3/b} (the float scan is information only), uniform_c and corollary_thresholds as in the text.
  tail.py: G_lo candidates (exact sums below 288, no primes in [284, 292], Mertens beyond), Π_α and V as products over p < q, ψ's
  branches on exact rationals, upper Riemann sums with Ψ at the right cell end, Rankin majorant on (0, 1/50], restricted integrals with
  t_i ≤ log Q₂/u_r, B_rel at the left end. Coverage: κ table rows ω = 1…10⁴ consecutive; cells start at u₁ = log p₁₀₀₀₁ and run past 60;
  small primes 3 ≤ p ≤ 10⁷ meet the analytic part p > 10⁷; Corollary E scans start at L_ω. No gap, every scan starts at its threshold.
- **C. Independent recomputation** (own code, no import of sieve.py/cert.py; mpmath 256 bits): indep_kappa.py (explicit-stack enumeration of
  D(y) in units 2⁻¹²⁸, exact rationals for ω ≤ 14; kernel constants and Lemma Z integrals by quadrature, sup_{L≥L_m} Z/L by scan): rows
  ω = 1, 2, 3, 5, 8, 10, 27, 37, 52, 100, 613, 1000, 4321, 10⁴ — C_ω identical (11,887,948 at ω = 10⁴), T_ω equal to all printed digits,
  margins positive at the table κ. zcmp.py: own (z_a, z_b, z_e) and sup Z/L equal EP's closed forms to 12–14 digits at the tail point.
  tail_margin.py: (5.6) margin 4.05·10⁻⁴ at c_T = 1.955 with the certificate's rounded τ, K_A, K_B, c_σ (cert 4.08·10⁻⁴; the difference is
  the rounding of those inputs), −2.3·10⁻⁴ at 1.954; own I = 0.882183, J = 3.27051, K₃ = 4.72091 (midpoint, 2·10⁵ cells) below the
  certified 0.882775, 3.27487, 4.72882; B_rel(60) = 0.888613 ≤ 0.889211. tail_check.py: c_σ ∈ [3.2795769, 3.2795802] ≤ 3.2795975;
  h_Q'(2.609) = 1.86·10⁻⁴ > 0 (argmax 2.60950); L₁ = 104403.7613; m_δ = 0.272694; K_B = 0.00906308; Ψ(1) = 1.75049.
  sweep.py (logs/sweep_rows1-2500.txt, logs/sweep_sample.txt): rows 1–2750 and a 25-step sample to 10⁴ — L_ω, margin (4.4) > 0, Theorem D
  (1.956^{3.4}L_ω^{2.4} ≥ κ_ω) and Corollary E (ω > 120 at L_ω) all pass; rows 2727 and 9995 (least sampled margins) re-done with own
  sieve at 256 bits: C 924,930 / 11,886,728 identical, own least κ = table κ (logs/indep_kappa2.txt). sizes.py: the page's and DERIVATION
  sec. 7's size table reproduced digit for digit (logs/sizes.txt).
- **D. Falsification.** gstar.py (segmented, vectorized; validated per prime against two brute-force methods for p ≤ 10⁵): all 5,761,454
  primes 3 ≤ p ≤ 10⁸ in 144 s, 68 MB. Max g* = 223 (p = 73,189,117; new record beyond the page's range), records ≤ 10⁷ identical to the
  page. Worst log g*/log bound: Theorem A 0.10866 (p = 4,022,911, ω = 4), B′ 0.10832 (p = 5), D 0.05881. Per ω(p − 1) = 1…8 the largest
  g* is 3, 61, 131, 211, 193, 223, 193, 137 — bounds are ≥ 10¹⁹. Tail lower bounds tested against truth: exact G_q(ξ_q) (enumeration of
  squarefree d < ξ, P⁺(d) < q) for all 9,577 primes q ∈ [50, 10⁵] at u = u₁, u₁ + 0.1: max G_lo/G = 0.964 (q = 53); true Part A(u₁) =
  0.0222 ≤ certified 0.0390; for the first 300 primes q > 10⁵ at u = u₁, u₁ + 1: max (Φ_q − 1)/ψ(a_q) = 0.765. targ.py (worst cases of the
  shifted sieve): least (probable) primes p = 1 + kP_m, m ≤ 80 (p up to 10^{166.5}), Fermat primes, p = 2q + 1, 4q + 1 near 10²⁰, 10⁴⁰,
  10⁶⁰: worst log g*/log (κ_ω log p)^{1/b} = 0.077 (p = 1 + 4P₉, g* = 149); largest g* 547 (p = 1 + 46P₆₀). No counterexample.


## Findings (no false step; cosmetic only)
1. **Printed c_w column of logs/kappa_table.txt is rounded to nearest (%.6f), not up.** In about half the rows the printed c_w has
   c_w³L_ω² < κ_ω (1,253 of rows 1–2500; max relative shortfall 3.73·10⁻⁶ over the sweep), so the header's "≤ (c_w log p)^{3/(1−B)}" holds for
   the certified value (C:156–158, Fr(cw) checked by pos) but not always for the printed digits. No stated theorem uses the column
   (B: 1.4312, B′: 0.814, D: 1.956, all certified directly). Fix: print ceil at 6 decimals.
2. **κ₂ = 11.1 is valid but not minimal** (own least κ with positive margin 8.01; kappa_fast, C:107–133, returns the first passing
   3-digit value above a non-converged 6-step fixed-point estimate). DERIVATION sec. 4 "cert finds the least κ" is inaccurate for ω = 2
   (all other 13 sampled rows are minimal). Only Theorem B's constant moves (1.4312 → 1.4301 = c₁); the page's items are unaffected.
3. "For every ω ≤ 10⁴" means 1 ≤ ω ≤ 10⁴ (ω = 0 ⇔ p = 2, excluded throughout). h_Q′(a_cap) > 0 is argued in DERIVATION but not
   certified in tail.py (a_cap from a float scan); recomputed here: 1.86·10⁻⁴ > 0, argmax 2.60950 > 2.609.
4. tail.py sets ctx.prec = 120 on import, so cert.py's later sections run at 120 bits (DERIVATION sec. 8 says 200); balls stay rigorous.
5. DERIVATION sec. 0 states R1/R2 with "<"; the published corollary is quoted with "≤" (2608.01498 l.1197); no step needs strictness.

## What survives / dependents
Everything survives as written. Dependents (rg): DISCOVERIES.md l.11 (half-plane consequences line), TREE.md l.15 (N11 record),
RECORDS.md l.15 and Riemann/Deep/HP141/ (θ₁₄₁ rerun: HP141/PrimRootNum/cert.py is identical, HP141/ExplProg/cert.py changes only B and its
core check, so this audit of the method carries over; the θ₁₄₁ numbers 26.85/23.69/7.896 rest on HP141's θ₁₄₁ input and were not rerun
here), Riemann/Tree/R409.lean l.5. COROLLARIES.md, NOGO.md: none (NOGO l.217–218 concern the OpenAI primitive-root preprint).

## Files (Riemann/Deep/PrimRootAudit/)
FULLPROOF.md (chain with line numbers); gstar.py + bf.py (exact g*(p); logs/gstar_1e7.txt, logs/gstar_1e8.txt); indep_kappa.py
(logs/indep_kappa.txt); zcmp.py (Lemma Z vs EP closed forms); tail_check.py (logs/tail_check.txt); tail_margin.py (logs/tail_margin.txt);
sweep.py (logs/sweep_rows1-2500.txt, logs/sweep_sample.txt); sizes.py (logs/sizes.txt); targ.py (logs/targ.txt); bprime.py (logs/bprime.txt). Footprints (memcap): gstar 1e8 68 MB, 144 s; indep_kappa 163 MB; tail_check 23 MB; tail_margin 129 MB; sweep
≤ 20 MB with the cache pruned (one unpruned attempt hit 292 MB from mpmath's per-interval quadrature cache and was killed by memcap; the
full 10⁴-row sweep was then cut to rows 1–2750 plus a 25-step sample because the machine's load (≈ 27) slowed it 5-fold). Rerun: cd Riemann/Deep/PrimRootAudit && uv run --no-project --with numpy python3 gstar.py 1e8; uv run --no-project
--with mpmath python3 indep_kappa.py 1 2 10000; uv run --no-project --with mpmath --with numpy --with python-flint python3 tail_check.py.


RH relevance: N11
