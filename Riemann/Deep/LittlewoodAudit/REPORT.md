AUDIT PASSED: every item (a)–(f) of littlewood-half-plane holds as written, apart from three cosmetic slips (Dusart theorem numbers off by one, two limit constants rounded up in the 6th digit, one Davenport chapter reference), none of which changes a bound.

# LittlewoodAudit: line-by-line referee of littlewood-half-plane (Riemann/Deep/Littlewood78), 19 Nov 2026

Agent lwaudit (Orca), node N11. Read-only web access (arXiv TeX of Dusart 1002.0442). Nobody contacted. No Lean, no commits.
Riemann/Deep/Littlewood78/ and discoveries/ untouched. The chain with source line numbers is in FULLPROOF.md.

## Verdict per item (hypothesis H(B), B = 0.87338196661320752, Lean-checked: HPAudit l.1, l.38)

| Item | Claim | Holds? | Evidence |
| --- | --- | --- | --- |
| (a) | \|L(1, χ)\| ≥ 0.11/max(3, log log q), every primitive χ mod q ≥ 3 | yes | proof refereed; independent adaptive cover of [3, 40] (152 cells, worst margin 0.0178), Λ ≤ 3 value F(e³) = −3.287948, tail c_40 = 0.115659 |
| (b) | 0.115/log log q for log log q ≥ 3.13; 0.116 for ≥ 4 | yes | independent covers from exactly 3.13 (worst margin 0.00094) and 4 (margin 0.119 at 4; binding end 99.7–100, 0.0203); tail c_100 = 0.116424 |
| (c) | liminf \|L(1, χ)\| log log q ≥ (π²/6)e^{−γ}(1 − B) | yes | Lemma T + limit h = 1/Λ_N rechecked; constant = 0.1169398360… (page's "= 0.11694" is rounded up, slip S2) |
| (d) | h(−d) ≥ 0.0350√d/max(3, log log d), d > 4; > 2000 for d ≥ 2.957·10¹⁰; > 10⁵ for d ≥ 7.976·10¹³ | yes | h = √d L(1)/π (w = 2), 0.11/π = 0.035014; all five thresholds (100 … 10⁵) re-certified by independent cells starting below log d₀ |
| (e) | \|L(1, χ)\| ≤ 15 log log q for log log q ≥ 3 | yes | Prop. 11 and tail T' refereed; independent cover of [3, 40] (worst margin 0.0192 at [3, 3.001]); tail 14.2223 |
| (f) | B = 7/8: 0.105/max(3, log log q), limit 0.11545 | yes | independent rerun with B = 7/8 (lwcheck_78.log): worst margin 0.0526, tail 0.114185; limit = 0.1154455 (rounded up, S2); 7/8 thresholds 3.047·10¹⁰, 8.137·10¹³ re-certified |

## Slips (wording only; no bound changes)

- S1 Dusart numbering. In arXiv 1002.0442 (only version v1, 2 Feb 2010) section 6 has a numbered Remark (TeX l.368) and a
  commented-out Corollary (l.529), so Σ_{p≤x} 1/p is Theorem 6.10 (l.660–670) and Σ_{p≤x} log p/p is Theorem 6.11
  (l.692–702); 6.12 is Mertens' product. PROOF.md (l.153, 154, 252), Littlewood78/REPORT.md (l.101, 102, 131) and
  dusart_check.py cite 6.11 and 6.12. The quoted statements, constants and ranges (x ≥ 10372, x ≥ 2974) are exact.
  PROOF.md l.252 also cites "Ramanujan J. 45 (2018)"; only the arXiv v1 statements were checked, not the journal version.
- S2 Rounding direction of limits. (π²/6)e^{−γ}(1 − B) = 0.11693984 and (π²/6)e^{−γ}/8 = 0.11544548; the page (line 3),
  Littlewood78/REPORT.md (l.1, l.132) and PROOF.md l.30 print 0.11694 and 0.11545, which overstate them by 1.6·10⁻⁷ and
  4.5·10⁻⁶. Read as decimal lower bounds they are not proved; write 0.116939 and 0.115445. PROOF.md l.18 and l.233 are right.
- S3 Davenport chapters. PROOF.md l.58 cites Ch. 16–17 for the bounds on L'/L used to move the contour; Ch. 16 (N(T, χ))
  has the strip bound, but Ch. 17 is the ζ case; the χ version for Re s ≤ −1 is Ch. 19. Qualitative input, standard. Also
  l.63 says H(B) gives L(σ, χ) ≠ 0 for σ ≥ 1, which is unconditional.

## Proof referee (details in FULLPROOF.md, sections 0–8)

- Explicit formula (Lemma 2): residues −L'/L(σ) at w = 0 and −m_ρV(ρ − σ) at the zeros; signs right; trivial zeros
  (s = 0 for even χ) included. log L(1, χ) = −∫_1^∞ L'/L gives P_0 + Σ_ρ G(ρ); interchange justified by Lemmas 5, 7.
- Hadamard (Lemma 4): Davenport Ch. 12 (17)–(18) at s = 1, real parts; ½ψ((1 + a)/2) ≤ −γ_E/2 for both parities;
  primitive χ, so the conductor is q. Pairing ρ ↔ 1 − ρ̄: the functional equation gives a zero 1 − ρ of L(s, χ̄),
  conjugation brings it back to L(s, χ) at the same height, multiplicity kept; so 1 − B ≤ β ≤ B for complex χ too.
- Lemma 6 (pairs): redone (linear in r, two endpoints, fixed points β = ½ covered). Loop: κ_1 < 1 is a hypothesis of
  Prop. 10, enforced in cert.py by the certain ball test k1 < 1 and on the tail by Lemma T(i) (κ_1 decreasing in L).
- Prime sums (Lemmas 8, 9): every inequality checked, including Re w/(1 − w) ≥ −r/(1 + r), the τ̄ majorant and its
  monotonicity ranges, Lemma 9(c) = (M − γ_E) + log ζ(2) = 0.181982, and Dusart's ranges (e^10 > 10372, e^8 > 2974).
- Lemma T (the repaired tail): (i) κ_1 decreasing; (ii) ω(8)|_{h=2} ≤ −1.7237; (iii) κ_0L, κ_0e^Λ, 1/(1 − κ_1), ε_D,
  1/(x − 1), τ̄, T_0 all decrease in Λ (b < ½ is what makes κ_0e^Λ decrease). Limit argument rechecked.
- Corollary 2: χ_{−d} primitive, real, odd, conductor d; w = 2 for d > 4; 0.11/π = 0.0350141 ≥ 0.0350.

## Code audit of cert.py

- Rigor: python-flint arb at 256 bits; B, 0.2, E_UP and the targets are balls around exact decimals; parameters are
  exact rationals n/1000; every PASS is a certain ball comparison (val > target, k1 < 1); floats only choose parameters.
- F (C:42–60) and Fup (C:105–117) equal Prop. 10 and Prop. 11 term by term; E_UP = −1.33258227573 lies above E.
- Coverage: cells [i/1000, j/1000] contiguous from round(1000Λ_lo) to round(1000Λ_hi); starts exactly at 3, 3.13
  (round(3130.0000000000005) = 3130) and 4; tails start exactly at 40 and 100; Λ ≤ 3 by one value at log q = e³ (F
  decreases in log q). No double-above-decimal gap. Class thresholds: cells start at the certified e₀ and the printed
  decimal is rounded up from e^{e₀} (2.957·10¹⁰ ≥ e^{24.11} = 2.9569·10¹⁰, 7.976·10¹³ ≥ e^{32.01} = 7.9757·10¹³).
- Independent recomputation (lwcheck.py, mpmath interval arithmetic at 256 bits, own formula from PROOF.md, own
  optimiser, own adaptive cells): the cert's tightest cells give the same margins to 15 digits (A 0.017384416090759,
  B 0.000471199722710, B2 0.022785041674912, upper 0.018771131531913); 60 random cert cells pass with own parameters;
  all covers, tails, thresholds and both B values pass (lwcheck_B.log, lwcheck_78.log).

## Falsification (no counterexample)

| Test | Range | Result | Bound | Log |
| --- | --- | --- | --- | --- |
| exact h(−d) by reduced forms, all fundamental −d | 4 < d ≤ 10⁷ (3,039,630 discriminants) | min h·max(3, log log d)/√d = 0.18985 (d = 991027, h = 63); largest d with h ≤ 100 is 2383747 (= Watkins) | 0.0350 | classno_1e7.log |
| L(1, χ_{−d}) = πh/√d, same range | d ≤ 10⁷ | min L(1)·max(3, log log d) = 0.596 | 0.11 | classno_1e7.log |
| all nonprincipal χ of prime conductor (FFT, digamma formula) | 3 ≤ q ≤ 20000, 21,166,667 characters | min \|L(1, χ)\| complex 0.2517 (q = 19543), real 0.2461 (q = 163) | 0.11/3 = 0.0367 | complexL.log |
| PARI qfbclassno = quadclassunit, just above each threshold, 200 consecutive + 25 with χ(p) = −1 for p ≤ 31 | d₀ = 6.973·10⁷ … 8.137·10¹³ | min h ≥ 550 ≫ H; min ratio 0.1976; no mismatch | H, 0.0350 | thresh.log |
| (e) where its bound is smallest (Λ = 3, bound 45), 25 d with χ(p) = +1 for p ≤ 31 | d ≥ e^{e³} | max L(1) = 7.00 | 45 | sidemax.log |
| explicit-formula residuals for complex χ (Lemmas 2–7) | q = 5, 7, 11, 13, 28 characters | all inside [−κZ − T, κZ + T], Z > 0 | — | cplxsanity.log |
| transitions | Λ = 3, 3.13, 4 | Λ ≤ 3 regime margin 0.0179 (even 0.11/Λ holds on [2.99, 3], 0.0146); 0.11/Λ from [3, 3.001], 0.0178; 0.115 fails at [3.125, 3.126] (−0.00007), holds at [3.129, 3.13] (+0.00074); 0.116 margin 0.116 at [3.99, 4] | — | trans.log |

The numerical tests sit far from the bounds (factor 5–7): the proved bound is 0.11/max(3, log log q), observed minima
are 0.60/max(3, log log q) (real) and 0.76/3 (complex) at these sizes. They exclude gross errors (sign, factor 2, wrong parity term), not small ones; the
small ones are excluded by the referee and the independent recomputation.


## Sensitivity to the one numerical input (Dusart)

Dusart's 2010 bounds use RH for Gourdon's first 10¹³ zeros (TeX l.133; height ≈ 2.445·10¹²), now rigorous by Platt–Trudgian
(2021, 3·10¹²); dusart_check.py confirms both bounds at every prime ≤ 10⁸ in their ranges (ratios ≤ 0.94).
If the error terms are weakened to the classical Rosser–Schoenfeld shapes, 1/(10L²) + 4/(15L³) → 1/(2L²) and 0.2/L₁ + 0.2/L₁² → 1/(2L₁) (lwweak.py, independent
covers): (a) still holds (worst margin 0.0152), 0.116 from Λ = 4 holds (0.0203), the Λ ≤ 3 value holds (0.0153), ω(8) < 0;
only the 0.115 threshold moves from 3.13 to 3.14. So only (b)'s decimal 3.13 needs Dusart's exact constants.

## Reviewer (one subagent, the most delicate step pasted inline)

Step: Hadamard identity (signs, ½ψ term), the involution ρ ↔ 1 − ρ̄ for complex χ, Lemma 6 (pairs) with fixed points and
multiplicities, and the loop (Z-bound, (★), where κ₁ < 1 and κ₀ ≥ 0 are used). Verdict: CORRECT ("pairing nonnegative sums
is legitimate", "κ₁ < 1 permits division without reversing that inequality; κ₀ ≥ 0 gives −κ₀Z ≥ −κ'(C − Re P₁)").

## Dependents (all unaffected: every bound they use holds)

discoveries/class-number-5000.md, class-number-10000.md, unit-generated-orders.md; DISCOVERIES.md l.11; COROLLARIES.md
l.37, 39, 40, 44; TREE.md l.15, 378; Riemann/Tree/R342, R358, R370, R407, R411; Riemann/Deep/ImqClass, ImqClass2, Exp48,
KoppLag, HP141, ExplProg, Harvest1–4 and 6, HPAudit (its table row). Text fixes only: S1 also in R358.lean l.8–9;
S2 also in HP141/REPORT.md l.25 ("limit 0.11694"). Full list: dependents.txt.

## Not checked

- Davenport Ch. 12 (17)–(18) and the contour bounds were checked at statement level from the standard form, not
  against the book's page; Dusart's proofs were not re-derived (statements only, plus the numeric range check above).
- The Lean half-plane itself: taken from HPAudit (AUDIT PASSED, kernels run in third-party CI).
- No Lean formalization of this chain (out of scope: "No Lean").

## Files (Riemann/Deep/LittlewoodAudit/)

FULLPROOF.md (the chain, line references P: PROOF.md, C: cert.py, D: Dusart TeX); lwcheck.py + lwcheck_B.log, lwcheck_78.log
(+ .err: /usr/bin/time and memcap peaks, ≤ 16 MB); trans.py + trans.log; lwweak.py + lwweak.log; classno.py +
classno_1e7.log (peak 55 MB); complexL.py + complexL.log (41 MB); thresh.gp + thresh.log; sidemax.gp + sidemax.log;
cplxsanity.py + cplxsanity.log; dependents.txt. Every run under .scratch/memcap.sh 290, one process at a time.

RH relevance: N11
