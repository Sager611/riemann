AUDIT PASSED: Λ_f > 218/125 holds as written — Thm mixed follows from the paper's shifted moment by its own Cauchy step for every fixed μ, ν, the binomial recombination gives the uniform measure exactly, the reduction (zero normalisation, multiplicities, Re/Im, telescoping, subconvex end terms) is correct (reviewer CORRECT), the certified a_j are the ones in H, the Arb cover reaches every gap in (0, π], H is exact, both scripts rerun byte-identically, (V₆) survives an independent exact-arithmetic falsification test, and Δ's mean values match the theorem's main term to ≤ 1%.

# gl2gapaudit: line-by-line audit of the GL(2) large-gap record Λ_f > 1.744 (gapfamily, 20 Nov 2026; nodes N8, N12)
Claim audited: Riemann/Deep/GapFamily/REPORT.md l.1–14 ('RECORD', 'Statement') and the RECORDS.md line '(N8/N12) GL(2) large gaps'.
Full chain with source lines: FULLPROOF.md (§0 objects, §1 Thm mixed, §2 ∫|g^{(j)}|², §3 reduction, §4 certificate, §5 numerics, §6 used as printed).
TeX: arXiv 1410.7765 (single file; line numbers as GapFamily's).

## Steps
| step | status | evidence |
|---|---|---|
| Λ_f and its normaliser (l.144–162) | holds | N(T,f) = (T/π)log(√qT/2πe) + O(log 𝔮) (IK Thm 5.8, degree 2, γ > 0); log(√qγ) = log T + O_q(1) and ≥ log T on [T, 2T]; Δ: 245 sign changes of Z on [200, 400] vs 244.89 predicted (zeros.gp) |
| Thm mixed hypotheses (l.239–247) | as stated | primitive f, GL(2)/Q, level q fixed, holomorphic or Maass (l.130–132); all μ, ν ≥ 0; error O(μ!ν!T(log T)^{μ+ν}); we use fixed μ, ν ≤ 6 only |
| Thm mixed from Thm shiftedmoment (l.400–478) | holds | re-derived (FULLPROOF §1): F(z;T) series, Cauchy on radius 1/log T (needs shifts ≤ 2/log T, complex: as printed), ∂^{μ+ν}F at 0 |
| ∫_T^{2T}|g^{(j)}|² = 2c_fTL^{2j+1}/(2j+1)(1 + O_j(1/L)), g = e^{itL}L, L = log T | holds | binomial expansion of g^{(j)}; Σ_{k,l}C(j,k)C(j,l)(−1)^{k+l}2^{k+l+1}/(k+l+1) = 2/(2j+1) exact for j ≤ 6 (indep.py (c)); j ≤ 1 = paper l.354–362 |
| scale L = log T (not log(√qT/2π)) | legitimate | both the zero density (1/π)ℓ(t) and the spectral half-width are ℓ(t) = log(√qt/2π) = L + O_q(1): leading moments identical; corrections O(1/L) vanish in the limit |
| reduction (gaps ≤ Aπ/L, distinct zeros, partial gaps, Re/Im, telescoping, end pieces) | holds | FULLPROOF §3; reviewer CORRECT; end pieces and end jets O(T^{1−2δ+ε}) by t-aspect subconvexity (l.330–337; Blomer–Harcos: all levels, Maass) |
| certified form = form in H | identical | cert5.py writes the rationals it certifies; gl2crit.py reads them; own rationalisation of SRC = certified a, φ (13/13) |
| cover starts at or below the smallest gap | yes, at 0 | g < 2.108718 = (π/√X*).lower(): every mode above X* = 2.2195, termwise PSD (C ≻ 0, leading coefficients > 0); [2.108718, π]: 110 Arb cells |
| H exact | yes | Fraction arithmetic (gl2crit.py) = 7.7046859308427·10⁻⁴; second formula −(1/Λ)∫_0^Λ p(v²)dv at 60 digits agrees to 4·10⁻⁶³; H(1.7445) = 2.6·10⁻⁴ > 0 (room for A) |
| reruns from a copy (run/) | identical | cert5.py (JJ=6) and gl2crit.py: both JSON byte-identical (cmp); CERTIFIED, CRITERION HOLDS; code sha256 = GapCertAudit REPORT2's audited files |
| falsify (V₆) directly | no counterexample | indep.py (d): form on u^{m+1}(1 − u), m < 22, exact rationals, 13 gap lengths 0.01 … 3.1415926: all LDLᵀ pivots > 0; min form/‖f‖² 4.77·10⁻⁴ at g ≈ π, larger below |
| mean values, Δ (level 1, weight 12) | normalisation right | PARI lfun on [T, 2T], μ, ν ≤ 2: I_μν/(recipe main term) = 0.975–0.988 (T = 200), 1.002–1.009 (T = 400); leading term with log T: 1.11 … 0.55 (T = 400), the expected 1/log T terms |

## Findings (none load-bearing)
- 1410.7765 slips: (i) the dual factor (t/2π)^{−2(α+β)} (l.261, l.1328–1365) lacks √q: the difference is O(T) for fixed q; (ii) R is
  defined with ∫_0^T (l.428), the proof is for [T, 2T] and the closing "dyadic intervals" sentence (l.477) is a leftover; (iii)
  "(t/2π)^{−2z} = T^{−2z}(1 + O(1/log T))" (l.405–406) times c_f/z is not O(1) as z → 0; the true bound O(|z|) makes it O(1)
  uniformly; (iv) "Λ_f ≥ limsup κ_T" (l.298) omits the normalisation κ_T log T/π. All four are wording; statements stand.
- GapFamily REPORT: (i) "−2c_f H(A)TL ≥ −o(TL)" omits the common positive factor λ = L/A (the sum is λ·2c_fTL(−H(A) + O(1/L)));
  (ii) the gap straddling T needs slack Λ₀ < B < A (its normaliser is at γ < T): reviewer's clarification, FULLPROOF §3. Wording only.
- Finite-T behaviour: with the twist log T the moments of g approach the asymptotic slowly (ratio 1.50 at j = 2, T = 400), because
  log T sits log 2π above the spectral centre ℓ(t) = log(t/2π); the full main term explains this to 0.2%. Irrelevant for the limit.
- The per-gap form at 218/125 has continuity room up to A ≈ 1.7447 (H(1.7445) > 0 > H(1.745)); the float method ceiling is 1.7452.

## Not re-refereed (used as printed)
- Thm shiftedmoment's error analysis (§7–8 of 1410.7765 v4 = J. Math. Anal. Appl. 2015, DOI 10.1016/j.jmaa.2015.04.007; l.1177–1544:
  approximate functional equation after Ramachandra, J_i terms, Lemmas weightsum, truncatedsum, stirling, GeneralMV). The printed
  error analysis is a sketch: one "representative" J_i (l.1374–1430), the others "follow similarly"; its absolute-value step writes
  n^{−1/2−α} for n^{−1/2−Re α}. Robustness: the Cauchy step only needs R = o(T log T) uniformly on |α|, |β| ≤ 2/log T (then
  ∂^μ∂^νR = o(TL^{μ+ν+1}) and every main term above survives), a standard mean-value statement (approximate functional equation +
  Montgomery–Vaughan; diagonal = the recipe's main terms). Checked: its main terms S_1, S_2 (l.1283–1365) hold for complex shifts ≪ 1/log T
  and their X-dependent parts are O(T) uniformly, incl. α + β → 0; Δ's data match its main term (≤ 1% at T = 400, zero shifts, μ, ν ≤ 2).
- Lemma 4's algebra (audited in GapCertAudit REPORT2 with sympy, same files by sha256); t-aspect subconvexity (cited theorems).

## Review
- reviewer (verdict-only; Steps 1–3 with inputs A (Thm shiftedmoment), B (V₆ + H), C (subconvexity) pasted inline): **CORRECT**.
  Points checked: Cauchy uniformity incl. α + β = 0 (removable), the log T scale (O(1) changes only lower-order terms), multiplicities,
  Re/Im and exact telescoping (the boundary products combine to Re(g^{(p)}conj g^{(q)})), end terms O(T^{1−2δ}). Clarification
  (straddling gap, B < A) adopted.

## Runs (from Riemann/Deep/GL2GapAudit; OMP etc. = 1; each under bash ../../../.scratch/memcap.sh 290 <s>)
| command | time | peak |
|---|---|---|
| cd run; JJ=6 SRC=logs/gl2_J6_1.744m.json uv run --no-project --with python-flint --with numpy python3 cert5.py gl2J6_1.744m 0 0 0 61 | 10 s | < memcap sampling (GapFamily: 22 MB) |
| cd run; uv run --no-project python3 gl2crit.py gl2J6_1.744m 218/125 | < 1 s | — |
| uv run --no-project --with mpmath python3 indep.py run/logs/gl2_J6_1.744m.json run/logs/cert5_gl2J6_1.744m.json 22 (also 8) | 18 s | 14 MB |
| gp -q -s 120000000 run_T200.gp, run_T400.gp (read meanval.gp; laurent.gp alone: 11 s, 74 MB) | 36 s, 149 s | 78 MB (RSS 87 MB) |
| gp -q -s 64000000 zeros.gp (zero count, T = 200) | 31 s | 12 MB |
| T = 1000 attempt: killed at the 900 s limit (lfun cost grows with height); not retried | — | 78 MB |

## Files
- Riemann/Deep/GL2GapAudit/: REPORT.md, FULLPROOF.md, indep.py, meanval.gp, laurent.gp, zeros.gp, run_T200.gp, run_T400.gp; run/ (copies of
  cert5.py, exact5.py, exact6.py, gl2crit.py; run/logs: float form, cert JSON, criterion JSON and logs); logs/ (indep_M22.log,
  laurent.log, meanval_T200.log, meanval_T400.log, zeros_T200.log). .scratch/gl2gapaudit/STATE.md.

RH relevance: N8, N12
