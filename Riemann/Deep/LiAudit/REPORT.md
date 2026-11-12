AUDIT PASSED — every LiRate claim (Lemma 0, Lemma A, Theorems 1–4, Corollaries 1–4, the Ω-result, RH ⟺ Σ m^δ y_m² < ∞) now has a complete proof (FULLPROOF.md) with no mathematical error; only documentation slips were found.

# LiAudit: line-by-line referee of the Li-coefficient remainder discovery (liaudit, 12 Nov 2026)

Audited: Riemann/Deep/LiRate/proofs.md (P), discoveries/li-remainder.md. Proofs: FULLPROOF.md (F0–F12). Code, logs: this folder.

## Verdict per claim
| Claim | Hypothesis | FULLPROOF | Verdict |
| --- | --- | --- | --- |
| Li's formula λ_n = Σ_γ 2(1 − cos nθ_γ), λ_{−n} = λ_n | RH | F0 | correct |
| Lemma 0: ¬RH ⇒ \|R(n)\| ≥ e^{cn} i.o. | none | F1 | correct |
| Lemma A: λ_n^{(0)} = t(n) − 1/4 + O(log n/n) | none | F2 | correct (−1/4 = 2∫r, r from θ_RS's −π/8) |
| Lemma A: exact decomposition with 7/4 | RH | F3 | correct, for every T > 0 not an ordinate |
| Thm 1: Σ_{n≤N}R² ≥ (c₁ − ε)N^{3/2}(log N)^{1/2}, c₁ = 0.055790 | RH only | F4 | correct |
| Cor 1: limsup \|R\|/(n log n)^{1/4} ≥ 0.28928; Ω((n log n)^{1/4}) unconditional | RH / none | F5 | correct |
| Cor 2: Σ m^{1/2}y_m² = ∞ | RH | F5 | correct |
| Thm 2: window sums ≪ N^{3/2}(log N)^{3/2}; ≪ N^{3/2}(log N)^{1/2} log log N | RH; RH + Montgomery | F6 | correct |
| Cor 3: Σ m^δ y_m² < ∞ ⟺ δ < 1/2; RH ⟺ (∃/∀ δ ∈ [0, ½)) convergence; δ = ½ never | RH / none | F7 | correct |
| Thm 3: \|R(n)\| ≤ (0.2021 + o(1))√n log n/√(log log n) | RH (+ CCM) | F8 | correct |
| Thm 4: smoothed mean of R is 7/4 + O_φ(log N/N) | RH | F9 | correct |
| Cor 4: Merlini–Sala–Sala's a·log n bound false for every a | none | F10 | correct |

## Where P was only a sketch (filled in FULLPROOF, no change of statement)
- Thm 1 step 4 (trend lemma): rigorous via the distributional FT of t(|x|), a symmetric cutoff around ξ and Poisson; error
  O_ψ(1 + N^{−2}θ^{−4} log(2/θ)), so O(1) per zero at the chosen T (P: O(log N), also enough).
- Thm 1 order of limits: ε → ψ (flat top, (∫ψ)³ ≥ 1 − ε/2c₁) → N₀(ψ, ε); no uniformity in ψ needed.
- Lemma A (A7): the log singularity of ρ at v = π needs a split at π − 1/n before integrating by parts.
- Thm 2(a): Montgomery's identity carries the weight 4/(4 + u²); pairs with |γ − γ'| > 1 are handled separately (harmless).
- Thm 4: the bound |G_N''| ≪ N³/t⁴ for t ≥ N needs Re Φ_+'(θ) ≪ N³θ (|sin x| ≤ |x|), not the crude N²t⁻³.

## Slips (none affects a theorem)
1. P Cor 4, "(Their evidence: n ≤ 30.)": Merlini–Sala–Sala computed n ≤ 31, quoted Maślanka's values to n = 3000 and read
   values off plots to n = 8000 (5080 → 16.83 vs 2 log n = 17.07). Our λ̃_n matches their computed values to 2·10⁻⁶.
2. Discovery Numerics: smoothed mean on (N/2, N) at N = 10⁴ is 1.7443, not 1.7515 (that is N = 3·10⁴, 1.7512).
3. P Thm 1 data "Σ_{n≤N}R² ≥ 1.7 c₁N^{3/2}(log N)^{1/2}, N ≤ 10⁵" is false at, e.g., N = 660 (ratio 1.472) and N = 1626 (1.476);
   over all 100 ≤ N ≤ 10⁵ the minimum is 1.472 (logs/minratio.txt; P sampled few N). Still > 1, so consistent with Theorem 1.
4. Lemma A's hypothesis γ₁ < T is unnecessary; Lemma 0's "Σ λ_n z^n has the singularities of log ξ(1/(1−z))" should read
   Σ λ_n z^{n−1} = d/dz log ξ(1/(1−z)) (same radius).
5. Discovery Method: "−1/4 from the Fejér kernel and the −π/8 constant": the Fejér kernel gives t(n) + O(1/n); −1/4 is the −π/8.
6. P's Barrier ("RH does not supply" cancellation in Z_T) is a heuristic remark, not a barrier theorem; the proved part is the range
   [1/4, 1/2] for the pointwise exponent.

## Numerical falsification (FULLPROOF F12; every process under memcap 290 MB)
- λ_n for n ≤ 10⁵ recomputed independently from Odlyzko's 10⁵ zeros (RH above them) with an explicit tail bound (zeros_lambda.py, 31 MB): agrees with the
  FFT rerun (97 MB) to 4.9·10⁻⁶ everywhere, within 2% of the error bound; λ_100000 = Johansson's value.
- Theorem 1 ratio Σ_{n≤N}R²/(c₁N^{3/2}√log N) over every 100 ≤ N ≤ 10⁵: 1.472–2.694; Thm 2 window ratio 0.11–0.32.
- Theorem 1's inequality chain evaluated with real zeros and λ_n at N = 10³–10⁵, two weights and the corner T = N^{1/2}/2: every
  inequality holds; actual Σ R² is 6.1× the finite-N bound at 10⁵. Lemma A constants, c₁, 0.289, κ₃ and Thm 2 consistent.
- Theorem 4: the error [Σφ(n/N)(R(n) − 7/4)]/log N stays bounded for N = 2·10³–10⁵ in every window tested, including narrow
  windows at the support edges: (0.5,1) |·| ≤ 1.4, (0.9,1) ≤ 26, (0.98,1) ≤ 65, (0.01,0.02) ≤ 1.4 (O_φ constant grows as the
  window narrows, as the proof predicts); logs/thm4_scan.txt.
- Cor 4 numbers reproduced exactly (first n = 5083, 7584 exceedances, max 4.1814 at 79302, a = 3, 4 from 46992, 79286).

## Inputs checked at source (FULLPROOF F11)
CCM 1309.1526 (S₁ bounds π/24, π/48; |S| ≤ ¼ log t/log log t); Lagarias math/0404394 (O(√n log n)); Voros math/0506326
(RH ⟺ λ_n = t(n) + o(n); 7/4 "should not count asymptotically"; (TS): λ̄_n − t(n) → 3/4) and math/0404213 (7/4);
Montgomery as restated by 2108.09258 (uniform |α| ≤ 1); Arias de Reyna: abstract only (OpenAlex), theorem numbers unchecked;
Merlini–Sala–Sala 1903.06233 (wording, λ_tiny = λ̃); Coffey math-ph/0505052, 0507042 (only an O(n^{1/2+ε}) conjecture; compatible).

## Review
reviewer subagent on Theorem 1 steps (3)–(7) pasted inline: VERDICT CORRECT (one slip, θ(T) ≥ 1/(1.01T) needs T ≥ 3: fixed).

## Not covered
Lean files (out of scope), heuristics H1–H3 of derivation.md, novelty.

## Files
FULLPROOF.md; zeros_lambda.py, audit_num.py, thm1_chain.py, lemmaA.py, mss_check.py, minratio.py, thm4_scan.py, fft_lambda.py (unmodified copy);
logs/: fft_r999.log, fft_trend.log, zeros_lambda.log, audit_num.txt, thm1_chain.txt, lemmaA.txt, mss_check.txt, minratio.txt, thm4_scan.txt.

RH relevance: N3
