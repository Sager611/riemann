import Riemann.Deep.LiRate.Lean.Asym

/-!
# LiRate: (H_τ) reduced to its ingredients (proofs.md Theorem 1, steps 4, 5, 7)

`tail_le` is the favourable-sign step of proofs.md step 5 (line 50): the zeros between `T` and
`X` enter `Σ_{γ'} B(θ_γ, θ_{γ'})` through `-Φ(θ_γ - θ_{γ'}) ≤ 0` and are dropped. `trendTail_of`
derives `TrendTail` from (i) the tiny terms (rapid decay of `Φ` off `0`, lines 48–49, 52) and
(ii) the trend estimate `-Σ_{γ≤T} τ(θ_γ) ≤ (1 + η/2) T³L²/(6π)` (step 4 + step 7: `τ(θ) =
log(2πθ)/θ² + O(log N)`, `τ_T ~ T³L²/(6π)` by partial summation with `N(T)`).
-/

open Finset Real Filter

namespace LiRate

/-- Favourable sign (proofs.md step 5, line 50): for `J ⊆ X`, since `Φ ≥ 0`,
`tail_J(a) ≤ Σ_{k ∈ X∖J} (2Φ(a) - Φ(a + θ_k)) + tail_X(a)`. -/
lemma tail_le {κ : Type*} [DecidableEq κ] (N : ℕ) (w : ℕ → ℝ) (hΦ : ∀ x, 0 ≤ Phi N w x)
    (θ : κ → ℝ) (lam : ℕ → ℝ) (J X : Finset κ) (hJX : J ⊆ X) (a : ℝ) :
    tail N w θ lam J a ≤
      ∑ k ∈ X \ J, (2 * Phi N w a - Phi N w (a + θ k)) + tail N w θ lam X a := by
  unfold tail
  have hsplit : ∑ k ∈ X \ J, B N w a (θ k) + ∑ k ∈ J, B N w a (θ k) = ∑ k ∈ X, B N w a (θ k) :=
    sum_sdiff hJX
  have hle : ∑ k ∈ X \ J, B N w a (θ k) ≤ ∑ k ∈ X \ J, (2 * Phi N w a - Phi N w (a + θ k)) :=
    sum_le_sum fun k _ => by unfold B; linarith [hΦ (a - θ k)]
  linarith

/-- (H_τ) from (i) the tiny terms and (ii) the trend estimate; `Xs N` = the zeros with `γ ≤ X`
(proofs.md: `X = N³`). -/
theorem trendTail_of {κ : Type*} [DecidableEq κ] (w : ℕ → ℕ → ℝ) (θ : κ → ℝ) (lam t : ℕ → ℝ)
    (Z : ℝ → Finset κ) (η : ℝ) (hΦ : ∀ N x, 0 ≤ Phi N (w N) x) (Xs : ℕ → Finset κ)
    (hsub : ∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 → Z T ⊆ Xs N)
    (htiny : ∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 →
      2 * (Z T).card * ∑ j ∈ Z T, Phi N (w N) (θ j) +
        ∑ j ∈ Z T, (∑ k ∈ Xs N \ Z T, (2 * Phi N (w N) (θ j) - Phi N (w N) (θ j + θ k)) +
          tail N (w N) θ lam (Xs N) (θ j)) ≤
        η / 2 * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π)))
    (htrend : ∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 →
      -∑ j ∈ Z T, trend N (w N) t (θ j) ≤
        (1 + η / 2) * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π))) :
    TrendTail w θ lam t Z η := by
  filter_upwards [hsub, htiny, htrend] with N h1 h2 h3
  intro T hT1 hT2
  have hs := h1 T hT1 hT2
  have htl : ∑ j ∈ Z T, tail N (w N) θ lam (Z T) (θ j) ≤
      ∑ j ∈ Z T, (∑ k ∈ Xs N \ Z T, (2 * Phi N (w N) (θ j) - Phi N (w N) (θ j + θ k)) +
        tail N (w N) θ lam (Xs N) (θ j)) :=
    sum_le_sum fun j _ => tail_le N (w N) (hΦ N) θ lam (Z T) (Xs N) hs (θ j)
  have h2' := h2 T hT1 hT2
  have h3' := h3 T hT1 hT2
  unfold tauPrime
  linarith

end LiRate
