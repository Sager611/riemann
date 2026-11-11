import Riemann.Deep.LiRate.Lean.Asym

/-!
# LiRate: the Fejér weights satisfy (H_w) (non-vacuity check of `GoodWeights`)

`fejer N n = (1 - n/(N+1))₊ = φ(n/(N+1))`, `φ(u) = (1 - |u|)₊ = ψ * ψ` for `ψ = 1_{(-1/2,1/2)}`: the
endpoint of proofs.md's family (line 36, `ψ → 1_{(-1/2,1/2)}`, `∫φ → 1`). Here `Φ` is the Fejér
kernel, `(N+1) Φ(x) = |Σ_{k≤N} e^{ikx}|² ≥ 0`, and `Φ(0) = N + 1`. (Its `Φ` decays only like
`1/(N x²)`, so it does not satisfy (H_tiny); proofs.md uses smooth `ψ` for that.)
-/

open Finset Real Filter

namespace LiRate

/-- Fejér weights. -/
noncomputable def fejer (N n : ℕ) : ℝ := max 0 (1 - n / (N + 1))

lemma fejer_eq (N n : ℕ) (hn : n ≤ N) : fejer N n = 1 - n / (N + 1) := by
  unfold fejer
  apply max_eq_right
  rw [sub_nonneg, div_le_one (by positivity)]
  exact_mod_cast (by omega : n ≤ N + 1)

lemma cos_reflect (N : ℕ) (x : ℝ) :
    ∑ k ∈ range (N + 1), cos (((N : ℝ) + 1 - k) * x) = ∑ n ∈ Icc 1 (N + 1), cos (n * x) := by
  apply sum_nbij' (fun k => N + 1 - k) (fun n => N + 1 - n)
  · intro k hk; simp only [mem_range] at hk; simp only [mem_Icc]; omega
  · intro n hn; simp only [mem_Icc] at hn; simp only [mem_range]; omega
  · intro k hk; simp only [mem_range] at hk; omega
  · intro n hn; simp only [mem_Icc] at hn; omega
  · intro k hk
    simp only [mem_range] at hk
    congr 1
    rw [Nat.cast_sub (by omega : k ≤ N + 1)]
    push_cast
    ring

/-- `(Σ_{k≤N} cos kx)² + (Σ_{k≤N} sin kx)² = (N+1) + 2 Σ_{n=1}^N (N+1-n) cos nx`. -/
lemma fejer_sq (N : ℕ) (x : ℝ) :
    (∑ k ∈ range (N + 1), cos (k * x)) ^ 2 + (∑ k ∈ range (N + 1), sin (k * x)) ^ 2 =
      (N + 1) + 2 * ∑ n ∈ Icc 1 N, ((N : ℝ) + 1 - n) * cos (n * x) := by
  induction N with
  | zero => simp
  | succ N ih =>
    rw [sum_range_succ (fun k : ℕ => cos (k * x)), sum_range_succ (fun k : ℕ => sin (k * x))]
    push_cast
    have hcross : cos (((N : ℝ) + 1) * x) * ∑ k ∈ range (N + 1), cos (k * x) +
        sin (((N : ℝ) + 1) * x) * ∑ k ∈ range (N + 1), sin (k * x) =
          ∑ n ∈ Icc 1 (N + 1), cos (n * x) := by
      rw [← cos_reflect, mul_sum, mul_sum, ← sum_add_distrib]
      refine sum_congr rfl fun k _ => ?_
      rw [← cos_sub]
      congr 1
      ring
    have hpy := sin_sq_add_cos_sq (((N : ℝ) + 1) * x)
    have hG : ∑ n ∈ Icc 1 (N + 1), ((N : ℝ) + 1 + 1 - n) * cos (n * x) =
        ∑ n ∈ Icc 1 N, ((N : ℝ) + 1 - n) * cos (n * x) + ∑ n ∈ Icc 1 (N + 1), cos (n * x) := by
      rw [sum_Icc_succ_top (by omega : 1 ≤ N + 1), sum_Icc_succ_top (by omega : 1 ≤ N + 1)]
      have h1 : ∑ n ∈ Icc 1 N, ((N : ℝ) + 1 + 1 - n) * cos (n * x) =
          ∑ n ∈ Icc 1 N, ((N : ℝ) + 1 - n) * cos (n * x) + ∑ n ∈ Icc 1 N, cos (n * x) := by
        rw [← sum_add_distrib]
        exact sum_congr rfl fun n _ => by ring
      rw [h1]
      push_cast
      ring
    linear_combination ih + 2 * hcross + hpy - 2 * hG

lemma phi_fejer (N : ℕ) (x : ℝ) :
    ((N : ℝ) + 1) * Phi N (fejer N) x =
      (N + 1) + 2 * ∑ n ∈ Icc 1 N, ((N : ℝ) + 1 - n) * cos (n * x) := by
  have hN1 : (N : ℝ) + 1 ≠ 0 := by positivity
  have h : ∑ n ∈ Icc 1 N, ((N : ℝ) + 1 - n) * cos (n * x) =
      ((N : ℝ) + 1) * ∑ n ∈ Icc 1 N, fejer N n * cos (n * x) := by
    rw [mul_sum]
    refine sum_congr rfl fun n hn => ?_
    rw [mem_Icc] at hn
    rw [fejer_eq N n hn.2]
    field_simp
  rw [h]
  unfold Phi
  rw [fejer_eq N 0 (Nat.zero_le _)]
  simp only [Nat.cast_zero, zero_div, sub_zero]
  ring

/-- Fejér positivity: `Φ ≥ 0`. -/
lemma fejer_phi_nonneg (N : ℕ) (x : ℝ) : 0 ≤ Phi N (fejer N) x := by
  have h := phi_fejer N x
  rw [← fejer_sq] at h
  have hpos : (0 : ℝ) < N + 1 := by positivity
  exact (mul_nonneg_iff_of_pos_left hpos).1 (h ▸ by positivity)

/-- `Φ(0) = N + 1`. -/
lemma fejer_phi_zero (N : ℕ) : Phi N (fejer N) 0 = N + 1 := by
  have h := phi_fejer N 0
  rw [← fejer_sq] at h
  simp only [mul_zero, cos_zero, sin_zero, sum_const, card_range, nsmul_eq_mul, mul_one] at h
  have hpos : (N : ℝ) + 1 ≠ 0 := by positivity
  push_cast at h
  field_simp at h
  nlinarith [h]

/-- The Fejér weights satisfy (H_w) for every `η ≥ 0`. -/
theorem fejer_goodWeights (η : ℝ) (hη : 0 ≤ η) : GoodWeights fejer η := by
  refine ⟨fun N n => le_max_left _ _, fun N n => max_le zero_le_one ?_, fejer_phi_nonneg,
    Eventually.of_forall fun N => ?_⟩
  · have : (0 : ℝ) ≤ n / (N + 1) := by positivity
    linarith
  · rw [fejer_phi_zero]
    nlinarith [(Nat.cast_nonneg N : (0 : ℝ) ≤ N)]

end LiRate
