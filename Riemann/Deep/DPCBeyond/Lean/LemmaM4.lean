import Riemann.Deep.DPCBeyond.Lean.LemmaM

/-!
# dpcbeyond Lemma M, part 2: M4 (discrete moments over the zeros) and the kernel sum

M4_count, M4_moment: the exponent κ(1 − 1/w)² with the (|m| + 2) factor; key_ineq, inv_sq_sum, kernel_sum:
Σ_{n<M}(n+2)/(n+1)^A ≤ 4 for natural A ≥ 3 (see LemmaM.lean for the statement list).
-/

namespace DPCBLean

open Finset DPCLean

variable {ex lg : ℝ → ℝ}

/-- **M4's count**: Nc ≤ |B_a|^{1−1/w}(Σ_{ρ'}Ñ^w)^{1/w} (named covering + Hölder), |B_a| ≤ C T L e^{−κa},
Σ_{ρ'}(Ñ^{(m)}_{ρ'})^w ≤ C (|m|+2)^w T L give Nc ≤ C (|m|+2) T L e^{−κa(1−1/w)} (mm = |m|). -/
theorem M4_count (h : ExpLog ex lg) {Nc B Sw C T L κ a w mm : ℝ} (hC : 0 < C) (hT : 0 < T)
    (hL : 0 < L) (hw : 1 < w) (hm : 0 ≤ mm) (hN : Nc ≤ rp ex lg B (1 - 1 / w) * rp ex lg Sw (1 / w))
    (hB : B ≤ C * (T * L) * ex (-(κ * a))) (hS : Sw ≤ C * rp ex lg (mm + 2) w * (T * L)) :
    Nc ≤ C * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w))) := by
  have hw0 : 0 < w := by linarith only [hw]
  have he1 : 0 ≤ 1 - 1 / w := by
    have h1 : 1 / w < 1 := (div_lt_one hw0).2 hw
    linarith only [h1]
  have hiw : 0 ≤ 1 / w := by positivity
  have hTL : 0 < T * L := mul_pos hT hL
  have hm2 : 0 < mm + 2 := by linarith only [hm]
  have r1 : rp ex lg B (1 - 1 / w) ≤
      rp ex lg C (1 - 1 / w) * rp ex lg (T * L) (1 - 1 / w) * ex (-(κ * a * (1 - 1 / w))) := by
    calc rp ex lg B (1 - 1 / w) ≤ rp ex lg (C * (T * L) * ex (-(κ * a))) (1 - 1 / w) :=
          rp_le_rp h hB he1
      _ = rp ex lg C (1 - 1 / w) * rp ex lg (T * L) (1 - 1 / w) * rp ex lg (ex (-(κ * a))) (1 - 1 / w) := by
          rw [rp_mul h (mul_pos hC hTL) (h.pos _), rp_mul h hC hTL]
      _ = rp ex lg C (1 - 1 / w) * rp ex lg (T * L) (1 - 1 / w) * ex (-(κ * a * (1 - 1 / w))) := by
          rw [rp_ex h, show (1 - 1 / w) * -(κ * a) = -(κ * a * (1 - 1 / w)) by ring]
  have r2 : rp ex lg Sw (1 / w) ≤ rp ex lg C (1 / w) * (mm + 2) * rp ex lg (T * L) (1 / w) := by
    calc rp ex lg Sw (1 / w) ≤ rp ex lg (C * rp ex lg (mm + 2) w * (T * L)) (1 / w) := rp_le_rp h hS hiw
      _ = rp ex lg C (1 / w) * rp ex lg (rp ex lg (mm + 2) w) (1 / w) * rp ex lg (T * L) (1 / w) := by
          rw [rp_mul h (mul_pos hC (rp_pos h hm2 w)) hTL, rp_mul h hC (rp_pos h hm2 w)]
      _ = rp ex lg C (1 / w) * (mm + 2) * rp ex lg (T * L) (1 / w) := by
          rw [rp_rp h hm2, mul_one_div_cancel hw0.ne', rp_one h hm2]
  have eC : rp ex lg C (1 - 1 / w) * rp ex lg C (1 / w) = C := by
    rw [← rp_add h hC, show 1 - 1 / w + 1 / w = 1 by ring, rp_one h hC]
  have eT : rp ex lg (T * L) (1 - 1 / w) * rp ex lg (T * L) (1 / w) = T * L := by
    rw [← rp_add h hTL, show 1 - 1 / w + 1 / w = 1 by ring, rp_one h hTL]
  calc Nc ≤ rp ex lg B (1 - 1 / w) * rp ex lg Sw (1 / w) := hN
    _ ≤ (rp ex lg C (1 - 1 / w) * rp ex lg (T * L) (1 - 1 / w) * ex (-(κ * a * (1 - 1 / w)))) *
        (rp ex lg C (1 / w) * (mm + 2) * rp ex lg (T * L) (1 / w)) :=
        mul_le_mul r1 r2 (rp_nonneg h _ _) (mul_nonneg (mul_nonneg (rp_nonneg h _ _)
          (rp_nonneg h _ _)) (h.pos _).le)
    _ = (rp ex lg C (1 - 1 / w) * rp ex lg C (1 / w)) * (mm + 2) *
        (rp ex lg (T * L) (1 - 1 / w) * rp ex lg (T * L) (1 / w)) * ex (-(κ * a * (1 - 1 / w))) := by
        ring
    _ = C * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w))) := by rw [eC, eT]

/-- **M4's moment**: Σ_ρ (n^{(a)}_m(x_ρ))^p ≤ Nc^{1−1/w}(Σ_ρ n_m(x_ρ)^{pw})^{1/w} (named Hölder), Nc from M4_count,
Σ_ρ n_m(x_ρ)^{pw} ≤ C' T L (Lz DoubleSumZeros) give ≤ (C+1)(C'+1)(|m|+2) T L e^{−κa(1−1/w)²}. -/
theorem M4_moment (h : ExpLog ex lg) {I Nc Sp C C' T L κ a w mm : ℝ} (hC : 0 ≤ C) (hC' : 0 ≤ C')
    (hT : 0 < T) (hL : 0 < L) (hw : 1 < w) (hm : 0 ≤ mm)
    (hI : I ≤ rp ex lg Nc (1 - 1 / w) * rp ex lg Sp (1 / w))
    (hN : Nc ≤ C * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w)))) (hS : Sp ≤ C' * (T * L)) :
    I ≤ (C + 1) * (C' + 1) * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w) ^ 2)) := by
  have hw0 : 0 < w := by linarith only [hw]
  have hiw1 : 1 / w ≤ 1 := by rw [div_le_one hw0]; exact hw.le
  have he1 : 0 ≤ 1 - 1 / w := by linarith only [hiw1]
  have he11 : 1 - 1 / w ≤ 1 := by
    have h1 : 0 ≤ 1 / w := by positivity
    linarith only [h1]
  have hiw : 0 ≤ 1 / w := by positivity
  have hTL : 0 < T * L := mul_pos hT hL
  have hm2 : 0 < mm + 2 := by linarith only [hm]
  have hK : 0 < (C + 1) * (mm + 2) := mul_pos (by linarith only [hC]) hm2
  have hK1 : 1 ≤ (C + 1) * (mm + 2) :=
    one_le_mul_of_one_le_of_one_le (by linarith only [hC]) (by linarith only [hm])
  have hC1 : 1 ≤ C' + 1 := by linarith only [hC']
  have r1 : rp ex lg Nc (1 - 1 / w) ≤
      (C + 1) * (mm + 2) * rp ex lg (T * L) (1 - 1 / w) * ex (-(κ * a * (1 - 1 / w) ^ 2)) := by
    have hN' : Nc ≤ (C + 1) * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w))) := by
      have h2 : C * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w))) ≤
          (C + 1) * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w))) := by
        apply mul_le_mul_of_nonneg_right _ (h.pos _).le
        apply mul_le_mul_of_nonneg_right _ hTL.le
        exact mul_le_mul_of_nonneg_right (by linarith only) hm2.le
      exact le_trans hN h2
    calc rp ex lg Nc (1 - 1 / w) ≤
          rp ex lg ((C + 1) * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w)))) (1 - 1 / w) :=
          rp_le_rp h hN' he1
      _ = rp ex lg ((C + 1) * (mm + 2)) (1 - 1 / w) * rp ex lg (T * L) (1 - 1 / w) *
          rp ex lg (ex (-(κ * a * (1 - 1 / w)))) (1 - 1 / w) := by
          rw [rp_mul h (mul_pos hK hTL) (h.pos _), rp_mul h hK hTL]
      _ ≤ (C + 1) * (mm + 2) * rp ex lg (T * L) (1 - 1 / w) * ex (-(κ * a * (1 - 1 / w) ^ 2)) := by
          rw [rp_ex h, show (1 - 1 / w) * -(κ * a * (1 - 1 / w)) = -(κ * a * (1 - 1 / w) ^ 2) by ring]
          apply mul_le_mul_of_nonneg_right _ (h.pos _).le
          apply mul_le_mul_of_nonneg_right _ (rp_nonneg h _ _)
          have h2 := rp_le_max h (x := (C + 1) * (mm + 2)) he1 he11
          rwa [max_eq_right hK1] at h2
  have r2 : rp ex lg Sp (1 / w) ≤ (C' + 1) * rp ex lg (T * L) (1 / w) := by
    have hS' : Sp ≤ (C' + 1) * (T * L) := by
      have h2 : C' * (T * L) ≤ (C' + 1) * (T * L) := mul_le_mul_of_nonneg_right (by linarith only) hTL.le
      linarith only [hS, h2]
    calc rp ex lg Sp (1 / w) ≤ rp ex lg ((C' + 1) * (T * L)) (1 / w) := rp_le_rp h hS' hiw
      _ = rp ex lg (C' + 1) (1 / w) * rp ex lg (T * L) (1 / w) := rp_mul h (by linarith only [hC']) hTL _
      _ ≤ (C' + 1) * rp ex lg (T * L) (1 / w) := by
          apply mul_le_mul_of_nonneg_right _ (rp_nonneg h _ _)
          have h2 := rp_le_max h (x := C' + 1) hiw hiw1
          rwa [max_eq_right hC1] at h2
  have eT : rp ex lg (T * L) (1 - 1 / w) * rp ex lg (T * L) (1 / w) = T * L := by
    rw [← rp_add h hTL, show 1 - 1 / w + 1 / w = 1 by ring, rp_one h hTL]
  calc I ≤ rp ex lg Nc (1 - 1 / w) * rp ex lg Sp (1 / w) := hI
    _ ≤ ((C + 1) * (mm + 2) * rp ex lg (T * L) (1 - 1 / w) * ex (-(κ * a * (1 - 1 / w) ^ 2))) *
        ((C' + 1) * rp ex lg (T * L) (1 / w)) :=
        mul_le_mul r1 r2 (rp_nonneg h _ _) (by
          have := rp_nonneg h (T * L) (1 - 1 / w)
          have := (h.pos (-(κ * a * (1 - 1 / w) ^ 2))).le
          positivity)
    _ = (C + 1) * (C' + 1) * (mm + 2) * (rp ex lg (T * L) (1 - 1 / w) * rp ex lg (T * L) (1 / w)) *
        ex (-(κ * a * (1 - 1 / w) ^ 2)) := by ring
    _ = (C + 1) * (C' + 1) * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w) ^ 2)) := by rw [eT]

lemma key_ineq {x : ℝ} (hx : 1 ≤ x) : 1 / (x + 1) ^ 2 ≤ 1 / x - 1 / (x + 1) := by
  have hx0 : 0 < x := by linarith only [hx]
  have e : 1 / x - 1 / (x + 1) = 1 / (x * (x + 1)) := by
    rw [div_sub_div _ _ hx0.ne' (by positivity)]
    congr 1
    ring
  rw [e]
  apply one_div_le_one_div_of_le (by positivity)
  have e2 : (x + 1) ^ 2 - x * (x + 1) = x + 1 := by ring
  linarith only [e2, hx]

/-- Σ_{n ≤ k} 1/(n+1)² ≤ 2 − 1/(k+1). -/
lemma inv_sq_sum (k : ℕ) : ∑ n ∈ range (k + 1), 1 / ((n : ℝ) + 1) ^ 2 ≤ 2 - 1 / ((k : ℝ) + 1) := by
  induction k with
  | zero => norm_num
  | succ k ih =>
    rw [sum_range_succ]
    have hk : (1 : ℝ) ≤ (k : ℝ) + 1 := by
      have := (Nat.cast_nonneg k : (0 : ℝ) ≤ k)
      linarith only [this]
    have key := key_ineq hk
    push_cast
    linarith only [ih, key]

/-- **M4's kernel sum**: Σ_{n<M} (n+2)/(n+1)^A ≤ 4 for natural A ≥ 3 (n = |m|; over m ∈ ℤ at most twice this). -/
theorem kernel_sum (A : ℕ) (hA : 3 ≤ A) (M : ℕ) :
    ∑ n ∈ range M, ((n : ℝ) + 2) / ((n : ℝ) + 1) ^ A ≤ 4 := by
  have hterm : ∀ n ∈ range M, ((n : ℝ) + 2) / ((n : ℝ) + 1) ^ A ≤ 2 * (1 / ((n : ℝ) + 1) ^ 2) := by
    intro n _
    have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have h1 : (0 : ℝ) < (n : ℝ) + 1 := by linarith only [hn]
    have h3 : ((n : ℝ) + 1) ^ 3 ≤ ((n : ℝ) + 1) ^ A := pow_le_pow_right₀ (by linarith only [hn]) hA
    rw [mul_one_div, div_le_div_iff₀ (by positivity) (by positivity)]
    have e4 : 2 * ((n : ℝ) + 1) ^ 3 - ((n : ℝ) + 2) * ((n : ℝ) + 1) ^ 2 = ((n : ℝ) + 1) ^ 2 * n := by
      ring
    have h5 := mul_nonneg (sq_nonneg ((n : ℝ) + 1)) hn
    have e6 : 1 * ((n : ℝ) + 1) ^ A = ((n : ℝ) + 1) ^ A := one_mul _
    linarith only [e4, h5, h3, e6]
  have hsum : ∑ n ∈ range M, 1 / ((n : ℝ) + 1) ^ 2 ≤ 2 := by
    rcases Nat.eq_zero_or_pos M with hM | hM
    · rw [hM]
      simp
    · obtain ⟨k, rfl⟩ : ∃ k, M = k + 1 := ⟨M - 1, by omega⟩
      have h1 := inv_sq_sum k
      have h2 : 0 ≤ 1 / ((k : ℝ) + 1) := by positivity
      linarith only [h1, h2]
  calc ∑ n ∈ range M, ((n : ℝ) + 2) / ((n : ℝ) + 1) ^ A ≤ ∑ n ∈ range M, 2 * (1 / ((n : ℝ) + 1) ^ 2) :=
        sum_le_sum hterm
    _ = 2 * ∑ n ∈ range M, 1 / ((n : ℝ) + 1) ^ 2 := by rw [mul_sum]
    _ ≤ 4 := by linarith only [hsum]

end DPCBLean

