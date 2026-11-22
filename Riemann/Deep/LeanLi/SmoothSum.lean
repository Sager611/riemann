module

public import Riemann.Deep.LeanLi.Thm2
public import Riemann.Deep.LeanLi.Smooth

/-!
# LeanLi: smoothed sums `E_j(x) = Σ_n φ(n/N) n^j e^{inx}` for `φ ∈ C_c^∞(0, ∞)`, and their decay at `x = θ(t)`
-/

@[expose] public section

open Real Filter Set Finset Complex
open scoped ContDiff

namespace LeanLi

/-- A smooth test function supported in `(0, B]`. -/
structure TestFn (φ : ℝ → ℝ) (B : ℕ) : Prop where
  smooth : ContDiff ℝ ∞ φ
  compact : HasCompactSupport φ
  pos : ∀ y, φ y ≠ 0 → 0 < y
  le : ∀ y, φ y ≠ 0 → y ≤ B

lemma testFn_of {φ : ℝ → ℝ} (hφ : ContDiff ℝ ∞ φ) (hc : HasCompactSupport φ) (hs : tsupport φ ⊆ Ioi 0) :
    ∃ B : ℕ, 1 ≤ B ∧ TestFn φ B := by
  obtain ⟨R, hR0, hR⟩ := exists_R_of_compact hc
  refine ⟨⌈R⌉₊ + 1, by omega, hφ, hc, fun y hy => hs (subset_tsupport _ hy), fun y hy => ?_⟩
  by_contra h
  push Not at h
  have h1 : (R : ℝ) ≤ ⌈R⌉₊ := Nat.le_ceil R
  push_cast at h
  have hy0 : 0 < y := by linarith
  exact hy (hR y (by rw [abs_of_pos hy0]; linarith))

lemma TestFn.eq_zero_of_nonpos {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) {y : ℝ} (hy : y ≤ 0) : φ y = 0 := by
  by_contra h'; exact absurd (h.pos y h') (not_lt.2 hy)

lemma TestFn.eq_zero_of_gt {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) {y : ℝ} (hy : (B : ℝ) < y) : φ y = 0 := by
  by_contra h'; exact absurd (h.le y h') (not_le.2 hy)

lemma TestFn.bound {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) : ∃ M : ℝ, 0 ≤ M ∧ ∀ y, |φ y| ≤ M := by
  obtain ⟨M, hM⟩ := exists_bound_iteratedDeriv h.smooth h.compact 0
  simp only [iteratedDeriv_zero] at hM
  exact ⟨max M 0, le_max_right _ _, fun y => (hM y).trans (le_max_left _ _)⟩

/-- Terms beyond `Bs` vanish. -/
lemma TestFn.vanish {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) {s : ℝ} (hs : 0 < s) {n : ℕ}
    (hn : ⌊(B : ℝ) * s⌋₊ < n) : φ (n / s) = 0 := by
  apply h.eq_zero_of_gt
  rw [lt_div_iff₀ hs]
  exact (Nat.floor_lt (by positivity)).1 hn

lemma TestFn.zero {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) (s : ℝ) : φ ((0 : ℕ) / s) = 0 :=
  h.eq_zero_of_nonpos (by simp)

lemma floor_cast_mul (B N : ℕ) : ⌊(B : ℝ) * (N : ℝ)⌋₊ = B * N := by
  rw [← Nat.cast_mul, Nat.floor_natCast]

/-- `E_j(x) = Σ_{n ≤ Bs} φ(n/s) n^j e^{inx}` at a real scale `s`. -/
noncomputable def Ej (φ : ℝ → ℝ) (B : ℕ) (s : ℝ) (j : ℕ) (x : ℝ) : ℂ :=
  ∑ n ∈ range (⌊(B : ℝ) * s⌋₊ + 1), ((φ (n / s) * (n : ℝ) ^ j : ℝ) : ℂ) * Complex.exp (I * ((n : ℝ) * x))

lemma Ej_eq_esum {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) {s : ℝ} (hs : 1 ≤ s) (j : ℕ) (x : ℝ) :
    Ej φ B s j x = (s : ℂ) ^ j * esum (fun n : ℤ => ((n : ℝ) / s) ^ j * φ ((n : ℝ) / s)) x := by
  have hs0 : 0 < s := by linarith
  have hN0 : s ≠ 0 := hs0.ne'
  unfold esum
  rw [tsum_eq_sum (s := (range (⌊(B : ℝ) * s⌋₊ + 1)).map Nat.castEmbedding) ?_]
  · rw [Finset.sum_map, Finset.mul_sum]
    refine Finset.sum_congr rfl fun n _ => ?_
    simp only [Nat.castEmbedding, Function.Embedding.coeFn_mk, Int.cast_natCast]
    have e : (φ (n / s) * (n : ℝ) ^ j : ℝ) = s ^ j * (((n : ℝ) / s) ^ j * φ (n / s)) := by
      rw [div_pow]
      have hNj : s ^ j ≠ 0 := pow_ne_zero _ hN0
      field_simp
    rw [e]; push_cast; ring
  · intro n hn
    have hz : φ ((n : ℝ) / s) = 0 := by
      rcases lt_or_ge n 0 with h0 | h0
      · apply h.eq_zero_of_nonpos
        apply div_nonpos_of_nonpos_of_nonneg _ hs0.le
        exact_mod_cast h0.le
      · obtain ⟨m, rfl⟩ := Int.eq_ofNat_of_zero_le h0
        have hm : ⌊(B : ℝ) * s⌋₊ < m := by
          by_contra hm; push Not at hm
          exact hn (Finset.mem_map.2 ⟨m, Finset.mem_range.2 (by omega), rfl⟩)
        have := h.vanish hs0 hm
        simpa using this
    simp [hz]

/-- Decay of `E_j`. -/
lemma Ej_decay {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) (j K : ℕ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ s : ℝ, 1 ≤ s → ∀ x : ℝ, |x| ≤ Real.pi →
      ‖Ej φ B s j x‖ * (1 + s * |x|) ^ K ≤ C * s ^ (j + 1) := by
  have hψ : ContDiff ℝ ∞ (fun y : ℝ => y ^ j * φ y) := (contDiff_id.pow j).mul h.smooth
  have hψc : HasCompactSupport (fun y : ℝ => y ^ j * φ y) := h.compact.mul_left
  obtain ⟨C, hC0, hC⟩ := esum_decay hψ hψc K
  refine ⟨C, hC0, fun s hs x hx => ?_⟩
  have hs0 : 0 < s := by linarith
  rw [Ej_eq_esum h hs j x, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hs0]
  have := hC s hs x hx
  calc s ^ j * ‖esum (fun n : ℤ => ((n : ℝ) / s) ^ j * φ ((n : ℝ) / s)) x‖ * (1 + s * |x|) ^ K
      = s ^ j * (‖esum (fun n : ℤ => ((n : ℝ) / s) ^ j * φ ((n : ℝ) / s)) x‖ * (1 + s * |x|) ^ K) := by ring
    _ ≤ s ^ j * (C * s) := mul_le_mul_of_nonneg_left this (by positivity)
    _ = C * s ^ (j + 1) := by ring

/-- `θ(t) ≥ 1/(2t + 1)` for `t ≥ 0`. -/
lemma liTheta_ge_inv {t : ℝ} (ht : 0 ≤ t) : 1 / (2 * t + 1) ≤ liTheta t := by
  rcases le_or_gt t (1 / 2) with h | h
  · have h1 : liTheta (1 / 2) ≤ liTheta t := liTheta_anti h
    have h2 : liTheta (1 / 2) = Real.pi / 2 := by
      unfold liTheta; rw [show 2 * (1 / 2 : ℝ) = 1 by norm_num, Real.arctan_one]; ring
    have h3 : 1 / (2 * t + 1) ≤ 1 := by rw [div_le_one (by positivity)]; linarith
    have := Real.pi_gt_three
    linarith
  · have h1 := liTheta_ge (show 0 < t by linarith)
    refine le_trans ?_ h1
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith

lemma liTheta_mem {t : ℝ} (ht : 0 ≤ t) : |liTheta t| ≤ Real.pi := by
  unfold liTheta
  have h1 : 0 ≤ arctan (2 * t) := arctan_nonneg.2 (by linarith)
  have h2 := arctan_lt_pi_div_two (2 * t)
  rw [abs_le]; constructor <;> linarith

/-- `‖E_j(θ(t))‖ ≤ C N^{j+1} (2t+1)⁴/(2t+1+N)⁴`. -/
lemma Ej_theta {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B) (j : ℕ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ N : ℝ, 1 ≤ N → ∀ t : ℝ, 0 ≤ t →
      ‖Ej φ B N j (liTheta t)‖ ≤ C * (N : ℝ) ^ (j + 1) * (2 * t + 1) ^ 4 / (2 * t + 1 + N) ^ 4 := by
  obtain ⟨C, hC0, hC⟩ := Ej_decay h j 4
  refine ⟨C, hC0, fun N hN t ht => ?_⟩
  have hN1 : (1 : ℝ) ≤ N := hN
  have h1 := hC N hN (liTheta t) (liTheta_mem ht)
  have hθ := liTheta_ge_inv ht
  have hθ0 : 0 ≤ liTheta t := le_trans (by positivity) hθ
  rw [abs_of_nonneg hθ0] at h1
  have hu : 0 < 2 * t + 1 := by linarith
  have hlow : (2 * t + 1 + N) / (2 * t + 1) ≤ 1 + N * liTheta t := by
    have : (N : ℝ) / (2 * t + 1) ≤ N * liTheta t := by
      rw [div_eq_mul_one_div]; exact mul_le_mul_of_nonneg_left hθ (by linarith)
    have e : (2 * t + 1 + N) / (2 * t + 1) = 1 + N / (2 * t + 1) := by field_simp
    linarith
  have hlow4 : ((2 * t + 1 + N) / (2 * t + 1)) ^ 4 ≤ (1 + N * liTheta t) ^ 4 :=
    pow_le_pow_left₀ (by positivity) hlow 4
  have hpos : 0 < ((2 * t + 1 + N) / (2 * t + 1)) ^ 4 := by positivity
  have h2 : ‖Ej φ B N j (liTheta t)‖ * ((2 * t + 1 + N) / (2 * t + 1)) ^ 4 ≤ C * (N : ℝ) ^ (j + 1) :=
    le_trans (mul_le_mul_of_nonneg_left hlow4 (norm_nonneg _)) h1
  rw [div_pow] at h2
  rw [le_div_iff₀ (by positivity)]
  have e : ‖Ej φ B N j (liTheta t)‖ * (2 * t + 1 + N) ^ 4 =
      ‖Ej φ B N j (liTheta t)‖ * ((2 * t + 1 + N) ^ 4 / (2 * t + 1) ^ 4) * (2 * t + 1) ^ 4 := by
    field_simp
  rw [e]
  exact mul_le_mul_of_nonneg_right h2 (by positivity)

end LeanLi
