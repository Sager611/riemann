/-
E1Lean.AlgebraC — Wang's c(θ) (E1.md §1, §4): x + cot x decreasing, c(θ) increasing and < 1 on (0, 1],
h(θ; k) increasing in θ, positivity ⟺ k > k*(θ), and the range bookkeeping of E1.md §1–§2.
-/
import Riemann.Deep.SantiAudit.Lean.AlgebraH

namespace E1Lean

variable (F : Fns)

theorem sin_pos' {x : ℝ} (hx : 0 < x) (hx1 : x ≤ 1) : 0 < F.sin x := by
  have h := F.sin_gt_sub_cube x hx
  have h2 : x * x ≤ 1 := by nlinarith
  have h3 : x ^ 3 = x * (x * x) := by ring
  have h4 : x * (x * x) ≤ x := by nlinarith
  linarith

/-- x + cot x decreases on (0, 1]: y − x < cot x − cot y (from sin_sub, sin_lt, sin_gt_sub_cube). -/
theorem cot_sub_gt {x y : ℝ} (hx : 0 < x) (hxy : x < y) (hy : y ≤ 1) : y - x < cot F x - cot F y := by
  have hy0 : 0 < y := by linarith
  have hsx := sin_pos' F hx (by linarith)
  have hsy := sin_pos' F hy0 hy
  have e : cot F x - cot F y = F.sin (y - x) / (F.sin x * F.sin y) := by
    rw [F.sin_sub]; unfold cot; field_simp
  rw [e, lt_div_iff₀ (mul_pos hsx hsy)]
  have hd : 0 < y - x := by linarith
  have hsd := F.sin_gt_sub_cube (y - x) hd
  have hxs := F.sin_lt x hx
  have hys := F.sin_lt y hy0
  have hp : F.sin x * F.sin y < x * y := by nlinarith
  have h1 : (y - x) * (F.sin x * F.sin y) < (y - x) * (x * y) := mul_lt_mul_of_pos_left hp hd
  have hx2 : x * x ≤ 1 := by nlinarith
  have hy2 : y * y ≤ 1 := by nlinarith
  have hxy1 : x * y ≤ 1 := by nlinarith
  have h3 : 0 ≤ 1 - (y - x) ^ 2 / 6 - x * y := by nlinarith
  have h2 : (y - x) * (x * y) ≤ (y - x) - (y - x) ^ 3 / 6 := by nlinarith [mul_nonneg hd.le h3]
  linarith

theorem sqrt2_facts : 0 < F.sqrt 2 ∧ F.sqrt 2 * F.sqrt 2 = 2 ∧ 1 ≤ F.sqrt 2 := by
  have h := F.mul_self_sqrt 2 (by norm_num)
  have h0 := F.sqrt_nonneg 2
  refine ⟨?_, h, ?_⟩
  · rcases h0.lt_or_eq with h1 | h1
    · exact h1
    · rw [← h1] at h; norm_num at h
  · nlinarith

theorem cW_lt_aux {s a b θ₁ θ₂ : ℝ} (hs : 0 < s) (hss : s * s = 2) (key : θ₂ / s - θ₁ / s < a - b) :
    2 - θ₁ / 2 - a / s < 2 - θ₂ / 2 - b / s := by
  have e : ∀ t : ℝ, t / s = t * s / 2 := fun t => by
    rw [div_eq_div_iff hs.ne' two_ne_zero, mul_assoc, hss]
  rw [e θ₁, e θ₂] at key
  rw [e a, e b]
  have hk := mul_lt_mul_of_pos_right key hs
  have e1 : (θ₂ * s / 2 - θ₁ * s / 2) * s = θ₂ - θ₁ := by
    have : (θ₂ * s / 2 - θ₁ * s / 2) * s = (θ₂ - θ₁) * (s * s) / 2 := by ring
    rw [this, hss]; ring
  rw [e1] at hk
  linarith

/-- c(θ) increases on (0, 1] (E1.md §4: c′(θ) = ½cot²(θ/√2) > 0). -/
theorem cW_lt {θ₁ θ₂ : ℝ} (h1 : 0 < θ₁) (h12 : θ₁ < θ₂) (h2 : θ₂ ≤ 1) : cW F θ₁ < cW F θ₂ := by
  obtain ⟨hs, hss, hs1⟩ := sqrt2_facts F
  have hx1 : 0 < θ₁ / F.sqrt 2 := div_pos h1 hs
  have hx12 : θ₁ / F.sqrt 2 < θ₂ / F.sqrt 2 := (div_lt_div_iff_of_pos_right hs).2 h12
  have hx2 : θ₂ / F.sqrt 2 ≤ 1 := by rw [div_le_iff₀ hs]; linarith
  exact cW_lt_aux hs hss (cot_sub_gt F hx1 hx12 hx2)

/-- c(θ) < 1 on (0, 1] (so 9 − k − 8c > 0 for k < 1: √ of a positive number in h). -/
theorem cW_lt_one {θ : ℝ} (h0 : 0 < θ) (h1 : θ ≤ 1) : cW F θ < 1 := by
  obtain ⟨hs, hss, hs1⟩ := sqrt2_facts F
  set x := θ / F.sqrt 2 with hxdef
  have hx : 0 < x := div_pos h0 hs
  have hx1 : x ≤ 1 := by rw [hxdef, div_le_iff₀ hs]; linarith
  have hsx : F.sqrt 2 * x = θ := by rw [hxdef]; field_simp
  have hxx : x * x * 2 = θ * θ := by
    have : (F.sqrt 2 * x) * (F.sqrt 2 * x) = θ * θ := by rw [hsx]
    nlinarith [this, hss]
  have hcb := F.cos_bound x (by rw [abs_of_pos hx]; exact hx1)
  rw [abs_of_pos hx] at hcb
  have hcl : 1 - x ^ 2 / 2 - x ^ 4 * (5 / 96) ≤ F.cos x := by linarith [(abs_le.1 hcb).1]
  have hsin := F.sin_lt x hx
  have hsin0 := sin_pos' F hx hx1
  have hx2 : x ^ 2 = θ ^ 2 / 2 := by nlinarith [hxx]
  have hx4 : x ^ 4 = θ ^ 4 / 4 := by rw [show x ^ 4 = (x ^ 2) ^ 2 by ring, hx2]; ring
  have t2 : θ ^ 2 ≤ 1 := by nlinarith
  have t4 : θ ^ 4 ≤ 1 := by nlinarith [mul_nonneg (sq_nonneg θ) (sub_nonneg.2 t2)]
  unfold cW
  rw [← hxdef]
  unfold cot
  rw [div_div]
  have hden : 0 < F.sin x * F.sqrt 2 := mul_pos hsin0 hs
  have key : (1 - θ / 2) * (F.sin x * F.sqrt 2) < F.cos x := by
    have h2 : F.sin x * F.sqrt 2 < θ := by nlinarith
    have h3 : 0 < 1 - θ / 2 := by linarith
    have h4 : (1 - θ / 2) * (F.sin x * F.sqrt 2) < (1 - θ / 2) * θ := mul_lt_mul_of_pos_left h2 h3
    have h5 : (1 - θ / 2) * θ ≤ 1 - x ^ 2 / 2 - x ^ 4 * (5 / 96) := by
      rw [hx2, hx4]; nlinarith [mul_nonneg (sub_nonneg.2 h1) (show (0 : ℝ) ≤ 3 - θ by linarith)]
    linarith
  have := (lt_div_iff₀ hden).2 key
  linarith

/-- h(θ; k) increases in θ on (0, 1] (k < 1): E1.md §4, SL §5. -/
theorem h_mono_theta {θ₁ θ₂ k : ℝ} (h1 : 0 < θ₁) (h12 : θ₁ ≤ θ₂) (h2 : θ₂ ≤ 1) (hk : k < 1) :
    h F θ₁ k ≤ h F θ₂ k := by
  unfold h
  rcases h12.lt_or_eq with hlt | heq
  · exact hc_mono_c F hk (cW_lt F h1 hlt h2).le (cW_lt_one F (by linarith) h2)
  · rw [heq]

/-- E1.md §4: for θ ∈ (0, 1], k < 1: h(θ; k) > 0 ⟺ k > k*(θ). -/
theorem h_pos_iff {θ k : ℝ} (h0 : 0 < θ) (h1 : θ ≤ 1) (hk : k < 1) : 0 < h F θ k ↔ kstar F θ < k :=
  hc_pos_iff F hk (cW_lt_one F h0 h1)

/-- E1.md §1: ν = 3/25 < (3θ − 1)/4 ⟺ θ > 37/75. -/
theorem nu_range_iff (θ : ℝ) : (3 / 25 : ℝ) < (3 * θ - 1) / 4 ↔ 37 / 75 < θ := by
  constructor <;> intro h <;> linarith

/-- E1.md §1: ν = 1249/10000 < (3θ − 1)/4 ⟺ θ > 3749/7500. -/
theorem nu2_range_iff (θ : ℝ) : (1249 / 10000 : ℝ) < (3 * θ - 1) / 4 ↔ 3749 / 7500 < θ := by
  constructor <;> intro h <;> linarith

/-- E1.md §2 Step 5: relative error T^{1/2 + 2ν − 3θ/2}, power-small ⟺ ν < (3θ − 1)/4. -/
theorem exponent_iff (θ ν : ℝ) : 1 / 2 + 2 * ν - 3 * θ / 2 < 0 ↔ ν < (3 * θ - 1) / 4 := by
  constructor <;> intro h <;> linarith

/-- E1.md §1: (3θ − 1)/4 < 1/2 for θ < 1. -/
theorem range_lt_half {θ : ℝ} (h : θ < 1) : (3 * θ - 1) / 4 < 1 / 2 := by linarith

end E1Lean
