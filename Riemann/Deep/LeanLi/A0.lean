module

public import Riemann.Deep.LeanLi.IClosed
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
public import Mathlib.Analysis.Real.Pi.Bounds

/-!
# LeanLi: bounds on `ℓ(v) = log cot(v/2) - log 4π` and `a₀(v) = ℓ(v)/(8π sin²(v/2))`
-/

@[expose] public section

open Real Set

namespace LeanLi

lemma sin_half_ge {v : ℝ} (h0 : 0 ≤ v) (h1 : v ≤ π) : v / π ≤ sin (v / 2) := by
  have := mul_le_sin (x := v / 2) (by linarith) (by linarith)
  have e : 2 / π * (v / 2) = v / π := by field_simp
  linarith

lemma cos_half_ge {v : ℝ} (h0 : 0 ≤ v) (h1 : v ≤ π) : (π - v) / π ≤ cos (v / 2) := by
  have := mul_le_sin (x := π / 2 - v / 2) (by linarith) (by linarith [pi_pos])
  rw [sin_pi_div_two_sub] at this
  have e : 2 / π * (π / 2 - v / 2) = (π - v) / π := by field_simp
  linarith

lemma sin_half_pos {v : ℝ} (h0 : 0 < v) (h1 : v ≤ π) : 0 < sin (v / 2) :=
  sin_pos_of_pos_of_lt_pi (by linarith) (by linarith [pi_pos])

lemma cos_half_pos {v : ℝ} (h0 : 0 ≤ v) (h1 : v < π) : 0 < cos (v / 2) :=
  cos_pos_of_mem_Ioo ⟨by linarith [pi_pos], by linarith⟩

lemma neg_log_le_of_ge {x y : ℝ} (hy : 0 < y) (hyx : y ≤ x) : -log x ≤ log (1 / y) := by
  rw [one_div, log_inv]; linarith [log_le_log hy hyx]

lemma abs_log_le_of_ge {x y : ℝ} (hy : 0 < y) (hyx : y ≤ x) (hx1 : x ≤ 1) : |log x| ≤ log (1 / y) := by
  have hx : 0 < x := lt_of_lt_of_le hy hyx
  rw [abs_of_nonpos (log_nonpos hx.le hx1)]
  exact neg_log_le_of_ge hy hyx

/-- `|ℓ(v)| ≤ log(π/v) + log(π/(π - v)) + log 4π`. -/
lemma abs_ell_le {v : ℝ} (h0 : 0 < v) (h1 : v < π) :
    |ell v| ≤ log (π / v) + log (π / (π - v)) + log (4 * π) := by
  have hπ := pi_pos
  have hs := abs_log_le_of_ge (div_pos h0 hπ) (sin_half_ge h0.le h1.le) (sin_le_one _)
  have hc := abs_log_le_of_ge (div_pos (by linarith) hπ) (cos_half_ge h0.le h1.le) (cos_le_one _)
  rw [one_div_div] at hs hc
  have h4 : 0 < log (4 * π) := log_pos (by nlinarith [pi_gt_three])
  unfold ell
  calc |log (cos (v / 2)) - log (sin (v / 2)) - log (4 * π)|
      ≤ |log (cos (v / 2))| + |log (sin (v / 2))| + |log (4 * π)| := by
        have := abs_sub (log (cos (v / 2)) - log (sin (v / 2))) (log (4 * π))
        have := abs_sub (log (cos (v / 2))) (log (sin (v / 2)))
        linarith
    _ ≤ log (π / (π - v)) + log (π / v) + log (4 * π) := by
        rw [abs_of_pos h4]; linarith
    _ = _ := by ring

/-- `|a₀(v)| ≤ π|ℓ(v)|/(8v²)`. -/
lemma abs_a0_le {v : ℝ} (h0 : 0 < v) (h1 : v ≤ π) : |a0 v| ≤ π * |ell v| / (8 * v ^ 2) := by
  have hπ := pi_pos
  have hs := sin_half_ge h0.le h1
  have hsp := sin_half_pos h0 h1
  unfold a0
  rw [abs_div, abs_of_pos (by positivity : 0 < 8 * π * sin (v / 2) ^ 2)]
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  have h2 : (v / π) ^ 2 ≤ sin (v / 2) ^ 2 := pow_le_pow_left₀ (by positivity) hs 2
  have h3 : v ^ 2 ≤ π ^ 2 * sin (v / 2) ^ 2 := by
    rw [div_pow, div_le_iff₀ (by positivity)] at h2; linarith
  have h4 : 0 ≤ |ell v| := abs_nonneg _
  nlinarith [mul_le_mul_of_nonneg_left h3 (mul_nonneg h4 (by positivity : (0:ℝ) ≤ 8))]

/-- Derivative of `ℓ`: `ℓ'(v) = -sin(v/2)/(2cos(v/2)) - cos(v/2)/(2 sin(v/2))`. -/
noncomputable def ell' (v : ℝ) : ℝ := -(sin (v / 2) / (2 * cos (v / 2))) - cos (v / 2) / (2 * sin (v / 2))

lemma hasDerivAt_ell {v : ℝ} (h0 : 0 < v) (h1 : v < π) : HasDerivAt ell (ell' v) v := by
  have hs := (sin_half_pos h0 h1.le).ne'
  have hc := (cos_half_pos h0.le h1).ne'
  have hd2 : HasDerivAt (fun v : ℝ => v / 2) (1 / 2) v := by
    simpa using (hasDerivAt_id v).div_const 2
  have hcos : HasDerivAt (fun v : ℝ => cos (v / 2)) (-sin (v / 2) * (1 / 2)) v := hd2.cos
  have hsin : HasDerivAt (fun v : ℝ => sin (v / 2)) (cos (v / 2) * (1 / 2)) v := hd2.sin
  have h := ((hcos.log hc).sub (hsin.log hs)).sub_const (log (4 * π))
  unfold ell ell'
  convert h using 1
  field_simp
  try ring

/-- Derivative of `a₀`. -/
noncomputable def a0' (v : ℝ) : ℝ :=
  (ell' v * sin (v / 2) - ell v * cos (v / 2)) / (8 * π * sin (v / 2) ^ 3)

lemma hasDerivAt_a0 {v : ℝ} (h0 : 0 < v) (h1 : v < π) : HasDerivAt a0 (a0' v) v := by
  have hs := (sin_half_pos h0 h1.le).ne'
  have hd2 : HasDerivAt (fun v : ℝ => v / 2) (1 / 2) v := by
    simpa using (hasDerivAt_id v).div_const 2
  have hsin : HasDerivAt (fun v : ℝ => sin (v / 2)) (cos (v / 2) * (1 / 2)) v := hd2.sin
  have hden : HasDerivAt (fun v : ℝ => 8 * π * sin (v / 2) ^ 2)
      (8 * π * (2 * sin (v / 2) * (cos (v / 2) * (1 / 2)))) v := by
    have := (hsin.pow 2).const_mul (8 * π)
    convert this using 1; push_cast; ring
  have hden0 : 8 * π * sin (v / 2) ^ 2 ≠ 0 := by positivity
  have h := (hasDerivAt_ell h0 h1).div hden hden0
  unfold a0 a0'
  convert h using 1
  field_simp
  try ring

/-- `|a₀'(v)| ≤ (1/(2c₀) + 1/2 + |ℓ(v)|)/(8π s₀³)` when `sin(v/2) ≥ s₀ > 0`, `cos(v/2) ≥ c₀ > 0`. -/
lemma abs_a0'_le {v s₀ c₀ : ℝ} (hs₀ : 0 < s₀) (hc₀ : 0 < c₀) (hs : s₀ ≤ sin (v / 2))
    (hc : c₀ ≤ cos (v / 2)) :
    |a0' v| ≤ (1 / (2 * c₀) + 1 / 2 + |ell v|) / (8 * π * s₀ ^ 3) := by
  have hπ := pi_pos
  have hsp : 0 < sin (v / 2) := lt_of_lt_of_le hs₀ hs
  have hcp : 0 < cos (v / 2) := lt_of_lt_of_le hc₀ hc
  have e : ell' v * sin (v / 2) - ell v * cos (v / 2) =
      -(sin (v / 2) ^ 2 / (2 * cos (v / 2))) - cos (v / 2) / 2 - ell v * cos (v / 2) := by
    unfold ell'; field_simp; try ring
  have hnum : |ell' v * sin (v / 2) - ell v * cos (v / 2)| ≤ 1 / (2 * c₀) + 1 / 2 + |ell v| := by
    rw [e]
    have hs1 : sin (v / 2) ^ 2 ≤ 1 := by nlinarith [sin_sq_add_cos_sq (v / 2)]
    have hc1 := cos_le_one (v / 2)
    have t1 : sin (v / 2) ^ 2 / (2 * cos (v / 2)) ≤ 1 / (2 * c₀) := by
      rw [div_le_div_iff₀ (by positivity) (by positivity)]; nlinarith
    have t1' : 0 ≤ sin (v / 2) ^ 2 / (2 * cos (v / 2)) := by positivity
    have t3 : |ell v * cos (v / 2)| ≤ |ell v| := by
      rw [abs_mul, abs_of_pos hcp]
      exact mul_le_of_le_one_right (abs_nonneg _) hc1
    calc |-(sin (v / 2) ^ 2 / (2 * cos (v / 2))) - cos (v / 2) / 2 - ell v * cos (v / 2)|
        ≤ |-(sin (v / 2) ^ 2 / (2 * cos (v / 2))) - cos (v / 2) / 2| + |ell v * cos (v / 2)| :=
          abs_sub _ _
      _ ≤ (sin (v / 2) ^ 2 / (2 * cos (v / 2)) + cos (v / 2) / 2) + |ell v| := by
          gcongr
          rw [show -(sin (v / 2) ^ 2 / (2 * cos (v / 2))) - cos (v / 2) / 2 =
            -(sin (v / 2) ^ 2 / (2 * cos (v / 2)) + cos (v / 2) / 2) by ring, abs_neg]
          exact le_of_eq (abs_of_nonneg (by positivity))
      _ ≤ 1 / (2 * c₀) + 1 / 2 + |ell v| := by linarith
  unfold a0'
  rw [abs_div, abs_of_pos (by positivity : 0 < 8 * π * sin (v / 2) ^ 3)]
  have hden : 8 * π * s₀ ^ 3 ≤ 8 * π * sin (v / 2) ^ 3 := by
    have := pow_le_pow_left₀ hs₀.le hs 3
    nlinarith
  calc |ell' v * sin (v / 2) - ell v * cos (v / 2)| / (8 * π * sin (v / 2) ^ 3)
      ≤ (1 / (2 * c₀) + 1 / 2 + |ell v|) / (8 * π * sin (v / 2) ^ 3) :=
        div_le_div_of_nonneg_right hnum (by positivity)
    _ ≤ (1 / (2 * c₀) + 1 / 2 + |ell v|) / (8 * π * s₀ ^ 3) :=
        div_le_div_of_nonneg_left (by positivity) (by positivity) hden

/-- **`2πa₀(a) ≤ -log(2πa)/a² + O(1 + |log 2πa|)`**:
`ℓ(a)/(4 sin²(a/2)) ≤ -log(2πa)/a² + 3 + |log(2πa)|` for `0 < a ≤ π`. -/
lemma ell_div_le {a : ℝ} (h0 : 0 < a) (h1 : a ≤ π) :
    ell a / (4 * sin (a / 2) ^ 2) ≤ -log (2 * π * a) / a ^ 2 + 3 + |log (2 * π * a)| := by
  have hπ := pi_pos
  have hsp := sin_half_pos h0 h1
  have hlc : log (cos (a / 2)) ≤ 0 := by
    rcases eq_or_lt_of_le (cos_nonneg_of_mem_Icc ⟨by linarith, by linarith⟩ : 0 ≤ cos (a / 2))
      with h | h
    · rw [← h, log_zero]
    · exact log_nonpos h.le (cos_le_one _)
  have hell : ell a ≤ -log (sin (a / 2)) - log (4 * π) := by unfold ell; linarith
  have hS0 : 0 < 4 * sin (a / 2) ^ 2 := by positivity
  rcases le_or_gt a 1 with ha1 | ha1
  · have hs1 : a / 2 * (1 - a ^ 2 / 24) ≤ sin (a / 2) := by
      have := sin_gt_sub_cube (x := a / 2) (by linarith); nlinarith
    have hs2 : sin (a / 2) ≤ a / 2 := sin_le (by linarith)
    have ha2 : a ^ 2 ≤ 1 := by nlinarith
    have hy : 0 < 1 - a ^ 2 / 24 := by nlinarith
    have hlog : -log (sin (a / 2)) ≤ log (2 / a) + a ^ 2 / 12 := by
      have h1' : log (a / 2 * (1 - a ^ 2 / 24)) ≤ log (sin (a / 2)) := log_le_log (by positivity) hs1
      rw [log_mul (by positivity) hy.ne'] at h1'
      have h2' := one_sub_inv_le_log_of_pos hy
      have h3' : -(a ^ 2 / 12) ≤ 1 - (1 - a ^ 2 / 24)⁻¹ := by
        have e : 1 - (1 - a ^ 2 / 24)⁻¹ = -(a ^ 2 / 24) / (1 - a ^ 2 / 24) := by
          rw [eq_div_iff hy.ne', sub_mul, inv_mul_cancel₀ hy.ne']; ring
        rw [e, le_div_iff₀ hy]; nlinarith
      have h4' : log (a / 2) = -log (2 / a) := by
        rw [← log_inv]; congr 1; field_simp
      linarith
    set X := -log (2 * π * a) + a ^ 2 / 12 with hXdef
    have hX : ell a ≤ X := by
      have e : log (2 / a) - log (4 * π) = -log (2 * π * a) := by
        rw [← log_div (by positivity) (by positivity), ← log_inv]; congr 1; field_simp; ring
      rw [hXdef]; linarith
    have hS : a ^ 2 * (1 - a ^ 2 / 12) ≤ 4 * sin (a / 2) ^ 2 := by
      have := pow_le_pow_left₀ (by positivity) hs1 2
      nlinarith
    have hS' : 4 * sin (a / 2) ^ 2 ≤ a ^ 2 := by nlinarith [sin_pos_of_pos_of_lt_pi (x := a / 2)]
    have hy12 : 0 < 1 - a ^ 2 / 12 := by nlinarith
    set Dd := 1 / (4 * sin (a / 2) ^ 2) - 1 / a ^ 2 with hDd
    have hD0 : 0 ≤ Dd := sub_nonneg.2 (one_div_le_one_div_of_le (by positivity) hS')
    have hD1 : Dd ≤ 1 / 11 := by
      have t : 1 / (4 * sin (a / 2) ^ 2) ≤ 1 / (a ^ 2 * (1 - a ^ 2 / 12)) :=
        one_div_le_one_div_of_le (by positivity) hS
      have h12 : (12 - a ^ 2) ≠ 0 := by nlinarith
      have e : 1 / (a ^ 2 * (1 - a ^ 2 / 12)) - 1 / a ^ 2 = 1 / (12 - a ^ 2) := by
        rw [show a ^ 2 * (1 - a ^ 2 / 12) = a ^ 2 * (12 - a ^ 2) / 12 by ring]
        field_simp
        ring
      have t2 : 1 / (12 - a ^ 2) ≤ 1 / 11 := one_div_le_one_div_of_le (by norm_num) (by linarith)
      rw [hDd]; linarith
    have hsplit : X / (4 * sin (a / 2) ^ 2) = X / a ^ 2 + X * Dd := by
      rw [hDd]; field_simp; ring
    have hXD : X * Dd ≤ |X| / 11 := by
      calc X * Dd ≤ |X| * Dd := mul_le_mul_of_nonneg_right (le_abs_self X) hD0
        _ ≤ |X| * (1 / 11) := mul_le_mul_of_nonneg_left hD1 (abs_nonneg _)
        _ = |X| / 11 := by ring
    have hXa : X / a ^ 2 = -log (2 * π * a) / a ^ 2 + 1 / 12 := by
      rw [hXdef]; field_simp
    have hXabs : |X| ≤ |log (2 * π * a)| + 1 / 12 := by
      rw [hXdef]
      calc |-log (2 * π * a) + a ^ 2 / 12| ≤ |-log (2 * π * a)| + |a ^ 2 / 12| := abs_add_le _ _
        _ ≤ |log (2 * π * a)| + 1 / 12 := by
          rw [abs_neg, abs_of_nonneg (by positivity : (0 : ℝ) ≤ a ^ 2 / 12)]; linarith
    calc ell a / (4 * sin (a / 2) ^ 2) ≤ X / (4 * sin (a / 2) ^ 2) :=
          div_le_div_of_nonneg_right hX hS0.le
      _ = X / a ^ 2 + X * Dd := hsplit
      _ ≤ -log (2 * π * a) / a ^ 2 + 1 / 12 + |X| / 11 := by rw [hXa]; linarith
      _ ≤ -log (2 * π * a) / a ^ 2 + 3 + |log (2 * π * a)| := by
          have := abs_nonneg (log (2 * π * a)); linarith
  · have hsa : 1 / π ≤ sin (a / 2) := by
      have := sin_half_ge h0.le h1
      have : 1 / π ≤ a / π := div_le_div_of_nonneg_right ha1.le hπ.le
      linarith
    have hl : -log (sin (a / 2)) ≤ log π := by
      have := neg_log_le_of_ge (by positivity) hsa
      rwa [one_div_one_div] at this
    have h4 : log (4 * π) = log 4 + log π := log_mul (by norm_num) hπ.ne'
    have hl4 : 0 < log 4 := log_pos (by norm_num)
    have hneg : ell a ≤ 0 := by linarith
    have hlhs : ell a / (4 * sin (a / 2) ^ 2) ≤ 0 := div_nonpos_of_nonpos_of_nonneg hneg hS0.le
    have hrhs : 0 ≤ -log (2 * π * a) / a ^ 2 + 3 + |log (2 * π * a)| := by
      have ha2 : 1 < a ^ 2 := by nlinarith
      have t : -log (2 * π * a) / a ^ 2 ≥ -|log (2 * π * a)| := by
        rw [ge_iff_le, le_div_iff₀ (by positivity)]
        have := le_abs_self (log (2 * π * a))
        have := mul_nonneg (abs_nonneg (log (2 * π * a))) (by linarith : (0 : ℝ) ≤ a ^ 2 - 1)
        nlinarith
      linarith
    linarith

end LeanLi
