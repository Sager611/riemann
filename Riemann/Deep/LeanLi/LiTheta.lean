module

public import Riemann.Deep.LeanLi.LiRH
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
public import Mathlib.Analysis.Real.Pi.Bounds

/-!
# LeanLi: the phase `θ(γ) = π - 2 arctan(2γ) = 2 arctan(1/(2γ))`
-/

@[expose] public section

open Real

namespace LeanLi

lemma liTheta_eq {γ : ℝ} (hγ : 0 < γ) : liTheta γ = 2 * arctan (1 / (2 * γ)) := by
  unfold liTheta
  rw [one_div, arctan_inv_of_pos (by positivity)]; ring

lemma liTheta_pos {γ : ℝ} (hγ : 0 < γ) : 0 < liTheta γ := by
  rw [liTheta_eq hγ]
  have : 0 < arctan (1 / (2 * γ)) := arctan_pos.2 (by positivity)
  linarith

lemma liTheta_lt_pi {γ : ℝ} (hγ : 0 < γ) : liTheta γ < π := by
  unfold liTheta
  have : 0 < arctan (2 * γ) := arctan_pos.2 (by positivity)
  linarith

lemma liTheta_le_inv {γ : ℝ} (hγ : 0 < γ) : liTheta γ ≤ 1 / γ := by
  rw [liTheta_eq hγ]
  have := arctan_le_self (x := 1 / (2 * γ)) (by positivity)
  have e : 2 * (1 / (2 * γ)) = 1 / γ := by field_simp
  linarith

lemma liTheta_anti {γ₁ γ₂ : ℝ} (h : γ₁ ≤ γ₂) : liTheta γ₂ ≤ liTheta γ₁ := by
  unfold liTheta
  have := arctan_mono (show 2 * γ₁ ≤ 2 * γ₂ by linarith)
  linarith

/-- `arctan y ≥ y/(1 + y²)` for `y ≥ 0`. -/
lemma arctan_ge_div {y : ℝ} (hy : 0 ≤ y) : y / (1 + y ^ 2) ≤ arctan y := by
  have hu : 0 ≤ arctan y := arctan_nonneg.2 hy
  have h1 : sin (arctan y) * cos (arctan y) = y / (1 + y ^ 2) := by
    rw [sin_arctan, cos_arctan]
    have hp : 0 < 1 + y ^ 2 := by positivity
    field_simp
    rw [Real.sq_sqrt hp.le]
  have h2 : sin (2 * arctan y) ≤ 2 * arctan y := sin_le (by positivity)
  rw [sin_two_mul] at h2
  linarith

lemma arctan_ge_half {y : ℝ} (hy : 0 ≤ y) (hy1 : y ≤ 1) : y / 2 ≤ arctan y := by
  have h := arctan_ge_div hy
  have : y / 2 ≤ y / (1 + y ^ 2) := div_le_div_of_nonneg_left hy (by positivity) (by nlinarith)
  linarith

/-- `θ(γ) ≥ 4γ/(4γ² + 1)`. -/
lemma liTheta_ge {γ : ℝ} (hγ : 0 < γ) : 4 * γ / (4 * γ ^ 2 + 1) ≤ liTheta γ := by
  rw [liTheta_eq hγ]
  have h := arctan_ge_div (y := 1 / (2 * γ)) (by positivity)
  have e : 2 * (1 / (2 * γ) / (1 + (1 / (2 * γ)) ^ 2)) = 4 * γ / (4 * γ ^ 2 + 1) := by
    field_simp; ring
  linarith

/-- `θ(γ) ≥ 1/√N` for `0 < γ ≤ √N/2`. -/
lemma liTheta_ge_inv_sqrt {γ : ℝ} {N : ℕ} (hN : 1 ≤ N) (hγ : 0 < γ) (hγN : γ ≤ √(N : ℝ) / 2) :
    1 / √(N : ℝ) ≤ liTheta γ := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hs : 0 < √(N : ℝ) := Real.sqrt_pos.2 hNpos
  have hs1 : 1 ≤ √(N : ℝ) := by rw [Real.one_le_sqrt]; exact_mod_cast hN
  rw [liTheta_eq hγ]
  have h1 : 1 / √(N : ℝ) ≤ 1 / (2 * γ) := one_div_le_one_div_of_le (by positivity) (by linarith)
  have h2 : arctan (1 / √(N : ℝ)) ≤ arctan (1 / (2 * γ)) := arctan_mono h1
  have h3 := arctan_ge_half (y := 1 / √(N : ℝ)) (by positivity) (by rw [div_le_one hs]; exact hs1)
  linarith

/-- `-log(2πθ(γ))/θ(γ)² ≤ γ² log(γ/2π) + log γ + 1` for `γ ≥ 2π`. -/
lemma neg_log_div_sq_le {γ : ℝ} (hγ : 2 * π ≤ γ) :
    -log (2 * π * liTheta γ) / liTheta γ ^ 2 ≤ γ ^ 2 * log (γ / (2 * π)) + log γ + 1 := by
  have hπ := pi_gt_three
  have hγ0 : 0 < γ := by linarith
  have hγ1 : 1 ≤ γ := by linarith
  set θ := liTheta γ with hθ
  have hθ0 : 0 < θ := liTheta_pos hγ0
  have hlow := liTheta_ge hγ0
  have hlg : 0 ≤ log (γ / (2 * π)) := log_nonneg (by rw [le_div_iff₀ (by positivity)]; linarith)
  have hlγ : 0 ≤ log γ := log_nonneg hγ1
  have hR : 0 ≤ γ ^ 2 * log (γ / (2 * π)) + log γ + 1 := by positivity
  rcases le_or_gt 1 (2 * π * θ) with h | h
  · have : -log (2 * π * θ) ≤ 0 := by linarith [log_nonneg h]
    have : -log (2 * π * θ) / θ ^ 2 ≤ 0 := div_nonpos_of_nonpos_of_nonneg this (sq_nonneg _)
    linarith
  · set u := γ + 1 / (4 * γ) with hu
    have hupos : 0 < u := by positivity
    have hθlow : 1 / u ≤ θ := by
      have e : 1 / u = 4 * γ / (4 * γ ^ 2 + 1) := by rw [hu]; field_simp
      rw [e]; exact hlow
    -- -log(2πθ) ≤ log(u/2π)
    have hA : -log (2 * π * θ) ≤ log (u / (2 * π)) := by
      have : 2 * π / u ≤ 2 * π * θ := by
        have := mul_le_mul_of_nonneg_left hθlow (by positivity : (0 : ℝ) ≤ 2 * π)
        rwa [mul_one_div] at this
      have h1 := log_le_log (by positivity) this
      have e : log (u / (2 * π)) = -log (2 * π / u) := by
        rw [← log_inv, inv_div]
      linarith
    have hA0 : 0 < -log (2 * π * θ) := by
      have := log_neg (by positivity) h; linarith
    -- 1/θ² ≤ u²
    have hB : 1 / θ ^ 2 ≤ u ^ 2 := by
      have h1 : 1 / θ ≤ u := by
        rw [div_le_iff₀ hθ0]; rw [div_le_iff₀ hupos] at hθlow; linarith
      have : (1 / θ) ^ 2 ≤ u ^ 2 := pow_le_pow_left₀ (by positivity) h1 2
      rwa [div_pow, one_pow] at this
    have hstep : -log (2 * π * θ) / θ ^ 2 ≤ log (u / (2 * π)) * u ^ 2 := by
      rw [div_eq_mul_one_div]
      exact mul_le_mul hA hB (by positivity) (le_trans hA0.le hA)
    -- log(u/2π) ≤ log(γ/2π) + 1/(4γ²), u² ≤ γ² + 1
    have hC : log (u / (2 * π)) ≤ log (γ / (2 * π)) + 1 / (4 * γ ^ 2) := by
      have e : u / (2 * π) = γ / (2 * π) * (1 + 1 / (4 * γ ^ 2)) := by rw [hu]; field_simp
      rw [e, log_mul (by positivity) (by positivity)]
      have := log_le_sub_one_of_pos (show 0 < 1 + 1 / (4 * γ ^ 2) by positivity)
      linarith
    have hD : u ^ 2 ≤ γ ^ 2 + 1 := by
      have e : u ^ 2 = γ ^ 2 + 1 / 2 + 1 / (16 * γ ^ 2) := by rw [hu]; field_simp; ring
      have : 1 / (16 * γ ^ 2) ≤ 1 / 2 := by
        rw [div_le_div_iff₀ (by positivity) (by norm_num)]; nlinarith
      linarith
    have hE : log (u / (2 * π)) * u ^ 2 ≤ (log (γ / (2 * π)) + 1 / (4 * γ ^ 2)) * (γ ^ 2 + 1) :=
      mul_le_mul hC hD (sq_nonneg _) (by positivity)
    have hF : (log (γ / (2 * π)) + 1 / (4 * γ ^ 2)) * (γ ^ 2 + 1) ≤
        γ ^ 2 * log (γ / (2 * π)) + log γ + 1 := by
      have e : (log (γ / (2 * π)) + 1 / (4 * γ ^ 2)) * (γ ^ 2 + 1) =
          γ ^ 2 * log (γ / (2 * π)) + log (γ / (2 * π)) + 1 / 4 + 1 / (4 * γ ^ 2) := by
        field_simp; ring
      have h1 : log (γ / (2 * π)) ≤ log γ := log_le_log (by positivity) (div_le_self hγ0.le (by linarith))
      have h2 : 1 / (4 * γ ^ 2) ≤ 1 / 4 := by
        rw [div_le_div_iff₀ (by positivity) (by norm_num)]; nlinarith
      linarith
    linarith

end LeanLi
