import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# Theorem 1′(a) (RCApps2 NOTE Prop 1, l.28–61): the exponent bookkeeping

e(σ) = 2(1−σ)/(3−2σ); W = 𝔔T, 𝔔 = 9|D| (= |D|(M₀+3)², M₀ = 0), V = 𝔔^{1/2}, Y = W^{1/(3−2σ)},
Z = 2Y log²Y (RCApps2 l.35–36).
-/

open Real

namespace RingClass

/-- e(σ) = 2(1 − σ)/(3 − 2σ) (RCApps2 l.20) -/
noncomputable def ee (σ : ℝ) : ℝ := 2 * (1 - σ) / (3 - 2 * σ)

lemma ee_nonneg (σ : ℝ) (h1 : σ ≤ 1) (_h0 : 0.6 ≤ σ) : 0 ≤ ee σ := by
  unfold ee; apply div_nonneg <;> linarith

lemma ee_le (σ : ℝ) (h1 : σ ≤ 1) (_h0 : 0.6 ≤ σ) : ee σ ≤ 1 := by
  unfold ee; rw [div_le_one (by linarith)]; linarith

/-- the balance (RCApps2 l.59): Y^{1/2−σ}(𝔔T)^{1/2} = (𝔔T)^{e(σ)}, in exponents -/
lemma balance (σ : ℝ) (h1 : σ ≤ 1) : 1 / (3 - 2 * σ) * (1 / 2 - σ) + 1 / 2 = ee σ := by
  have h : (3 - 2 * σ) ≠ 0 := by linarith
  unfold ee; field_simp; ring

/-- V = 𝔔^{1/2} ≤ Y = W^{1/(3−2σ)} (RCApps2 l.35, "as 3 − 2σ ≤ 2") and V ≥ 2 -/
lemma V_le_Y (Q T σ : ℝ) (hQ : 4 ≤ Q) (hT : 1 ≤ T) (h0 : 0.6 ≤ σ) (h1 : σ ≤ 1) :
    2 ≤ Q ^ (1 / 2 : ℝ) ∧ Q ^ (1 / 2 : ℝ) ≤ (Q * T) ^ (1 / (3 - 2 * σ)) := by
  have hy : 1 / 2 ≤ 1 / (3 - 2 * σ) := one_div_le_one_div_of_le (by linarith) (by linarith)
  have hQT : Q ≤ Q * T := le_mul_of_one_le_right (by linarith) hT
  constructor
  · have h4 : (4 : ℝ) ^ (1 / 2 : ℝ) = 2 := by
      have e : (4 : ℝ) = 2 ^ (2 : ℝ) := by rw [Real.rpow_two]; norm_num
      have e2 : (2 : ℝ) * (1 / 2) = 1 := by norm_num
      rw [e, ← Real.rpow_mul (by norm_num : (0 : ℝ) ≤ 2), e2, Real.rpow_one]
    calc (2 : ℝ) = (4 : ℝ) ^ (1 / 2 : ℝ) := h4.symm
      _ ≤ Q ^ (1 / 2 : ℝ) := Real.rpow_le_rpow (by norm_num) hQ (by norm_num)
  · calc Q ^ (1 / 2 : ℝ) ≤ (Q * T) ^ (1 / 2 : ℝ) := Real.rpow_le_rpow (by linarith) hQT (by norm_num : (0 : ℝ) ≤ 1 / 2)
      _ ≤ (Q * T) ^ (1 / (3 - 2 * σ)) :=
          Real.rpow_le_rpow_of_exponent_le (one_le_mul_of_one_le_of_one_le (by linarith) hT) hy

end RingClass
