import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Positivity

/-!
# dpcuncond (PROOF.md), new steps (b) and (c): admissible θ (Steps 0, 8) and the Hölder bookkeeping (Step 6e)

(b) step0_nonempty, step0_signs, step8_empty, step0_iff: for φ > 0 the interval (1/φ − 1, 1/2) is nonempty iff
φ > 2/3, and every θ in it satisfies DPC's sign condition 0 < φθ < 1 < φ(1 + θ) < 2 (DPC 2610.28527 l.744, 753);
step8_threshold: φ(1 + θ) > 1 iff φ > 1/(1 + θ) (the method's exact limit at length T^θ, Step 8).
(c) Eexp w r = E = (1 − 1/w)(1/2 − 1/(3r)); exp_mean: the σ-exponents of A₁, A₂, A₃ sum to 3E;
holder_bookkeeping: the three class bounds of 6e give 𝓗_{T,a} ≤ K T B g σ^E (σ^e abstract: any s with
s(x + y) = s x · s y); exists_params: explicit δ, w, r with θ < (1 − δ)E for every θ < 1/2.
-/

namespace DPCLean

/-! ## (b) Steps 0 and 8 -/

/-- Step 0: for φ > 2/3 the interval (1/φ − 1, 1/2) is nonempty. -/
theorem step0_nonempty {φ : ℝ} (h1 : 2 / 3 < φ) : ∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < 1 / 2 := by
  have hφ : 0 < φ := by linarith
  have h : 1 / φ < 3 / 2 := by
    rw [div_lt_iff₀ hφ]
    linarith
  exact ⟨(1 / φ - 1 + 1 / 2) / 2, by linarith, by linarith⟩

/-- Step 0: every θ ∈ (1/φ − 1, 1/2), φ ∈ (2/3, 1), satisfies 0 < φθ < 1 < φ(1 + θ) < 2 (DPC l.753). -/
theorem step0_signs {φ θ : ℝ} (h1 : 2 / 3 < φ) (h2 : φ < 1) (h3 : 1 / φ - 1 < θ) (h4 : θ < 1 / 2) :
    0 < φ * θ ∧ φ * θ < 1 ∧ 1 < φ * (1 + θ) ∧ φ * (1 + θ) < 2 := by
  have hφ : 0 < φ := by linarith
  have h5 : 1 < φ * (1 + θ) := by
    have h6 : 1 / φ < 1 + θ := by linarith
    rw [div_lt_iff₀ hφ] at h6
    linarith
  have hθ : 0 < θ := by
    have h7 : 1 < 1 / φ := by
      rw [lt_div_iff₀ hφ]
      linarith
    linarith
  refine ⟨mul_pos hφ hθ, ?_, h5, ?_⟩
  · nlinarith
  · nlinarith

/-- Step 8: the interval (1/φ − 1, 1/2) is empty for 0 < φ ≤ 2/3. -/
theorem step8_empty {φ : ℝ} (h0 : 0 < φ) (h : φ ≤ 2 / 3) : ¬ ∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < 1 / 2 := by
  rintro ⟨θ, h1, h2⟩
  have h3 : 3 / 2 ≤ 1 / φ := by
    rw [le_div_iff₀ h0]
    linarith
  linarith

theorem step0_iff {φ : ℝ} (h0 : 0 < φ) : (∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < 1 / 2) ↔ 2 / 3 < φ := by
  constructor
  · intro h
    by_contra hc
    exact step8_empty h0 (not_lt.1 hc) h
  · exact step0_nonempty

/-- Step 8: DPC's sign condition φ(1 + θ) > 1 is φ > 1/(1 + θ). -/
theorem step8_threshold {φ θ : ℝ} (hθ : 0 < 1 + θ) : 1 < φ * (1 + θ) ↔ 1 / (1 + θ) < φ :=
  (div_lt_iff₀ hθ).symm

/-- Step 8: at length T^θ with 0 ≤ θ < 1/2 the threshold 1/(1 + θ) lies above 2/3. -/
theorem step8_above {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ < 1 / 2) : 2 / 3 < 1 / (1 + θ) := by
  rw [lt_div_iff₀ (by linarith)]
  linarith

/-! ## (c) Step 6e: exponent arithmetic and the Hölder bookkeeping -/

/-- E = (1 − 1/w)(1/2 − 1/(3r)) (PROOF 6e). -/
noncomputable def Eexp (w r : ℝ) : ℝ := (1 - 1 / w) * (1 / 2 - 1 / (3 * r))

/-- 6e: the σ-exponents of A₁, A₂, A₃ sum to 3E. -/
theorem exp_mean (w r : ℝ) :
    (1 - 1 / w) / 4 + (1 - 1 / w) / 4 + (1 - 1 / w) * (1 - 1 / r) = Eexp w r + Eexp w r + Eexp w r := by
  unfold Eexp
  ring

/-- 6e, Hölder with exponents (3, 3, 3), as exponent bookkeeping on abstract nonnegative reals. With s e = σ^e
(any s with s(x + y) = s x · s y, s ≥ 0), B = (a + 2)^C, g = e^{θ(a+1)}:
A₁ ≤ K T B σ^{(1−1/w)/4}, A₂ ≤ K T B g³ σ^{(1−1/w)/4}, A₃ ≤ K T B σ^{(1−1/w)(1−1/r)} and 𝓗³ ≤ A₁A₂A₃
give 𝓗 ≤ K T B g σ^E. -/
theorem holder_bookkeeping {s : ℝ → ℝ} (hs : ∀ x y, s (x + y) = s x * s y) (hs0 : ∀ x, 0 ≤ s x)
    {w r K T B g A₁ A₂ A₃ H : ℝ} (hK : 0 ≤ K) (hT : 0 ≤ T) (hB : 0 ≤ B) (hg : 0 ≤ g)
    (hA₂ : 0 ≤ A₂) (hA₃ : 0 ≤ A₃)
    (h1 : A₁ ≤ K * T * B * s ((1 - 1 / w) / 4))
    (h2 : A₂ ≤ K * T * B * g ^ 3 * s ((1 - 1 / w) / 4))
    (h3 : A₃ ≤ K * T * B * s ((1 - 1 / w) * (1 - 1 / r)))
    (hH : H ^ 3 ≤ A₁ * A₂ * A₃) :
    H ≤ K * T * B * g * s (Eexp w r) := by
  have hsE := hs0 (Eexp w r)
  have hs1 := hs0 ((1 - 1 / w) / 4)
  have hs3 := hs0 ((1 - 1 / w) * (1 - 1 / r))
  have hX : 0 ≤ K * T * B * g * s (Eexp w r) := by positivity
  have hb1 : 0 ≤ K * T * B * s ((1 - 1 / w) / 4) := by positivity
  have hb2 : 0 ≤ K * T * B * g ^ 3 * s ((1 - 1 / w) / 4) := by positivity
  have hprod : A₁ * A₂ * A₃ ≤ (K * T * B * s ((1 - 1 / w) / 4)) *
      (K * T * B * g ^ 3 * s ((1 - 1 / w) / 4)) * (K * T * B * s ((1 - 1 / w) * (1 - 1 / r))) :=
    mul_le_mul (mul_le_mul h1 h2 hA₂ hb1) h3 hA₃ (mul_nonneg hb1 hb2)
  have hexp : s ((1 - 1 / w) / 4) * s ((1 - 1 / w) / 4) * s ((1 - 1 / w) * (1 - 1 / r))
      = s (Eexp w r) * s (Eexp w r) * s (Eexp w r) := by
    calc s ((1 - 1 / w) / 4) * s ((1 - 1 / w) / 4) * s ((1 - 1 / w) * (1 - 1 / r))
        = s ((1 - 1 / w) / 4 + (1 - 1 / w) / 4 + (1 - 1 / w) * (1 - 1 / r)) := by rw [hs, hs]
      _ = s (Eexp w r + Eexp w r + Eexp w r) := by rw [exp_mean]
      _ = s (Eexp w r) * s (Eexp w r) * s (Eexp w r) := by rw [hs, hs]
  have hcube : H ^ 3 ≤ (K * T * B * g * s (Eexp w r)) ^ 3 := by
    calc H ^ 3 ≤ A₁ * A₂ * A₃ := hH
      _ ≤ _ := hprod
      _ = (K * T * B) ^ 3 * g ^ 3 * (s ((1 - 1 / w) / 4) * s ((1 - 1 / w) / 4)
            * s ((1 - 1 / w) * (1 - 1 / r))) := by ring
      _ = (K * T * B) ^ 3 * g ^ 3 * (s (Eexp w r) * s (Eexp w r) * s (Eexp w r)) := by rw [hexp]
      _ = (K * T * B * g * s (Eexp w r)) ^ 3 := by ring
  exact le_of_pow_le_pow_left₀ (by norm_num) hX hcube

/-- 6e–6f: for every θ < 1/2, δ = x, w = 1/x, r = 1/(3x) with x = (1/2 − max θ 0)/4 give 0 < δ < 1, w, r > 1 and
θ < (1 − δ)E, i.e. κE > θ for Jutila's slope κ = 1 − δ. -/
theorem exists_params {θ : ℝ} (hθ : θ < 1 / 2) :
    ∃ δ w r : ℝ, 0 < δ ∧ δ < 1 ∧ 1 < w ∧ 1 < r ∧ θ < (1 - δ) * Eexp w r := by
  have ht0 : (0 : ℝ) ≤ max θ 0 := le_max_right _ _
  have htθ : θ ≤ max θ 0 := le_max_left _ _
  have ht1 : max θ 0 < 1 / 2 := max_lt hθ (by norm_num)
  have hx0 : 0 < (1 / 2 - max θ 0) / 4 := by linarith
  have hx1 : (1 / 2 - max θ 0) / 4 ≤ 1 / 8 := by linarith
  refine ⟨(1 / 2 - max θ 0) / 4, 1 / ((1 / 2 - max θ 0) / 4), 1 / (3 * ((1 / 2 - max θ 0) / 4)),
    hx0, by linarith, ?_, ?_, ?_⟩
  · rw [lt_div_iff₀ hx0]
    linarith
  · rw [lt_div_iff₀ (by linarith)]
    linarith
  · have hne : (1 / 2 - max θ 0) / 4 ≠ 0 := hx0.ne'
    have hE : Eexp (1 / ((1 / 2 - max θ 0) / 4)) (1 / (3 * ((1 / 2 - max θ 0) / 4)))
        = (1 - (1 / 2 - max θ 0) / 4) * (1 / 2 - (1 / 2 - max θ 0) / 4) := by
      unfold Eexp
      rw [one_div_one_div]
      congr 2
      field_simp
    rw [hE]
    nlinarith [mul_nonneg (sq_nonneg ((1 / 2 - max θ 0) / 4))
      (show (0 : ℝ) ≤ 5 / 2 - (1 / 2 - max θ 0) / 4 by linarith)]

end DPCLean

