import Riemann.Deep.LeanNBRate.Bernstein
import Riemann.Deep.LeanNBRate.Nikolskii

/-!
# [N] Lemma B for one group

`p(u) = e^{icu} r(u)`, `r` with nodes in `[−w, w]`; `p′ = e^{icu}(icr + r′) = e^{icu}·(r.dmod c)`.
If `T₂ W² T₁² ≤ c²/2` (Turán `T₁ = (2κ)^{M−1}`, Nikolskii `T₂ = 4(2κ)^{2(M−1)}`, Bernstein `W`), then
`|r(u)|² ≤ 4T₁²T₂/(ℓc²) ∫_{ℓ/4}^{3ℓ/4} |p′|²` on `[0, ℓ]`.
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem ExpSum.dmod_eval_add (r : ExpSum) (c u : ℝ) :
    (r.dmod c).eval u = I * c * r.eval u + (r.dmod 0).eval u := by
  rw [ExpSum.dmod_eval_eq, ExpSum.dmod_eval_eq]
  unfold ExpSum.eval
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  push_cast; ring

/-- Turán between `[0, ℓ]` and `[ℓ/4, 3ℓ/4]`. -/
theorem turan_middle (hT : Nazarov1993_Turan) : ∃ κ : ℝ, 0 < κ ∧
    ∀ (r : ExpSum) (ℓ X : ℝ), 0 < ℓ → 1 ≤ r.M →
      (∀ u ∈ Icc (ℓ / 4) (3 * ℓ / 4), ‖r.eval u‖ ≤ X) →
      ∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ (2 * κ) ^ (r.M - 1) * X := by
  obtain ⟨κ, hκ, hT⟩ := turan_expSum hT
  refine ⟨κ, hκ, fun r ℓ X hℓ hM hX u hu => ?_⟩
  have hvol : (volume (Icc (ℓ / 4) (3 * ℓ / 4))).toReal = ℓ / 2 := by
    rw [Real.volume_Icc, ENNReal.toReal_ofReal (by linarith)]; ring
  have hpos : 0 < volume (Icc (ℓ / 4) (3 * ℓ / 4)) := by
    rw [Real.volume_Icc]; exact ENNReal.ofReal_pos.mpr (by linarith)
  have h := hT r 0 ℓ (Icc (ℓ / 4) (3 * ℓ / 4)) X hM hℓ measurableSet_Icc
    (Icc_subset_Icc (by linarith) (by linarith)) hpos hX u hu
  rw [hvol, show κ * (ℓ - 0) / (ℓ / 2) = 2 * κ by field_simp; ring] at h
  exact h

/-- **[N] Lemma B** (one group). -/
theorem lemmaB (hT : Nazarov1993_Turan) : ∃ κ : ℝ, 0 < κ ∧
    ∀ (r : ExpSum) (c w η ℓ : ℝ), 0 ≤ w → 0 < η → 0 < ℓ → 1 ≤ r.M → c ≠ 0 →
      (∀ j, r.node j ∈ Icc (-w) w) → 4 * A₀ / η * ((r.M : ℝ) + 1) ^ 2 ≤ ℓ →
      4 * (2 * κ) ^ (2 * (r.M - 1)) * ((2 * κ) ^ (r.M - 1) * CB * ((r.M : ℝ) + 1) ^ 3 *
        ((w + η) ^ 2 / η)) ^ 2 * ((2 * κ) ^ (r.M - 1)) ^ 2 ≤ c ^ 2 / 2 →
      ∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ^ 2 ≤
        4 * ((2 * κ) ^ (r.M - 1)) ^ 2 * (4 * (2 * κ) ^ (2 * (r.M - 1))) / (ℓ * c ^ 2) *
          ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖(r.dmod c).eval x‖ ^ 2 := by
  obtain ⟨κ₁, hκ₁, hTm⟩ := turan_middle hT
  obtain ⟨κ₂, hκ₂, hN⟩ := nikolskii hT
  obtain ⟨κ₃, hκ₃, hB⟩ := bernstein hT
  set κ := max κ₁ (max κ₂ κ₃)
  have hκ : 0 < κ := lt_max_of_lt_left hκ₁
  refine ⟨κ, hκ, fun r c w η ℓ hw hη hℓ hM hc hnodes hℓbig hsmall u hu => ?_⟩
  set M := r.M
  set T₁ := (2 * κ) ^ (M - 1)
  set T₂ := 4 * (2 * κ) ^ (2 * (M - 1))
  set W := (2 * κ) ^ (M - 1) * CB * ((M : ℝ) + 1) ^ 3 * ((w + η) ^ 2 / η)
  have hle1 : 2 * κ₁ ≤ 2 * κ := by linarith [le_max_left κ₁ (max κ₂ κ₃)]
  have hle2 : 2 * κ₂ ≤ 2 * κ := by
    linarith [le_max_left κ₂ κ₃, le_max_right κ₁ (max κ₂ κ₃)]
  have hle3 : 2 * κ₃ ≤ 2 * κ := by
    linarith [le_max_right κ₂ κ₃, le_max_right κ₁ (max κ₂ κ₃)]
  have hCB : 0 ≤ CB := by unfold CB; have := B₀_pos; have := A₀_pos; positivity
  refine lemmaB_abstract (r := r.eval) (rd := (r.dmod 0).eval) (q := (r.dmod c).eval)
    (T₁ := T₁) (T₂ := T₂) (W := W) hℓ (by positivity) (by positivity) (by positivity)
    r.continuous_eval (r.dmod 0).continuous_eval (r.dmod c).continuous_eval ?_ ?_ ?_ ?_ hsmall hc u hu
  · intro X hX u hu
    refine (hTm r ℓ X hℓ hM hX u hu).trans ?_
    have hX0 : 0 ≤ X := (norm_nonneg _).trans (hX (ℓ / 2) ⟨by linarith, by linarith⟩)
    exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (by linarith) hle1 _) hX0
  · intro u hu
    have h := hN r (ℓ / 4) (3 * ℓ / 4) hM (by linarith) u hu
    rw [show (3 * ℓ / 4 - ℓ / 4) = ℓ / 2 by ring] at h
    have hI : 0 ≤ ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖r.eval x‖ ^ 2 :=
      setIntegral_nonneg measurableSet_Icc fun x _ => by positivity
    have hp : (2 * κ₂) ^ (2 * (M - 1)) ≤ (2 * κ) ^ (2 * (M - 1)) := pow_le_pow_left₀ (by linarith) hle2 _
    calc ‖r.eval u‖ ^ 2 ≤ 2 * (2 * κ₂) ^ (2 * (M - 1)) *
          ((ℓ / 2)⁻¹ * ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖r.eval x‖ ^ 2) := h
      _ ≤ 2 * (2 * κ) ^ (2 * (M - 1)) *
          ((ℓ / 2)⁻¹ * ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖r.eval x‖ ^ 2) := by gcongr
      _ = T₂ / ℓ * ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖r.eval x‖ ^ 2 := by
          simp only [T₂]; field_simp; ring
  · intro X hX u hu
    have h := hB r w η ℓ X hw hη hℓ hM hnodes hℓbig hX u hu
    refine h.trans ?_
    have hX0 : 0 ≤ X := (norm_nonneg _).trans (hX 0 ⟨le_rfl, hℓ.le⟩)
    have hp : (2 * κ₃) ^ (M - 1) ≤ (2 * κ) ^ (M - 1) := pow_le_pow_left₀ (by linarith) hle3 _
    simp only [W]
    gcongr
  · intro u _; exact r.dmod_eval_add c u

end Riemann.Deep.LeanNBRate
