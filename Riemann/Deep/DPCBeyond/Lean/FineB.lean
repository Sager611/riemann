import Riemann.Deep.DPCBeyond.Lean.LemmaM3
import Riemann.Deep.DPCBeyond.Lean.FinalB

/-!
# dpcbeyond: Theorem U′ and D(λ) > 0 for λ > 7/11 with Lemma M3 from its named inputs (fine version)

InputsFineB: InputsB with MomentM replaced by M3Inputs (Lemma Z + Hölder, M1, Minkowski, the bins of Off, M2);
inputsB_of_fine uses momentM_of_inputs (LemmaM3.lean). theoremUB_fine, D_pos_seven_elevenths_fine: the main results
from these finer inputs.
-/

namespace DPCBLean

open PPLean DPCLean

variable {ex lg : ℝ → ℝ}

/-- The named inputs for one Step-0 construction, with Lemma M3 from its own named inputs. -/
structure InputsFineB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ φ θ : ℝ) (D : DataB) : Prop where
  step3 : Step3 D.base
  main : MainTerm D.base
  step6 : Step6ClassesB Z β ex lg κ θ D.Hp
  step1 : Step1Pos D.base
  holder : HolderB ex lg D
  m3 : M3Inputs ex lg κ θ D
  step5c : Step5cB D
  cover : Covering Z φ D.base

/-- Step 0′ with the fine inputs. -/
def DPCInputsFineB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ : ℝ) : Prop :=
  ∀ φ θ : ℝ, 0 < φ → φ < 1 → 0 < θ → θ < min (κ / 2) 1 → 0 < φ * θ → φ * θ < 1 →
    1 < φ * (1 + θ) → φ * (1 + θ) < 2 → ∃ D : DataB, InputsFineB Z β ex lg κ φ θ D

theorem inputsB_of_fine (hEL : ExpLog ex lg) {Z : Zeros} {β : ℕ → ℝ} {κ : ℝ}
    (hin : DPCInputsFineB Z β ex lg κ) : DPCInputsB Z β ex lg κ := by
  intro φ θ h1 h2 h3 h4 h5 h6 h7 h8
  obtain ⟨D, inp⟩ := hin φ θ h1 h2 h3 h4 h5 h6 h7 h8
  exact ⟨D, ⟨inp.step3, inp.main, inp.step6, inp.step1, inp.holder,
    momentM_of_inputs hEL h3 inp.m3, inp.step5c, inp.cover⟩⟩

/-- **Theorem U′** from the fine inputs. -/
theorem theoremUB_fine {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    {κ : ℝ} (hκ : 0 < κ) (hden : DensityK Z β ex lg κ) (hin : DPCInputsFineB Z β ex lg κ) :
    TheoremUB Z κ :=
  theoremUB hZ hβ hEL htp hsc hκ hden (inputsB_of_fine hEL hin)

/-- **D(λ) > 0 for every λ > 7/11** from the fine inputs at every κ = 8/7 − δ. -/
theorem D_pos_seven_elevenths_fine {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    (hR : RvMUp Z) (hden : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DensityK Z β ex lg (8 / 7 - δ))
    (hin : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DPCInputsFineB Z β ex lg (8 / 7 - δ)) :
    ∀ lam : ℝ, 7 / 11 < lam → 0 < Dlam Z lam :=
  D_pos_seven_elevenths hZ hβ hEL htp hsc hR hden (fun δ h1 h2 => inputsB_of_fine hEL (hin δ h1 h2))

end DPCBLean

