/-
E1PLean.Parity — H3′ (SL Thm 5.1 on (0, 1)) derived from its inputs (E1Lean.SL51Inputs: counts, parity product,
Wang's pair-sum limit, optimum; E1Lean.sl51), and Theorem E1′ with SL Thm 5.1 replaced by them.
-/
import Riemann.Deep.B3Cert.Lean.Window
import Riemann.Deep.SantiAudit.Lean.Parity

namespace E1PLean
open E1Lean

variable (F : Fns)

/-- H3′ follows from SL Thm 5.1's inputs at every θ ∈ (0, 1) (FULLPROOF §5: Wang's Thm local-pair for 0 < λ < θ < 1). -/
theorem H3P_of_inputs (I : Inputs) (hin : ∀ θ : ℝ, 0 < θ → θ < 1 → SL51Inputs F I θ) : H3P F I := by
  intro θ k hθ0 hθ1 _ hk1 hO
  obtain ⟨ι, Qf, C, hcount, hpar, hwang, hopt⟩ := hin θ hθ0 hθ1
  exact sl51 F I Qf C hθ0 hθ1.le hk1 hcount hpar hwang hopt hO

/-- **Theorem E1′** with SL Thm 5.1 replaced by its inputs. -/
theorem theoremE1P_parity (I : Inputs) (h1 : H1P I) (h2 : H2P F I) (hin : ∀ θ : ℝ, 0 < θ → θ < 1 → SL51Inputs F I θ)
    (h4 : H4P I) (θ : ℝ) (hθ1 : 371 / 750 < θ) (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₂ F I)) ∧ h F (371 / 750) (κ₂ F I) ≤ h F θ (κ₂ F I) ∧
      (0.0000335 : ℝ) < h F (371 / 750) (κ₂ F I) :=
  theoremE1P F I h1 h2 (H3P_of_inputs F I hin) h4 θ hθ1 hθ2

/-- **Theorem E1′**, finest inputs: Gaussian E1-M′, compact windows, SL Thm 1.3, SL Thm 5.1's inputs, Arb. -/
theorem theoremE1P_fine (I : Inputs) (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) (hG : H1G G) (hW : HCW I G)
    (h2 : H2P F I) (hin : ∀ θ : ℝ, 0 < θ → θ < 1 → SL51Inputs F I θ) (h4 : H4P I) (θ : ℝ) (hθ1 : 371 / 750 < θ)
    (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₂ F I)) ∧ h F (371 / 750) (κ₂ F I) ≤ h F θ (κ₂ F I) ∧
      (0.0000335 : ℝ) < h F (371 / 750) (κ₂ F I) :=
  theoremE1P F I (H1P_of_gauss I G hG hW) h2 (H3P_of_inputs F I hin) h4 θ hθ1 hθ2

end E1PLean
