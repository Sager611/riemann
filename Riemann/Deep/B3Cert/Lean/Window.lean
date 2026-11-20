/-
E1PLean.Window — the compact-window step of Theorem E1′ (FULLPROOF.md §5 l. 106–110; SL TeX l. 730–758) as a named
hypothesis, and its arithmetic: from the Gaussian moment at width T^θ′ (FULLPROOF §§2–4) to SL Thm 1.2 at window T^θ.
-/
import Riemann.Deep.B3Cert.Lean.Arch

namespace E1PLean
open E1Lean

variable (I : Inputs)

/-- (H1G) E1-M′ in its Gaussian form (FULLPROOF §§2–4: Lemma B, Lemma R, branches; §5 l. 108–109): for θ′ ∈ (1/3, 1) and
  0 < ν < (3θ′ − 1)/4, the Gaussian-weighted moment G θ′ ν R P Q (Gaussian of width T^θ′, uniformly in the centre
  T′ ≍ T) holds for every admissible (R, P, Q). G is opaque; its meaning is this docstring. -/
def H1G (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) : Prop :=
  ∀ θ ν : ℝ, 1 / 3 < θ → θ < 1 → 0 < ν → ν < (3 * θ - 1) / 4 → ∀ R P Q, Adm R P Q → G θ ν R P Q

/-- (HCW) SL TeX l. 730–758, compact windows (profile-free; FULLPROOF §5 l. 106–110): the Gaussian moment at width
  T^(θ − η), 0 < η < θ, gives SL Thm 1.2's conclusion at window T^θ (w_D = φ_{s0} * Σ_{j<D} …, D fixed). -/
def HCW (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) : Prop :=
  ∀ θ η ν : ℝ, 0 < η → η < θ → θ < 1 → 0 < ν → ∀ R P Q, Adm R P Q → G (θ - η) ν R P Q → I.Moment θ ν R P Q

/-- FULLPROOF §5's choice at ν = 121/1000: η := (θ − 371/750)/2 > 0, θ′ := θ − η > 371/750, so 121/1000 < (3θ′ − 1)/4. -/
theorem eta_edge {θ : ℝ} (hθ : 371 / 750 < θ) :
    0 < (θ - 371 / 750) / 2 ∧ (θ - 371 / 750) / 2 < θ ∧ 371 / 750 < θ - (θ - 371 / 750) / 2 ∧
      (121 / 1000 : ℝ) < (3 * (θ - (θ - 371 / 750) / 2) - 1) / 4 := by
  refine ⟨by linarith, by linarith, by linarith, ?_⟩
  rw [nu_edge_iff]
  linarith

/-- The same step for every ν > 0: η := (θ − (4ν + 1)/3)/2. -/
theorem eta_step {θ ν : ℝ} (hν0 : 0 < ν) (hν : ν < (3 * θ - 1) / 4) :
    0 < (θ - (4 * ν + 1) / 3) / 2 ∧ (θ - (4 * ν + 1) / 3) / 2 < θ ∧
      1 / 3 < θ - (θ - (4 * ν + 1) / 3) / 2 ∧ ν < (3 * (θ - (θ - (4 * ν + 1) / 3) / 2) - 1) / 4 := by
  refine ⟨by linarith, by linarith, by linarith, by linarith⟩

/-- H1G + HCW ⟹ H1′ (E1-M′ in SL Thm 1.2's form). -/
theorem H1P_of_gauss (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) (hG : H1G G) (hW : HCW I G) : H1P I := by
  intro θ ν hθ1 hθ2 hν0 hν R P Q hadm
  obtain ⟨e1, e2, e3, e4⟩ := eta_step hν0 hν
  exact hW θ _ ν e1 e2 hθ2 hν0 R P Q hadm (hG _ ν e3 (by linarith) hν0 e4 R P Q hadm)

/-- **Theorem E1′** from the Gaussian form of E1-M′ and the compact-window step. -/
theorem theoremE1P_gauss (F : Fns) (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) (hG : H1G G) (hW : HCW I G)
    (h2 : H2P F I) (h3 : H3P F I) (h4 : H4P I) (θ : ℝ) (hθ1 : 371 / 750 < θ) (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₂ F I)) ∧ h F (371 / 750) (κ₂ F I) ≤ h F θ (κ₂ F I) ∧
      (0.0000335 : ℝ) < h F (371 / 750) (κ₂ F I) :=
  theoremE1P F I (H1P_of_gauss I G hG hW) h2 h3 h4 θ hθ1 hθ2

/-- Non-vacuity of H1G + HCW (G ≡ True, Moment ≡ True). -/
theorem gauss_consistent : H1G (fun _ _ _ _ _ => True) ∧ HCW trivialInputsP (fun _ _ _ _ _ => True) :=
  ⟨fun _ _ _ _ _ _ _ _ _ _ => trivial, fun _ _ _ _ _ _ _ _ _ _ _ _ => trivial⟩

end E1PLean
