import Riemann.Deep.EvoAlg.Lean.Defs

/-!
# Lemma 2 (NOTE §4): the norm identity inside a twin class

Square form (the NOTE's statement): for `P ∈ M_k(𝕜)` (entries `P w a`, row `w`, column `a`),
`k(γ² − α) = 𝟙ᵀ(PᵀP − αI)𝟙 − (P𝟙 − γ𝟙)ᵀ(P𝟙 + γ𝟙)`, so `PᵀP = αI` and `P𝟙 = γ𝟙` give
`k(γ² − α) = 0`. Block form (as used in Theorem I: columns `S`, rows `S' = τS`, entry
`P a w = (φ e_a)_w`): `|S| α = |S'| γ²`.
-/

open Finset

namespace EvoAlg

variable {𝕜 : Type*} [Field 𝕜]

/-- The polynomial identity of Lemma 2. -/
theorem lemma2_identity {ι : Type*} [Fintype ι] [DecidableEq ι] (P : ι → ι → 𝕜) (α γ : 𝕜) :
    (Fintype.card ι : 𝕜) * (γ ^ 2 - α)
      = ∑ a, ∑ b, ((∑ w, P w a * P w b) - if a = b then α else 0)
        - ∑ w, ((∑ a, P w a) - γ) * ((∑ a, P w a) + γ) := by
  have hsq : ∑ a : ι, ∑ b : ι, ∑ w, P w a * P w b = ∑ w, (∑ a, P w a) ^ 2 := by
    calc ∑ a : ι, ∑ b : ι, ∑ w, P w a * P w b = ∑ a : ι, ∑ w, ∑ b : ι, P w a * P w b :=
          Finset.sum_congr rfl (fun a _ => Finset.sum_comm)
      _ = ∑ w, ∑ a : ι, ∑ b : ι, P w a * P w b := Finset.sum_comm
      _ = ∑ w, (∑ a, P w a) ^ 2 := Finset.sum_congr rfl (fun w _ => by rw [sq, Finset.sum_mul_sum])
  have h1 : ∑ a : ι, ∑ b : ι, ((∑ w, P w a * P w b) - if a = b then α else 0)
      = ∑ w, (∑ a, P w a) ^ 2 - (Fintype.card ι : 𝕜) * α := by
    simp only [Finset.sum_sub_distrib, Finset.sum_ite_eq, Finset.mem_univ, ↓reduceIte,
      Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    rw [hsq]
  have h2 : ∑ w : ι, ((∑ a, P w a) - γ) * ((∑ a, P w a) + γ)
      = ∑ w, (∑ a, P w a) ^ 2 - (Fintype.card ι : 𝕜) * γ ^ 2 := by
    have : ∀ w : ι, ((∑ a, P w a) - γ) * ((∑ a, P w a) + γ) = (∑ a, P w a) ^ 2 - γ ^ 2 :=
      fun w => by ring
    simp only [this, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  rw [h1, h2]
  ring

/-- **Lemma 2** (NOTE §4): `PᵀP = αI` and `P𝟙 = γ𝟙` imply `k(γ² − α) = 0`. -/
theorem lemma2 {ι : Type*} [Fintype ι] [DecidableEq ι] (P : ι → ι → 𝕜) (α γ : 𝕜)
    (ha : ∀ a b, ∑ w, P w a * P w b = if a = b then α else 0) (hb : ∀ w, ∑ a, P w a = γ) :
    (Fintype.card ι : 𝕜) * (γ ^ 2 - α) = 0 := by
  rw [lemma2_identity P α γ]
  simp [ha, hb]

/-- Lemma 2, block form: columns `S`, rows `S'`, `Σ_{w∈S'} P a w P b w = δ_{ab} α` and
`Σ_{a∈S} P a w = γ` (`w ∈ S'`) give `|S| α = |S'| γ²`. -/
theorem lemma2_blocks {V : Type*} [DecidableEq V] (S S' : Finset V) (P : V → V → 𝕜) (α γ : 𝕜)
    (ha : ∀ a ∈ S, ∀ b ∈ S, ∑ w ∈ S', P a w * P b w = if a = b then α else 0)
    (hb : ∀ w ∈ S', ∑ a ∈ S, P a w = γ) : (#S : 𝕜) * α = #S' * γ ^ 2 := by
  have h1 : ∑ a ∈ S, ∑ b ∈ S, ∑ w ∈ S', P a w * P b w = #S * α := by
    calc ∑ a ∈ S, ∑ b ∈ S, ∑ w ∈ S', P a w * P b w
        = ∑ a ∈ S, ∑ b ∈ S, (if a = b then α else 0) :=
          Finset.sum_congr rfl (fun a ha' => Finset.sum_congr rfl (fun b hb' => ha a ha' b hb'))
      _ = ∑ a ∈ S, α := Finset.sum_congr rfl (fun a ha' => by rw [Finset.sum_ite_eq, ite_eq_left ha'])
      _ = #S * α := by simp
  have h2 : ∑ a ∈ S, ∑ b ∈ S, ∑ w ∈ S', P a w * P b w = #S' * γ ^ 2 := by
    calc ∑ a ∈ S, ∑ b ∈ S, ∑ w ∈ S', P a w * P b w
        = ∑ a ∈ S, ∑ w ∈ S', ∑ b ∈ S, P a w * P b w :=
          Finset.sum_congr rfl (fun a _ => Finset.sum_comm)
      _ = ∑ w ∈ S', ∑ a ∈ S, ∑ b ∈ S, P a w * P b w := Finset.sum_comm
      _ = ∑ w ∈ S', (∑ a ∈ S, P a w) * (∑ b ∈ S, P b w) :=
          Finset.sum_congr rfl (fun w _ => (Finset.sum_mul_sum _ _ _ _).symm)
      _ = ∑ w ∈ S', γ ^ 2 := Finset.sum_congr rfl (fun w hw => by rw [hb w hw, sq])
      _ = #S' * γ ^ 2 := by simp
  rw [← h1, ← h2]

end EvoAlg

