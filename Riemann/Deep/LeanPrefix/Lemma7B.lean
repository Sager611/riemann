import Riemann.Deep.LeanPrefix.Lemma7A
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-!
# Lemma 7 (involution form), part B: counting helpers, quadratic Hilbert 90, |F| = |F₀|²
-/

namespace Riemann.Deep.LeanPrefix

open Equiv
open scoped LinearAlgebra.Projectivization

section

variable {F : Type} [Field F] [Fintype F] [DecidableEq F] {d : ℕ}

theorem card_vec : Fintype.card (Fin d → F) = Fintype.card F ^ d := by simp

/-- If the sum map A × B → V is injective, |A| |B| ≤ |V|. -/
theorem card_mul_le_of_add_inj (A B : Finset (Fin d → F)) (u : F)
    (h : ∀ a ∈ A, ∀ b ∈ B, ∀ a' ∈ A, ∀ b' ∈ B, a + u • b = a' + u • b' → a = a' ∧ b = b') :
    A.card * B.card ≤ Fintype.card (Fin d → F) := by
  rw [← Finset.card_product, ← Finset.card_univ]
  apply Finset.card_le_card_of_injOn (fun p => p.1 + u • p.2) (fun _ _ => Finset.mem_coe.mpr (Finset.mem_univ _))
  intro p hp p' hp' hpp
  simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe] at hp hp'
  obtain ⟨h1, h2⟩ := h p.1 hp.1 p.2 hp.2 p'.1 hp'.1 p'.2 hp'.2 hpp
  exact Prod.ext h1 h2

/-- A proper subspace has at most |V|/q elements. -/
theorem card_mul_q_le (S : Finset (Fin d → F)) (hsub : ∀ a ∈ S, ∀ b ∈ S, a - b ∈ S)
    (hsmul : ∀ (c : F), ∀ a ∈ S, c • a ∈ S) {w : Fin d → F} (hw : w ∉ S) :
    S.card * Fintype.card F ≤ Fintype.card (Fin d → F) := by
  rw [← Finset.card_univ (α := F), ← Finset.card_product, ← Finset.card_univ (α := Fin d → F)]
  apply Finset.card_le_card_of_injOn (fun p => p.1 + p.2 • w) (fun _ _ => Finset.mem_coe.mpr (Finset.mem_univ _))
  intro p hp p' hp' hpp
  simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe] at hp hp'
  simp only at hpp
  have hc : p.2 = p'.2 := by
    by_contra hne
    apply hw
    have h1 : (p.2 - p'.2) • w = p'.1 - p.1 := by
      rw [sub_smul]; rw [← sub_eq_zero] at hpp ⊢; rw [← hpp]; abel
    have h2 : w = (p.2 - p'.2)⁻¹ • (p'.1 - p.1) := by
      rw [← h1, smul_smul, inv_mul_cancel₀ (sub_ne_zero.mpr hne), one_smul]
    rw [h2]; exact hsmul _ _ (hsub _ hp'.1 _ hp.1)
  have ha : p.1 = p'.1 := by rw [hc] at hpp; exact add_right_cancel hpp
  exact Prod.ext ha hc

/-- A subspace containing v ≠ 0 has at least q elements. -/
theorem q_le_card (S : Finset (Fin d → F)) (hsmul : ∀ (c : F), ∀ a ∈ S, c • a ∈ S)
    {v : Fin d → F} (hv : v ∈ S) (hv0 : v ≠ 0) : Fintype.card F ≤ S.card := by
  rw [← Finset.card_univ]
  apply Finset.card_le_card_of_injOn (fun c => c • v) (fun c _ => Finset.mem_coe.mpr (hsmul c v hv))
  intro c _ c' _ h
  simp only at h
  by_contra hne
  have : (c - c') • v = 0 := by rw [sub_smul, h, sub_self]
  rcases smul_eq_zero.mp this with h1 | h1
  · exact hne (sub_eq_zero.mp h1)
  · exact hv0 h1

/-- Quadratic Hilbert 90: if σ² = id, σ ≠ id and ν σ(ν) = 1, then ν = σ(c)/c for some c ≠ 0. -/
theorem hilbert90_quad {σ : F ≃+* F} (hσ2 : ∀ c, σ (σ c) = c) {ω : F} (hω : σ ω ≠ ω) {ν : F}
    (hν : ν * σ ν = 1) : ∃ c : F, c ≠ 0 ∧ σ c = ν * c := by
  have hν0 : ν ≠ 0 := by rintro rfl; rw [zero_mul] at hν; exact zero_ne_one hν
  have hσν : σ ν = ν⁻¹ := eq_inv_of_mul_eq_one_right hν
  by_cases h1 : 1 + ν⁻¹ = 0
  · have hm : ν = -1 := by
      have : ν⁻¹ = -1 := by linear_combination h1
      rw [← inv_inv ν, this, inv_neg, inv_one]
    refine ⟨ω - σ ω, sub_ne_zero.mpr (Ne.symm hω), ?_⟩
    rw [map_sub, hσ2, hm]; ring
  · refine ⟨1 + ν⁻¹, h1, ?_⟩
    rw [map_add, map_one, map_inv₀, hσν, inv_inv, mul_add, mul_one, mul_inv_cancel₀ hν0, add_comm]

/-- |F| = |F₀|² for the fixed field F₀ of an involution σ ≠ id. -/
theorem card_eq_sq_fixed {σ : F ≃+* F} (hσ2 : ∀ c, σ (σ c) = c) {ω : F} (hω : σ ω ≠ ω) :
    Fintype.card F = (Finset.univ.filter (fun c : F => σ c = c)).card ^ 2 := by
  set F0 := Finset.univ.filter (fun c : F => σ c = c) with hF0
  have hmem : ∀ c, c ∈ F0 ↔ σ c = c := fun c => by simp [hF0]
  have hωω : ω - σ ω ≠ 0 := sub_ne_zero.mpr (Ne.symm hω)
  -- the map F0 × F0 → F, (a, b) ↦ a + b ω
  have hinj : (F0 ×ˢ F0).card ≤ Fintype.card F := by
    rw [← Finset.card_univ]
    apply Finset.card_le_card_of_injOn (fun p => p.1 + p.2 * ω)
      (fun _ _ => Finset.mem_coe.mpr (Finset.mem_univ _))
    intro p hp p' hp' hpp
    simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, hmem] at hp hp'
    simp only at hpp
    have hb : p.2 = p'.2 := by
      by_contra hne
      apply hω
      have h1 : ω = (p'.1 - p.1) / (p.2 - p'.2) := by
        rw [eq_div_iff (sub_ne_zero.mpr hne)]; linear_combination hpp
      rw [h1, map_div₀, map_sub, map_sub, hp.1, hp.2, hp'.1, hp'.2]
    have ha : p.1 = p'.1 := by rw [hb] at hpp; exact add_right_cancel hpp
    exact Prod.ext ha hb
  have hsurj : Fintype.card F ≤ (F0 ×ˢ F0).card := by
    rw [← Finset.card_univ]
    apply Finset.card_le_card_of_injOn
      (fun x => ((x - (x - σ x) / (ω - σ ω) * ω), (x - σ x) / (ω - σ ω)))
    · intro x _
      simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, hmem]
      have hb : σ ((x - σ x) / (ω - σ ω)) = (x - σ x) / (ω - σ ω) := by
        rw [map_div₀, map_sub, map_sub, hσ2, hσ2]
        rw [div_eq_div_iff (sub_ne_zero.mpr hω) hωω]; ring
      refine ⟨?_, hb⟩
      rw [map_sub, map_mul, hb]
      field_simp
      ring
    · intro x _ x' _ h
      simp only [Prod.mk.injEq] at h
      have := h.1
      rw [h.2] at this
      exact sub_left_injective this
  rw [Finset.card_product] at hinj hsurj
  rw [sq]; omega

end

end Riemann.Deep.LeanPrefix

