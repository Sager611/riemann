import Riemann.Fayed.Defs

/-!
# Inputs.strip and Inputs.reflect hold in Mathlib

Both fields of Inputs that concern only ζ's zeros are theorems of Mathlib
(functional equation riemannZeta_one_sub, riemannZeta_conj, non-vanishing on Re s ≥ 1).
The fields stay in Inputs (fixed interface); these lemmas show they are not real assumptions.
-/

open Complex
open scoped Real ComplexConjugate

namespace Fayed

/-- §2: ζ(s) = 0 with 0 < Re s < 1 ⇒ ζ(1 - s̄) = 0. -/
theorem reflect_holds :
    ∀ s : ℂ, riemannZeta s = 0 → 0 < s.re → s.re < 1 → riemannZeta (1 - conj s) = 0 := by
  intro s hs h0 h1
  have hn : ∀ n : ℕ, conj s ≠ -n := by
    intro n h
    have := congrArg Complex.re h
    simp at this
    linarith [n.cast_nonneg (α := ℝ)]
  have h1' : conj s ≠ 1 := by
    intro h
    have := congrArg Complex.re h
    simp at this
    linarith
  rw [riemannZeta_one_sub hn h1', riemannZeta_conj, hs, map_zero, mul_zero]

/-- §2: nontrivial zeros lie in 0 < Re s < 1. -/
theorem strip_holds : ∀ s : ℂ, riemannZeta s = 0 → (¬∃ n : ℕ, s = -2 * (n + 1)) → s ≠ 1 →
    0 < s.re ∧ s.re < 1 := by
  intro s hs hnt _
  have hlt : s.re < 1 := by
    by_contra h
    exact riemannZeta_ne_zero_of_one_le_re (not_lt.1 h) hs
  refine ⟨?_, hlt⟩
  by_contra hle
  rw [not_lt] at hle
  have hs0 : s ≠ 0 := by
    rintro rfl
    rw [riemannZeta_zero] at hs
    norm_num at hs
  have hn : ∀ n : ℕ, 1 - s ≠ -n := by
    intro n h
    have := congrArg Complex.re h
    simp at this
    linarith [n.cast_nonneg (α := ℝ)]
  have h1 : 1 - s ≠ 1 := fun h => hs0 (by linear_combination -h)
  have hfe := riemannZeta_one_sub hn h1
  rw [sub_sub_cancel, hs] at hfe
  have hre : 1 ≤ (1 - s).re := by simp; linarith
  have hζ := riemannZeta_ne_zero_of_one_le_re hre
  have hΓ : Gamma (1 - s) ≠ 0 := Complex.Gamma_ne_zero_of_re_pos (by linarith)
  have hπ : (π : ℂ) ≠ 0 := ofReal_ne_zero.2 Real.pi_ne_zero
  have h2π : (2 * (π : ℂ)) ^ (-(1 - s)) ≠ 0 := by
    rw [Ne, Complex.cpow_eq_zero_iff]
    exact fun h => mul_ne_zero two_ne_zero hπ h.1
  have hcos : Complex.cos (π * (1 - s) / 2) = 0 := by
    have := hfe.symm
    rcases mul_eq_zero.1 this with h | h
    · rcases mul_eq_zero.1 h with h | h
      · rcases mul_eq_zero.1 h with h | h
        · rcases mul_eq_zero.1 h with h | h
          · norm_num at h
          · exact absurd h h2π
        · exact absurd h hΓ
      · exact h
    · exact absurd h hζ
  obtain ⟨k, hk⟩ := Complex.cos_eq_zero_iff.1 hcos
  have hsk : s = -2 * (k : ℂ) := by
    have h' : (π : ℂ) * (s + 2 * k) = 0 := by linear_combination (-2) * hk
    rcases mul_eq_zero.1 h' with h | h
    · exact absurd h hπ
    · linear_combination h
  have hk0 : (0 : ℤ) ≤ k := by
    have := congrArg Complex.re hsk
    simp at this
    have : (0 : ℝ) ≤ k := by linarith
    exact_mod_cast this
  have hk1 : k ≠ 0 := by
    rintro rfl
    apply hs0
    rw [hsk]
    simp
  obtain ⟨n, rfl⟩ : ∃ n : ℕ, k = n + 1 := ⟨(k - 1).toNat, by omega⟩
  exact hnt ⟨n, by rw [hsk]; push_cast; ring⟩

end Fayed
