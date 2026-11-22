import Riemann.Deep.LeanPrefix.WreathC
import Riemann.Deep.LeanPrefix.OrbitsB
import Mathlib.Data.Set.Card

/-!
# Theorem B.2: n even, gcd(m+1, k+1) = 1: G is transitive and imprimitive with two blocks

The parity classes of the cyclic coordinate phi are blocks (r_n flips the parity of phi since
n + 1 is odd; r_k and P keep it); r_n maps one class onto the other, so each has n/2 points.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv

theorem G_le_parity {n m k θ : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (hθ : (k + 1) * θ % (n + 1) = 1) (hev : n % 2 = 0) :
    G n m k ≤ fibreStab (fun x : Fin n => phi n θ x % 2) := by
  have hn : rev n n ∈ fibreStab (fun x : Fin n => phi n θ x % 2) := by
    apply mem_fibreStab_of_invol _ _ (rev_rev n n) (fun c => (c + 1) % 2)
    intro x
    show phi n θ (rev n n x) % 2 = (phi n θ x % 2 + 1) % 2
    have h1 := phi_le (θ := θ) x
    rw [phi_rn hline hkm hθ x]; omega
  have hk : rev n k ∈ fibreStab (fun x : Fin n => phi n θ x % 2) := by
    apply mem_fibreStab_of_invol _ _ (rev_rev n k) id
    intro x
    show phi n θ (rev n k x) % 2 = phi n θ x % 2
    have h1 := phi_le (θ := θ) x
    rw [phi_rk hline hkm hθ x]; split_ifs <;> omega
  have hP : P n m k ∈ fibreStab (fun x : Fin n => phi n θ x % 2) := by
    intro x y
    show phi n θ x % 2 = phi n θ y % 2 ↔ phi n θ (P n m k x) % 2 = phi n θ (P n m k y) % 2
    have h1 := phi_le (θ := θ) x
    have h2 := phi_le (θ := θ) y
    rw [phi_P hline hkm hθ x, phi_P hline hkm hθ y]
    split_ifs <;> omega
  apply G_le hn _ hk
  rw [rm_eq]
  exact (fibreStab _).mul_mem ((fibreStab _).mul_mem hn hk) hP

/-- Theorem B.2: n even, gcd = 1: G is transitive and imprimitive, with a block of n/2 points. -/
theorem line_B2 {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m) (hline : n = m + k + 1)
    (hev : n % 2 = 0) (hg : Nat.Coprime (k + 1) (m + 1)) :
    MulAction.IsPretransitive (G n m k) (Fin n) ∧ ¬ MulAction.IsPreprimitive (G n m k) (Fin n) ∧
      ∃ B : Set (Fin n), MulAction.IsBlock (G n m k) B ∧ 2 * B.ncard = n := by
  obtain ⟨θ, hθ⟩ := exists_theta hline hg
  have htr : MulAction.IsPretransitive (G n m k) (Fin n) := (isPretransitive_iff_gcd hline hkm).mpr hg
  set B : Set (Fin n) := {x | phi n θ x % 2 = 0} with hBdef
  have hB : MulAction.IsBlock (G n m k) B := by
    rw [isBlock_iff_isBlockOf]
    exact isBlockOf_fibre _ _ (G_le_parity hline hkm hθ hev) 0
  have himg : ⇑(rev n n) '' B = Bᶜ := by
    ext y
    simp only [Set.mem_image, Set.mem_compl_iff, hBdef, Set.mem_setOf_eq]
    constructor
    · rintro ⟨x, hx, rfl⟩
      have h1 := phi_le (θ := θ) x
      rw [phi_rn hline hkm hθ x]; omega
    · intro hy
      refine ⟨rev n n y, ?_, rev_rev n n y⟩
      have h1 := phi_le (θ := θ) y
      rw [phi_rn hline hkm hθ y]; omega
  have hcard : 2 * B.ncard = n := by
    have h1 := Set.ncard_add_ncard_compl B
    rw [← himg, Set.ncard_image_of_injective _ (rev n n).injective, Nat.card_eq_fintype_card,
      Fintype.card_fin] at h1
    omega
  refine ⟨htr, fun hp => ?_, B, hB, hcard⟩
  rcases hp.isTrivialBlock_of_isBlock hB with h | h
  · have : B.ncard ≤ 1 := (Set.ncard_le_one_iff_subsingleton).mpr h
    omega
  · have : B.ncard = n := by rw [h, Set.ncard_univ, Nat.card_eq_fintype_card, Fintype.card_fin]
    omega

end Riemann.Deep.LeanPrefix

