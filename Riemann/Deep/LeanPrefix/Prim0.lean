import Riemann.Deep.LeanPrefix.BlocksB

/-!
# Blocks, primitivity (elementary form) and small facts about r_k and P

IsBlockOf / IsPrimitive are the elementary forms; Bridge.lean proves they agree with Mathlib's
MulAction.IsBlock / MulAction.IsPreprimitive.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv

/-- B is a block of H: every h ∈ H maps B onto B or onto a set disjoint from B. -/
def IsBlockOf {n : ℕ} (H : Subgroup (Perm (Fin n))) (B : Set (Fin n)) : Prop :=
  ∀ h ∈ H, (⇑h '' B = B) ∨ Disjoint (⇑h '' B) B

/-- H is primitive: transitive, and every block is trivial. -/
def IsPrimitive {n : ℕ} (H : Subgroup (Perm (Fin n))) : Prop :=
  (∀ x y : Fin n, ∃ h ∈ H, h x = y) ∧
    ∀ B : Set (Fin n), IsBlockOf H B → B.Subsingleton ∨ B = Set.univ

/-- r_k fixes the n - k points x ≥ k. -/
theorem fix_rk_card {n k : ℕ} (hk : k ≤ n) :
    n - k ≤ (Finset.univ.filter (fun x : Fin n => rev n k x = x)).card := by
  have h := Finset.card_le_card_of_injOn (s := (Finset.univ : Finset (Fin (n - k))))
    (t := Finset.univ.filter (fun x : Fin n => rev n k x = x))
    (fun i => (⟨k + i.val, by have := i.isLt; omega⟩ : Fin n))
    (by
      intro i _
      have hi := i.isLt
      rw [Finset.mem_coe, Finset.mem_filter]
      refine ⟨Finset.mem_univ _, ?_⟩
      apply Fin.ext
      rw [rev_val]
      dsimp only
      rw [ite_eq_right (by omega)])
    (by
      intro i _ j _ hij
      have := congrArg Fin.val hij
      exact Fin.ext (by simp only at this; omega))
  rwa [Finset.card_univ, Fintype.card_fin] at h

theorem rk_ne_one {n k : ℕ} (hk : 2 ≤ k) (hkn : k ≤ n) : rev n k ≠ 1 := by
  intro h
  have h0 : (rev n k ⟨0, by omega⟩).val = k - 1 := by
    rw [rev_val, ite_eq_left (show (0 : ℕ) < k ∧ k ≤ n from ⟨by omega, hkn⟩)]
    dsimp only
    omega
  rw [h, Perm.one_apply] at h0
  dsimp only at h0
  omega

theorem P_mem_G (n m k : ℕ) : P n m k ∈ G n m k :=
  (G n m k).mul_mem ((G n m k).mul_mem (rk_mem n m k) (rn_mem n m k)) (rm_mem n m k)

theorem even_of_coprime {n m k : ℕ} (hline : n = m + k + 1) (hodd : n % 2 = 1)
    (hg : Nat.Coprime (m + 1) (k + 1)) : m % 2 = 0 ∧ k % 2 = 0 := by
  by_contra h
  have h2 : 2 ∣ Nat.gcd (m + 1) (k + 1) := Nat.dvd_gcd (by omega) (by omega)
  rw [hg] at h2
  exact absurd (Nat.le_of_dvd (by omega) h2) (by omega)


end Riemann.Deep.LeanPrefix
