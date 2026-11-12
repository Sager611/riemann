import Mathlib.GroupTheory.Perm.Sign

/-!
# Prefix reversals on Fin n: definitions

Positions are 0, ..., n-1 (the NOTE's position x is x.val + 1 here).
rev n i reverses the first i positions: x maps to i - 1 - x for x < i, the identity on x >= i
(NOTE: r_i(x) = i + 1 - x for x <= i, 1-based). For i > n we let it be the identity
(never used: always i <= n). Products act right to left: (s * t) x = s (t x).
-/

namespace PrefixRev

open Equiv

/-- The underlying function of the prefix reversal r_i on Fin n. -/
def revFun (n i : ℕ) (x : Fin n) : Fin n :=
  if h : x.val < i ∧ i ≤ n then ⟨i - 1 - x.val, by omega⟩ else x

theorem revFun_val (n i : ℕ) (x : Fin n) :
    (revFun n i x).val = if x.val < i ∧ i ≤ n then i - 1 - x.val else x.val := by
  unfold revFun; split_ifs <;> rfl

theorem revFun_revFun (n i : ℕ) (x : Fin n) : revFun n i (revFun n i x) = x := by
  apply Fin.ext
  rw [revFun_val, revFun_val]
  by_cases h : x.val < i ∧ i ≤ n
  · rw [ite_eq_left h, ite_eq_left (by omega)]; omega
  · rw [ite_eq_right h, ite_eq_right h]

/-- The prefix reversal r_i as a permutation of Fin n. -/
def rev (n i : ℕ) : Perm (Fin n) where
  toFun := revFun n i
  invFun := revFun n i
  left_inv := revFun_revFun n i
  right_inv := revFun_revFun n i

theorem rev_val (n i : ℕ) (x : Fin n) :
    (rev n i x).val = if x.val < i ∧ i ≤ n then i - 1 - x.val else x.val :=
  revFun_val n i x

theorem rev_rev (n i : ℕ) (x : Fin n) : rev n i (rev n i x) = x := revFun_revFun n i x

theorem rev_mul_self (n i : ℕ) : rev n i * rev n i = 1 := by
  ext x; simp [Perm.mul_apply, rev_rev]

theorem rev_inv (n i : ℕ) : (rev n i)⁻¹ = rev n i := by
  rw [inv_eq_iff_mul_eq_one]; exact rev_mul_self n i

/-- G(n; m, k) = ⟨r_n, r_m, r_k⟩ ≤ S_n. -/
def G (n m k : ℕ) : Subgroup (Perm (Fin n)) :=
  Subgroup.closure {rev n n, rev n m, rev n k}

theorem rn_mem (n m k : ℕ) : rev n n ∈ G n m k := Subgroup.subset_closure (by simp)
theorem rm_mem (n m k : ℕ) : rev n m ∈ G n m k := Subgroup.subset_closure (by simp)
theorem rk_mem (n m k : ℕ) : rev n k ∈ G n m k := Subgroup.subset_closure (by simp)

/-- G ≤ H as soon as H contains the three generators. -/
theorem G_le {n m k : ℕ} {H : Subgroup (Perm (Fin n))} (hn : rev n n ∈ H) (hm : rev n m ∈ H)
    (hk : rev n k ∈ H) : G n m k ≤ H := by
  rw [G, Subgroup.closure_le]
  intro g hg
  simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hg
  rcases hg with rfl | rfl | rfl <;> assumption

/-- Permutations fixing every fibre of b (setwise). -/
def fibreFix {n : ℕ} {β : Type} (b : Fin n → β) : Subgroup (Perm (Fin n)) where
  carrier := {g | ∀ x, b (g x) = b x}
  one_mem' := fun _ => rfl
  mul_mem' := by
    intro g h hg hh x
    simp only [Set.mem_ofPred_eq] at hg hh ⊢
    rw [Perm.mul_apply, hg, hh]
  inv_mem' := by
    intro g hg x
    simp only [Set.mem_ofPred_eq] at hg ⊢
    have := hg (g⁻¹ x)
    rw [Perm.inv_def, Equiv.apply_symm_apply] at this
    exact this.symm

/-- Permutations preserving the partition into fibres of b (a block system). -/
def fibreStab {n : ℕ} {β : Type} (b : Fin n → β) : Subgroup (Perm (Fin n)) where
  carrier := {g | ∀ x y, b x = b y ↔ b (g x) = b (g y)}
  one_mem' := fun _ _ => Iff.rfl
  mul_mem' := by
    intro g h hg hh x y
    simp only [Set.mem_ofPred_eq] at hg hh ⊢
    rw [Perm.mul_apply, Perm.mul_apply, hh x y, hg (h x) (h y)]
  inv_mem' := by
    intro g hg x y
    simp only [Set.mem_ofPred_eq] at hg ⊢
    have := hg (g⁻¹ x) (g⁻¹ y)
    rw [Perm.inv_def, Equiv.apply_symm_apply, Equiv.apply_symm_apply] at this
    exact this.symm

end PrefixRev
