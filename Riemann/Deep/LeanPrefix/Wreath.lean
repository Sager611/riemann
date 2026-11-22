import Mathlib.GroupTheory.SpecificGroups.Dihedral
import Mathlib.GroupTheory.SemidirectProduct
import Mathlib.Tactic.Ring
import Mathlib.Tactic.LinearCombination

/-!
# The imprimitive wreath product D_{2α} ≀ D_{2q}

D_{2q} is Mathlib's DihedralGroup q (order 2q), acting on the vertices ZMod q of the q-gon by
dihedralPerm q: r i : x ↦ x + i, sr i : x ↦ -x - i. The wreath product D_{2α} ≀ D_{2q} (with
respect to this action) is the semidirect product (ZMod q → DihedralGroup α) ⋊ DihedralGroup q,
d acting on the base by f ↦ f ∘ (d⁻¹ · ); it acts imprimitively on ZMod q × ZMod α by
(f, d) · (i, u) = (d i, f (d i) u). Both actions are faithful for α, q ≥ 3.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv DihedralGroup

/-- The underlying permutations of the q-gon action. -/
def dihedralPermFun (q : ℕ) : DihedralGroup q → Perm (ZMod q)
  | r i => Equiv.addRight i
  | sr i => (Equiv.neg (ZMod q)).trans (Equiv.subRight i)

theorem dihedralPermFun_r {q : ℕ} (i x : ZMod q) : dihedralPermFun q (r i) x = x + i := rfl

theorem dihedralPermFun_sr {q : ℕ} (i x : ZMod q) : dihedralPermFun q (sr i) x = -x - i := rfl

/-- D_{2q} = DihedralGroup q acting on the vertices ZMod q of the q-gon:
r i : x ↦ x + i (rotations), sr i : x ↦ -x - i (reflections). -/
def dihedralPerm (q : ℕ) : DihedralGroup q →* Perm (ZMod q) where
  toFun := dihedralPermFun q
  map_one' := by
    rw [DihedralGroup.one_def]; ext x; rw [dihedralPermFun_r, Perm.one_apply, add_zero]
  map_mul' := by
    rintro (i | i) (j | j) <;> ext x <;>
      simp only [r_mul_r, r_mul_sr, sr_mul_r, sr_mul_sr, Perm.mul_apply, dihedralPermFun_r,
        dihedralPermFun_sr] <;> ring

theorem dihedralPerm_r {q : ℕ} (i x : ZMod q) : dihedralPerm q (r i) x = x + i := rfl

theorem dihedralPerm_sr {q : ℕ} (i x : ZMod q) : dihedralPerm q (sr i) x = -x - i := rfl

/-- The q-gon action is faithful for q ≥ 3. -/
theorem dihedralPerm_injective {q : ℕ} (hq : 3 ≤ q) : Function.Injective (dihedralPerm q) := by
  rw [injective_iff_map_eq_one]
  rintro (i | i) h
  · have h0 := congrArg (fun p => p 0) h
    simp only [dihedralPerm_r, Perm.one_apply, zero_add] at h0
    rw [h0, DihedralGroup.one_def]
  · have h0 := congrArg (fun p => p 0) h
    simp only [dihedralPerm_sr, Perm.one_apply, neg_zero, zero_sub, neg_eq_zero] at h0
    have h1 := congrArg (fun p => p 1) h
    simp only [dihedralPerm_sr, Perm.one_apply, h0, sub_zero] at h1
    have h2 : ((2 : ℕ) : ZMod q) = 0 := by
      rw [Nat.cast_ofNat]; linear_combination -h1
    rw [ZMod.natCast_eq_zero_iff] at h2
    have := Nat.le_of_dvd (by norm_num) h2
    omega

/-- DihedralGroup q acting on the base ZMod q → DihedralGroup α by f ↦ f ∘ (d⁻¹ · ). -/
def wreathAct (α q : ℕ) : DihedralGroup q →* MulAut (ZMod q → DihedralGroup α) where
  toFun d :=
    { toFun := fun f i => f ((dihedralPerm q d)⁻¹ i)
      invFun := fun f i => f (dihedralPerm q d i)
      left_inv := fun f => by funext i; simp
      right_inv := fun f => by funext i; simp
      map_mul' := fun f g => rfl }
  map_one' := by ext f i; simp
  map_mul' := by intro a b; ext f i; simp [mul_inv_rev, Perm.mul_apply]

theorem wreathAct_apply {α q : ℕ} (d : DihedralGroup q) (f : ZMod q → DihedralGroup α) (i : ZMod q) :
    wreathAct α q d f i = f ((dihedralPerm q d)⁻¹ i) := rfl

/-- The wreath product D_{2α} ≀ D_{2q} for the q-gon action. -/
abbrev Wreath (α q : ℕ) := (ZMod q → DihedralGroup α) ⋊[wreathAct α q] DihedralGroup q

/-- The imprimitive action of D_{2α} ≀ D_{2q} on ZMod q × ZMod α:
(f, d) · (i, u) = (d i, f (d i) u). -/
def wreathPerm (α q : ℕ) : Wreath α q →* Perm (ZMod q × ZMod α) where
  toFun w := Equiv.prodShear (dihedralPerm q w.right) (fun i => dihedralPerm α (w.left (dihedralPerm q w.right i)))
  map_one' := by
    ext ⟨i, u⟩ <;> simp [Equiv.prodShear]
  map_mul' := by
    intro a b
    ext ⟨i, u⟩ <;> simp [Equiv.prodShear, SemidirectProduct.mul_left, SemidirectProduct.mul_right,
      wreathAct_apply, Perm.mul_apply]

theorem wreathPerm_apply {α q : ℕ} (w : Wreath α q) (i : ZMod q) (u : ZMod α) :
    wreathPerm α q w (i, u) =
      (dihedralPerm q w.right i, dihedralPerm α (w.left (dihedralPerm q w.right i)) u) := rfl

/-- The imprimitive action is faithful for α, q ≥ 3. -/
theorem wreathPerm_injective {α q : ℕ} (hα : 3 ≤ α) (hq : 3 ≤ q) :
    Function.Injective (wreathPerm α q) := by
  rw [injective_iff_map_eq_one]
  intro w hw
  have h1 : ∀ i u, wreathPerm α q w (i, u) = (i, u) := fun i u => by rw [hw]; rfl
  have hd : dihedralPerm q w.right = 1 := by
    ext i; exact congrArg Prod.fst (h1 i 0)
  have hd1 : w.right = 1 := dihedralPerm_injective hq (by rw [hd, map_one])
  have hf : ∀ i, w.left i = 1 := by
    intro i
    apply dihedralPerm_injective hα
    rw [map_one]
    ext u
    have := congrArg Prod.snd (h1 i u)
    simp only [wreathPerm_apply, hd, Perm.one_apply] at this
    exact this
  ext i
  · rw [SemidirectProduct.one_left]; exact hf i
  · rw [hd1, SemidirectProduct.one_right]

/-- |D_{2α} ≀ D_{2q}| = (2α)^q · 2q. -/
theorem card_wreath (α q : ℕ) [NeZero q] : Nat.card (Wreath α q) = (2 * α) ^ q * (2 * q) := by
  rw [SemidirectProduct.card, Nat.card_fun, DihedralGroup.nat_card, DihedralGroup.nat_card,
    Nat.card_zmod]

end Riemann.Deep.LeanPrefix

