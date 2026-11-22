import Riemann.Deep.LeanPrefix.Wreath
import Riemann.Deep.LeanPrefix.BlocksA

/-!
# Lemma 6 (B.3), part 1: block coordinates and the generators as wreath-product elements

α = k + 1, n = q α, m = (q - 1) α. Block coordinates x ↦ (x / α, x mod α) ∈ ZMod q × ZMod α
(block index, offset; 0-based). In these coordinates (NOTE Lemma 6):
r_n = (const sr 1, sr 1): (i, u) ↦ (-i - 1, -u - 1);
r_m = (j ↦ if j = -1 then 1 else sr 1, sr 2): (i, u) ↦ (-i - 2, -u - 1), block q - 1 fixed;
r_k = (mulSingle 0 (sr 2), 1): (0, u) ↦ (0, -u - 2), other blocks fixed.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv DihedralGroup

/-- Block coordinates (block index, offset). -/
def blockFun (q α : ℕ) {n : ℕ} (x : Fin n) : ZMod q × ZMod α :=
  (((x.val / α : ℕ) : ZMod q), ((x.val % α : ℕ) : ZMod α))

theorem blockFun_injective {n q α : ℕ} (hn : n = q * α) (hα : 0 < α) :
    Function.Injective (blockFun q α (n := n)) := by
  intro x y h
  simp only [blockFun, Prod.mk.injEq] at h
  obtain ⟨h1, h2⟩ := h
  rw [ZMod.natCast_eq_natCast_iff'] at h1 h2
  have hx := x.isLt; have hy := y.isLt
  have hxq : x.val / α < q := (Nat.div_lt_iff_lt_mul hα).mpr (by rw [← hn]; exact hx)
  have hyq : y.val / α < q := (Nat.div_lt_iff_lt_mul hα).mpr (by rw [← hn]; exact hy)
  rw [Nat.mod_eq_of_lt hxq, Nat.mod_eq_of_lt hyq] at h1
  rw [Nat.mod_mod, Nat.mod_mod] at h2
  apply Fin.ext
  rw [← Nat.div_add_mod x.val α, ← Nat.div_add_mod y.val α, h1, h2]

/-- The block-coordinate bijection Fin n ≃ ZMod q × ZMod α (n = q α). -/
noncomputable def blockEquiv {n q α : ℕ} [NeZero q] [NeZero α] (hn : n = q * α) :
    Fin n ≃ ZMod q × ZMod α :=
  Equiv.ofBijective (blockFun q α) ((Fintype.bijective_iff_injective_and_card _).mpr
    ⟨blockFun_injective hn (Nat.pos_of_ne_zero (NeZero.ne α)), by
      rw [Fintype.card_prod, ZMod.card, ZMod.card, Fintype.card_fin, hn]⟩)

theorem blockEquiv_apply {n q α : ℕ} [NeZero q] [NeZero α] (hn : n = q * α) (x : Fin n) :
    blockEquiv hn x = (((x.val / α : ℕ) : ZMod q), ((x.val % α : ℕ) : ZMod α)) := rfl

theorem mod_reflect {a M x : ℕ} (ha : 0 < a) (hx : x < M * a) :
    (M * a - 1 - x) % a = a - 1 - x % a := by
  have h1 := Nat.div_add_mod (M * a - 1 - x) a
  rw [div_reflect ha hx] at h1
  have h2 := Nat.div_add_mod x a
  have h3 : x / a < M := (Nat.div_lt_iff_lt_mul ha).mpr hx
  have h4 := Nat.mod_lt x ha
  have h5 : a * (M - 1 - x / a) + a * (x / a) + a = M * a := by
    rw [← Nat.mul_add, ← Nat.mul_succ, Nat.mul_comm]; congr 1; omega
  omega

theorem cast_sub_sub (M a b : ℕ) (h : a + b ≤ M) :
    ((M - b - a : ℕ) : ZMod M) = -(a : ZMod M) - b := by
  rw [Nat.cast_sub (by omega), Nat.cast_sub (by omega), ZMod.natCast_self]; ring

theorem natCast_ne_zero_of_lt {M a : ℕ} (h1 : 0 < a) (h2 : a < M) : (a : ZMod M) ≠ 0 := by
  rw [Ne, ZMod.natCast_eq_zero_iff]
  intro h; have := Nat.le_of_dvd h1 h; omega

/-- The generators as elements of D_{2α} ≀ D_{2q}. -/
def wN (α q : ℕ) : Wreath α q := ⟨fun _ => sr 1, sr 1⟩

def wM (α q : ℕ) : Wreath α q := ⟨fun j => if j = -1 then 1 else sr 1, sr 2⟩

def wK (α q : ℕ) : Wreath α q := ⟨Pi.mulSingle 0 (sr 2), 1⟩

section gens

variable {n m k q α : ℕ} [NeZero q] [NeZero α] (hn : n = q * α) (hm : m = (q - 1) * α)
  (hk : k + 1 = α) (hq : 3 ≤ q) (hα : 3 ≤ α)
include hn hq hα

theorem block_rn (x : Fin n) :
    blockEquiv hn (rev n n x) = wreathPerm α q (wN α q) (blockEquiv hn x) := by
  have hx := x.isLt
  have hxn : x.val < q * α := by rw [← hn]; exact hx
  have hxq : x.val / α < q := (Nat.div_lt_iff_lt_mul (by omega)).mpr hxn
  have hxα := Nat.mod_lt x.val (show 0 < α by omega)
  have hv : (rev n n x).val = q * α - 1 - x.val := by
    rw [rev_val, ite_eq_left ⟨hx, le_rfl⟩]; omega
  rw [blockEquiv_apply, blockEquiv_apply, wreathPerm_apply, hv, div_reflect (by omega) hxn,
    mod_reflect (by omega) hxn]
  simp only [wN, dihedralPerm_sr]
  rw [cast_sub_sub q _ 1 (by omega), cast_sub_sub α _ 1 (by omega)]
  simp

include hm in
theorem block_rm (x : Fin n) :
    blockEquiv hn (rev n m x) = wreathPerm α q (wM α q) (blockEquiv hn x) := by
  have hx := x.isLt
  have hxn : x.val < q * α := by rw [← hn]; exact hx
  have hxq : x.val / α < q := (Nat.div_lt_iff_lt_mul (by omega)).mpr hxn
  have hxα := Nat.mod_lt x.val (show 0 < α by omega)
  rw [blockEquiv_apply, blockEquiv_apply, wreathPerm_apply]
  simp only [wM, dihedralPerm_sr]
  by_cases hxm : x.val < m
  · have hxm' : x.val < (q - 1) * α := by rw [← hm]; exact hxm
    have hxq1 : x.val / α < q - 1 := (Nat.div_lt_iff_lt_mul (by omega)).mpr hxm'
    have hv : (rev n m x).val = (q - 1) * α - 1 - x.val := by
      have hmn : m ≤ n := by
        have : (q - 1) * α ≤ q * α := Nat.mul_le_mul_right α (by omega)
        omega
      rw [rev_val, ite_eq_left ⟨hxm, hmn⟩]; omega
    rw [hv, div_reflect (by omega) hxm', mod_reflect (by omega) hxm']
    have hne : -((x.val / α : ℕ) : ZMod q) - 2 ≠ -1 := by
      intro h
      have h1 : (((x.val / α + 1 : ℕ)) : ZMod q) = 0 := by push_cast; linear_combination -h
      exact natCast_ne_zero_of_lt (M := q) (a := x.val / α + 1) (Nat.succ_pos _) (by omega) h1
    rw [if_neg hne, dihedralPerm_sr]
    rw [show q - 1 - 1 - x.val / α = q - 2 - x.val / α by omega, cast_sub_sub q _ 2 (by omega),
      cast_sub_sub α _ 1 (by omega)]
    push_cast; rfl
  · have hv : rev n m x = x := Fin.ext (by rw [rev_val, ite_eq_right (by omega)])
    rw [hv]
    have hb : x.val / α = q - 1 := by
      have h1 : (q - 1) * α ≤ x.val := by rw [← hm]; omega
      have h2 : q - 1 ≤ x.val / α := (Nat.le_div_iff_mul_le (by omega)).mpr h1
      omega
    have hc : ((x.val / α : ℕ) : ZMod q) = -1 := by
      rw [hb, Nat.cast_sub (by omega), ZMod.natCast_self]; ring
    have hc2 : -((x.val / α : ℕ) : ZMod q) - 2 = -1 := by rw [hc]; ring
    rw [if_pos hc2, map_one, Perm.one_apply, hc2, hc]

include hk in
theorem block_rk (x : Fin n) :
    blockEquiv hn (rev n k x) = wreathPerm α q (wK α q) (blockEquiv hn x) := by
  have hx := x.isLt
  have hxα := Nat.mod_lt x.val (show 0 < α by omega)
  rw [blockEquiv_apply, blockEquiv_apply, wreathPerm_apply]
  simp only [wK, map_one, Perm.one_apply]
  by_cases hxk : x.val < k
  · have hv : (rev n k x).val = α - 2 - x.val := by
      have hkn : k ≤ n := by
        have : 1 * α ≤ q * α := Nat.mul_le_mul_right α (by omega)
        omega
      rw [rev_val, ite_eq_left ⟨hxk, hkn⟩]; omega
    have h0 : x.val / α = 0 := Nat.div_eq_of_lt (by omega)
    have h0' : (α - 2 - x.val) / α = 0 := Nat.div_eq_of_lt (by omega)
    rw [hv, h0, h0', Nat.mod_eq_of_lt (show α - 2 - x.val < α by omega), Nat.cast_zero,
      Pi.mulSingle_eq_same, dihedralPerm_sr, Nat.mod_eq_of_lt (show x.val < α by omega),
      cast_sub_sub α _ 2 (by omega)]
    push_cast; rfl
  · have hv : rev n k x = x := Fin.ext (by rw [rev_val, ite_eq_right (by omega)])
    rw [hv]
    by_cases hxa : x.val < α
    · have hx1 : x.val = α - 1 := by omega
      have h0 : x.val / α = 0 := Nat.div_eq_of_lt hxa
      rw [h0, Nat.cast_zero, Pi.mulSingle_eq_same, dihedralPerm_sr, Nat.mod_eq_of_lt hxa, hx1,
        Nat.cast_sub (by omega), ZMod.natCast_self]
      ext <;> simp only; ring
    · have hxn : x.val < q * α := by rw [← hn]; exact hx
      have h1 : 1 ≤ x.val / α := (Nat.le_div_iff_mul_le (by omega)).mpr (by omega)
      have h2 : x.val / α < q := (Nat.div_lt_iff_lt_mul (by omega)).mpr hxn
      rw [Pi.mulSingle_eq_of_ne (natCast_ne_zero_of_lt (M := q) (a := x.val / α) (by omega) h2),
        map_one, Perm.one_apply]

end gens

end Riemann.Deep.LeanPrefix

