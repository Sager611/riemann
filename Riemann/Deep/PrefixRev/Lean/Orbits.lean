import Riemann.Deep.PrefixRev.Lean.BlocksA
import Riemann.Deep.PrefixRev.Lean.CycleA

/-!
# Theorem B.1, orbit classes: each O_c = {x : x ≡ ±c (mod g)} is G-invariant

g = gcd(k+1, m+1), positions 1-based (x.val + 1). cls g x = c ∈ [0, g/2] with
x.val + 1 ≡ ±c (mod g). Every generator preserves cls (r_i maps x to C - x, g ∣ C = i + 1),
so G does; every class 0 ≤ c ≤ g/2 is nonempty. (That each O_c is a single orbit is not
formalized here; for g = 1 it is P_cycle.)
-/

namespace PrefixRev

open Equiv

/-- The orbit label: c ∈ [0, g/2] with x + 1 ≡ ±c (mod g) (1-based position x + 1). -/
def cls (g : ℕ) {n : ℕ} (x : Fin n) : ℕ :=
  if 2 * ((x.val + 1) % g) ≤ g then (x.val + 1) % g else g - (x.val + 1) % g

theorem sub_mod_of_dvd {g C a : ℕ} (hg : 0 < g) (hC : g ∣ C) (ha : a ≤ C) :
    (C - a) % g = (g - a % g) % g := by
  have hs := Nat.mod_lt (C - a) hg
  have hr := Nat.mod_lt a hg
  have hsum : ((C - a) % g + a % g) % g = 0 := by
    rw [← Nat.add_mod, show C - a + a = C by omega]; exact Nat.mod_eq_zero_of_dvd hC
  rcases mod_cases ((C - a) % g + a % g) g (by omega) with ⟨_, h⟩ | ⟨_, h⟩
  · have h1 : a % g = 0 := by omega
    have h2 : (C - a) % g = 0 := by omega
    rw [h1, h2, Nat.sub_zero, Nat.mod_self]
  · rw [Nat.mod_eq_of_lt (show g - a % g < g by omega)]; omega

theorem cls_rev {n g i : ℕ} (hg : 0 < g) (hi : i ≤ n) (hgi : g ∣ i + 1) (x : Fin n) :
    cls g (rev n i x) = cls g x := by
  have key : ((rev n i x).val + 1) % g = (x.val + 1) % g ∨
      ((rev n i x).val + 1) % g = (g - (x.val + 1) % g) % g := by
    rw [rev_val]
    by_cases h : x.val < i
    · right
      rw [ite_eq_left ⟨h, hi⟩, show i - 1 - x.val + 1 = (i + 1) - (x.val + 1) by omega,
        sub_mod_of_dvd hg hgi (by omega)]
    · left; rw [ite_eq_right (by omega)]
  have hr := Nat.mod_lt (x.val + 1) hg
  unfold cls
  rcases key with h | h
  · rw [h]
  · rw [h]
    rcases Nat.eq_zero_or_pos ((x.val + 1) % g) with h0 | h0
    · rw [h0, Nat.sub_zero, Nat.mod_self]
      split_ifs <;> omega
    · rw [Nat.mod_eq_of_lt (show g - (x.val + 1) % g < g by omega)]
      split_ifs <;> omega

/-- B.1 (orbit classes): every element of G preserves cls g, g = gcd(k+1, m+1). -/
theorem G_le_cls {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) :
    G n m k ≤ fibreFix (fun x : Fin n => cls (Nat.gcd (k + 1) (m + 1)) x) := by
  have hpos : 0 < Nat.gcd (k + 1) (m + 1) := Nat.gcd_pos_of_pos_left _ (by omega)
  have hdk := Nat.gcd_dvd_left (k + 1) (m + 1)
  have hdm := Nat.gcd_dvd_right (k + 1) (m + 1)
  have hdn : Nat.gcd (k + 1) (m + 1) ∣ n + 1 := by
    rw [hline, show m + k + 1 + 1 = (m + 1) + (k + 1) by omega]; exact Nat.dvd_add hdm hdk
  apply G_le <;> apply mem_fibreFix_of <;> intro x
  · exact cls_rev hpos le_rfl hdn x
  · exact cls_rev hpos (by omega) hdm x
  · exact cls_rev hpos (by omega) hdk x

/-- Every label 0 ≤ c ≤ g/2 occurs. -/
theorem cls_surj {n m k : ℕ} (hline : n = m + k + 1) (c : ℕ)
    (hc : 2 * c ≤ Nat.gcd (k + 1) (m + 1)) :
    ∃ x : Fin n, cls (Nat.gcd (k + 1) (m + 1)) x = c := by
  have hpos : 0 < Nat.gcd (k + 1) (m + 1) := Nat.gcd_pos_of_pos_left _ (by omega)
  have hle : Nat.gcd (k + 1) (m + 1) ≤ k + 1 := Nat.le_of_dvd (by omega) (Nat.gcd_dvd_left _ _)
  rcases Nat.eq_zero_or_pos c with h0 | h0
  · refine ⟨⟨Nat.gcd (k + 1) (m + 1) - 1, by omega⟩, ?_⟩
    unfold cls
    rw [show Nat.gcd (k + 1) (m + 1) - 1 + 1 = Nat.gcd (k + 1) (m + 1) by omega, Nat.mod_self,
      ite_eq_left (by omega), h0]
  · refine ⟨⟨c - 1, by omega⟩, ?_⟩
    unfold cls
    rw [show c - 1 + 1 = c by omega, Nat.mod_eq_of_lt (by omega), ite_eq_left hc]

end PrefixRev
