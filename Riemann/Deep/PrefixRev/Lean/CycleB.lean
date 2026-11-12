import Riemann.Deep.PrefixRev.Lean.CycleA

/-!
# Theorem B.1, part B: Lemma 2 (cyclic coordinates)

With θ an inverse of α = k + 1 modulo N = n + 1 and phi x = θ (x+1) mod N (values 1..n):
P is the shift j ↦ j + 1 (j < n), n ↦ 1; r_n is j ↦ n + 1 - j; r_k is j ↦ n + 2 - j on
E = phi [0, k) and the identity off E.
-/

namespace PrefixRev

open Equiv

section coords

variable {n m k θ : ℕ} (hline : n = m + k + 1) (hkm : k < m) (hθ : (k + 1) * θ % (n + 1) = 1)
include hline hkm hθ

/-- Lemma 2, P is the shift: j ↦ j + 1 (j < n), n ↦ 1. -/
theorem phi_P (x : Fin n) :
    phi n θ (P n m k x) = if phi n θ x < n then phi n θ x + 1 else 1 := by
  have hx := x.isLt
  have hPv := P_val hline hkm x
  have hne := phi_ne_zero hθ (P n m k x)
  have hle : phi n θ x ≤ n := phi_le x
  have hpos := phi_ne_zero hθ x
  by_cases hxm : x.val = m
  · have h1 : phi n θ (P n m k x) = 1 := by
      unfold phi
      rw [hPv, ite_eq_right (by omega), ite_eq_left hxm, Nat.mul_comm, hθ]
    have h2 : (phi n θ x + 1) % (n + 1) = 0 := by
      unfold phi
      rw [← theta_add hθ, hxm, show m + 1 + (k + 1) = n + 1 by omega, Nat.mul_mod_left]
    have h3 : phi n θ x = n := by
      rcases mod_cases (phi n θ x + 1) (n + 1) (by omega) with ⟨_, h4⟩ | ⟨_, h4⟩ <;> omega
    rw [h1, ite_eq_right (by omega)]
  · have hc : ((P n m k x).val + 1) % (n + 1) = (x.val + 1 + (k + 1)) % (n + 1) := by
      rw [hPv]
      by_cases hxl : x.val < m
      · rw [ite_eq_left hxl]; congr 1; omega
      · rw [ite_eq_right hxl, ite_eq_right hxm, Nat.mod_eq_sub_mod (show n + 1 ≤ x.val + 1 + (k + 1) by omega)]
        congr 1; omega
    have h1 : phi n θ (P n m k x) = (phi n θ x + 1) % (n + 1) := by
      unfold phi
      rw [theta_congr hc, theta_add hθ]
    rcases mod_cases (phi n θ x + 1) (n + 1) (by omega) with ⟨h5, h4⟩ | ⟨h5, h4⟩
    · rw [h1, h4, ite_eq_left (by omega)]
    · exfalso; rw [h1, h4] at hne; omega

/-- Lemma 2, r_n is the reflection j ↦ n + 1 - j. -/
theorem phi_rn (x : Fin n) : phi n θ (rev n n x) = n + 1 - phi n θ x := by
  have hx := x.isLt
  have hv : (rev n n x).val = n - 1 - x.val := by rw [rev_val, ite_eq_left (by omega)]
  have hs : (phi n θ (rev n n x) + phi n θ x) % (n + 1) = 0 := by
    unfold phi
    rw [← Nat.add_mod, ← Nat.mul_add, hv, show n - 1 - x.val + 1 + (x.val + 1) = n + 1 by omega,
      Nat.mul_mod_left]
  have h1 := phi_ne_zero hθ (rev n n x); have h2 := phi_ne_zero hθ x
  have h3 : phi n θ (rev n n x) ≤ n := phi_le _
  have h4 : phi n θ x ≤ n := phi_le x
  rcases mod_cases (phi n θ (rev n n x) + phi n θ x) (n + 1) (by omega) with ⟨_, h5⟩ | ⟨_, h5⟩ <;>
    omega

/-- Lemma 2, r_k is j ↦ n + 2 - j on E = phi [0, k) and the identity off E. -/
theorem phi_rk (x : Fin n) :
    phi n θ (rev n k x) = if x.val < k then n + 2 - phi n θ x else phi n θ x := by
  have hx := x.isLt
  by_cases hxk : x.val < k
  · rw [ite_eq_left hxk]
    have hv : (rev n k x).val = k - 1 - x.val := by rw [rev_val, ite_eq_left (by omega)]
    have hs : (phi n θ (rev n k x) + phi n θ x) % (n + 1) = 1 := by
      unfold phi
      rw [← Nat.add_mod, ← Nat.mul_add, hv, show k - 1 - x.val + 1 + (x.val + 1) = k + 1 by omega,
        Nat.mul_comm, hθ]
    have h1 := phi_ne_zero hθ (rev n k x); have h2 := phi_ne_zero hθ x
    have h3 : phi n θ (rev n k x) ≤ n := phi_le _
    have h4 : phi n θ x ≤ n := phi_le x
    rcases mod_cases (phi n θ (rev n k x) + phi n θ x) (n + 1) (by omega) with ⟨_, h5⟩ | ⟨_, h5⟩ <;>
      omega
  · rw [ite_eq_right hxk]
    have : rev n k x = x := Fin.ext (by rw [rev_val, ite_eq_right (by omega)])
    rw [this]

end coords

end PrefixRev
