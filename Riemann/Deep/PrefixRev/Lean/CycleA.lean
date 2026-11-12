import Riemann.Deep.PrefixRev.Lean.Defs
import Mathlib.Data.Int.GCD

/-!
# Theorem B.1, part A: P = r_k r_n r_m on the line n = m + k + 1 (Lemma 1)

Explicit form of P, its first-return description (rotation by α = k + 1 on Z/(n+1) minus 0),
modular helpers, and the cyclic coordinate phi x = θ (x+1) mod (n+1) (θ an inverse of α).
-/

namespace PrefixRev

open Equiv

/-- P = r_k r_n r_m (r_m acts first). -/
def P (n m k : ℕ) : Perm (Fin n) := rev n k * rev n n * rev n m

/-- Lemma 1 (0-based): P x = x + α (x < m), P m = k, P x = x - β (x > m). -/
theorem P_val {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) (x : Fin n) :
    (P n m k x).val =
      if x.val < m then x.val + (k + 1) else if x.val = m then k else x.val - (m + 1) := by
  have hx := x.isLt
  simp only [P, Perm.mul_apply]
  rw [rev_val, rev_val, rev_val]
  split_ifs <;> omega

theorem mod_cases (s N : ℕ) (h : s < 2 * N) : (s < N ∧ s % N = s) ∨ (N ≤ s ∧ s % N = s - N) := by
  by_cases h1 : s < N
  · exact Or.inl ⟨h1, Nat.mod_eq_of_lt h1⟩
  · refine Or.inr ⟨by omega, ?_⟩
    rw [Nat.mod_eq_sub_mod (by omega), Nat.mod_eq_of_lt (by omega)]

/-- Lemma 1, first-return form. Identify x with y = x.val + 1 in Z/(n+1) minus 0 and let
R y = y + α. Then P y = R y if R y ≠ 0, and P y = R (R y) if R y = 0. -/
theorem P_firstReturn {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) (x : Fin n) :
    (P n m k x).val + 1 =
      if (x.val + 1 + (k + 1)) % (n + 1) = 0 then (x.val + 1 + 2 * (k + 1)) % (n + 1)
      else (x.val + 1 + (k + 1)) % (n + 1) := by
  have hx := x.isLt
  rw [P_val hline hkm x]
  rcases lt_trichotomy x.val m with h | h | h
  · have e1 : (x.val + 1 + (k + 1)) % (n + 1) = x.val + 1 + (k + 1) :=
      Nat.mod_eq_of_lt (by omega)
    rw [e1, ite_eq_left h, ite_eq_right (by omega)]; omega
  · have e1 : (x.val + 1 + (k + 1)) % (n + 1) = 0 := by
      rw [show x.val + 1 + (k + 1) = n + 1 by omega, Nat.mod_self]
    have e2 : (x.val + 1 + 2 * (k + 1)) % (n + 1) = k + 1 := by
      rw [show x.val + 1 + 2 * (k + 1) = (k + 1) + (n + 1) by omega, Nat.add_mod_right,
        Nat.mod_eq_of_lt (by omega)]
    rw [e1, e2, ite_eq_right (by omega), ite_eq_left h, ite_eq_left rfl]
  · have e1 : (x.val + 1 + (k + 1)) % (n + 1) = x.val - m := by
      rw [show x.val + 1 + (k + 1) = (x.val - m) + (n + 1) by omega, Nat.add_mod_right,
        Nat.mod_eq_of_lt (by omega)]
    rw [e1, ite_eq_right (by omega), ite_eq_right (by omega), ite_eq_right (by omega)]; omega

/-! Modular helpers (θ is an inverse of α modulo N). -/

theorem mulmod_one {α θ N : ℕ} (h : α * θ % N = 1) (a : ℕ) : α * (θ * a) % N = a % N := by
  rw [← Nat.mul_assoc, Nat.mul_mod, h, Nat.one_mul, Nat.mod_mod]

theorem theta_cancel {α θ N : ℕ} (h : α * θ % N = 1) {a b : ℕ} (hab : θ * a % N = θ * b % N) :
    a % N = b % N := by
  rw [← mulmod_one h a, ← mulmod_one h b, Nat.mul_mod, hab, ← Nat.mul_mod]

theorem theta_congr {θ N a b : ℕ} (hab : a % N = b % N) : θ * a % N = θ * b % N := by
  rw [Nat.mul_mod, hab, ← Nat.mul_mod]

theorem theta_add {α θ N : ℕ} (h : α * θ % N = 1) (a : ℕ) :
    θ * (a + α) % N = (θ * a % N + 1) % N := by
  rw [Nat.mul_add, Nat.add_mod, Nat.mul_comm θ α, h]

/-- Cyclic coordinate (1-based, as in the NOTE): phi x = θ (x+1) mod (n+1). -/
def phi (n θ : ℕ) (x : Fin n) : ℕ := θ * (x.val + 1) % (n + 1)

theorem phi_ne_zero {n θ α : ℕ} (hθ : α * θ % (n + 1) = 1) (x : Fin n) : phi n θ x ≠ 0 := by
  intro h0
  have hx := x.isLt
  have : (x.val + 1) % (n + 1) = 0 % (n + 1) :=
    theta_cancel hθ (by rw [Nat.mul_zero, Nat.zero_mod]; exact h0)
  rw [Nat.mod_eq_of_lt (by omega), Nat.zero_mod] at this
  omega

theorem phi_le {n θ : ℕ} (x : Fin n) : phi n θ x ≤ n := by
  have := Nat.mod_lt (θ * (x.val + 1)) (show 0 < n + 1 by omega)
  unfold phi; omega

theorem phi_inj {n θ α : ℕ} (hθ : α * θ % (n + 1) = 1) {x y : Fin n}
    (h : phi n θ x = phi n θ y) : x = y := by
  have hx := x.isLt; have hy := y.isLt
  have := theta_cancel hθ h
  rw [Nat.mod_eq_of_lt (by omega), Nat.mod_eq_of_lt (by omega)] at this
  exact Fin.ext (by omega)

end PrefixRev
