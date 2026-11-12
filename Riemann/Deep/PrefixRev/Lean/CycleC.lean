import Riemann.Deep.PrefixRev.Lean.CycleB

/-!
# Theorem B.1, part C: P is an n-cycle when gcd(m+1, k+1) = 1

In the 0-based coordinate psi = phi - 1, P is x ↦ x + 1 mod n; so P^t moves any point to
any point.
-/

namespace PrefixRev

open Equiv

section coords

variable {n m k θ : ℕ} (hline : n = m + k + 1) (hkm : k < m) (hθ : (k + 1) * θ % (n + 1) = 1)
include hline hkm hθ

theorem psi_P (x : Fin n) : phi n θ (P n m k x) - 1 = (phi n θ x - 1 + 1) % n := by
  have h1 := phi_ne_zero hθ x
  have h2 : phi n θ x ≤ n := phi_le x
  rw [phi_P hline hkm hθ x]
  split_ifs with h
  · rw [Nat.mod_eq_of_lt (by omega)]; omega
  · rw [show phi n θ x - 1 + 1 = n by omega, Nat.mod_self]

theorem psi_P_pow (x : Fin n) (t : ℕ) :
    phi n θ ((P n m k ^ t) x) - 1 = (phi n θ x - 1 + t) % n := by
  induction t with
  | zero =>
    have h1 := phi_ne_zero hθ x
    have h2 : phi n θ x ≤ n := phi_le x
    rw [pow_zero, Perm.one_apply, Nat.add_zero, Nat.mod_eq_of_lt (by omega)]
  | succ t ih =>
    rw [pow_succ', Perm.mul_apply, psi_P hline hkm hθ, ih, Nat.mod_add_mod, Nat.add_assoc]

/-- Every point reaches every point under powers of P. -/
theorem P_transitive (x y : Fin n) : ∃ t : ℕ, (P n m k ^ t) x = y := by
  have h1 := phi_ne_zero hθ x
  have h2 : phi n θ x ≤ n := phi_le x
  have h3 := phi_ne_zero hθ y
  have h4 : phi n θ y ≤ n := phi_le y
  refine ⟨phi n θ y - 1 + n - (phi n θ x - 1), phi_inj hθ ?_⟩
  have h5 := psi_P_pow hline hkm hθ x (phi n θ y - 1 + n - (phi n θ x - 1))
  rw [show phi n θ x - 1 + (phi n θ y - 1 + n - (phi n θ x - 1)) = phi n θ y - 1 + n by omega,
    Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)] at h5
  have h6 := phi_ne_zero hθ ((P n m k ^ (phi n θ y - 1 + n - (phi n θ x - 1))) x)
  omega

end coords

/-- θ exists when gcd(m+1, k+1) = 1 (then gcd(k+1, n+1) = 1 on the line). -/
theorem exists_theta {n m k : ℕ} (hline : n = m + k + 1) (hg : Nat.Coprime (k + 1) (m + 1)) :
    ∃ θ, (k + 1) * θ % (n + 1) = 1 := by
  have hc : Nat.Coprime (k + 1) (n + 1) := by
    rw [show n + 1 = (m + 1) + (k + 1) by omega]
    show Nat.gcd (k + 1) ((m + 1) + (k + 1)) = 1
    rw [Nat.gcd_add_self_right]; exact hg
  obtain ⟨θ, -, h⟩ := Nat.exists_mul_mod_eq_one_of_coprime hc (by omega)
  exact ⟨θ, h⟩

/-- Theorem B.1 (n-cycle): on the line n = m + k + 1 with gcd(m+1, k+1) = 1, the powers of
P = r_k r_n r_m move every point to every point, so P is an n-cycle. -/
theorem P_cycle {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (hg : Nat.Coprime (k + 1) (m + 1)) (x y : Fin n) : ∃ t : ℕ, (P n m k ^ t) x = y := by
  obtain ⟨θ, hθ⟩ := exists_theta hline hg
  exact P_transitive hline hkm hθ x y

end PrefixRev
