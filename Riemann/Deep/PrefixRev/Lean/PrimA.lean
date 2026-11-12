import Riemann.Deep.PrefixRev.Lean.Main

/-!
# Lemma 3, part A: blocks of G are residue classes of the cyclic coordinate

On the line with θ an inverse of α = k + 1 mod n + 1, let psi = phi - 1 (values 0..n-1);
P acts as psi ↦ psi + 1 mod n. A block B of G containing x0 is a singleton, everything, or
{y : psi y ≡ psi x0 (mod d)} for a divisor 1 < d < n of n.
-/

namespace PrefixRev

open Equiv

theorem block_eq_of_mem {n : ℕ} {H : Subgroup (Perm (Fin n))} {B : Set (Fin n)}
    (hB : IsBlockOf H B) {g : Perm (Fin n)} (hg : g ∈ H) {x : Fin n} (hx : x ∈ B)
    (hgx : g x ∈ B) : ⇑g '' B = B := by
  rcases hB g hg with h | h
  · exact h
  · exact absurd hgx (Set.disjoint_left.mp h ⟨x, hx, rfl⟩)

theorem block_mem_of {n : ℕ} {H : Subgroup (Perm (Fin n))} {B : Set (Fin n)}
    (hB : IsBlockOf H B) {g : Perm (Fin n)} (hg : g ∈ H) {x y : Fin n} (hx : x ∈ B)
    (hgx : g x ∈ B) (hy : y ∈ B) : g y ∈ B := by
  rw [← block_eq_of_mem hB hg hx hgx]; exact ⟨y, hy, rfl⟩

section blocks

variable {n m k θ : ℕ} (hline : n = m + k + 1) (hkm : k < m) (hθ : (k + 1) * θ % (n + 1) = 1)
include hline hkm hθ

theorem P_pow_eq_iff (x y : Fin n) (t : ℕ) :
    (P n m k ^ t) x = y ↔ (phi n θ x - 1 + t) % n = phi n θ y - 1 := by
  constructor
  · intro h; rw [← h, psi_P_pow hline hkm hθ]
  · intro h
    apply phi_inj hθ
    have := psi_P_pow hline hkm hθ x t
    rw [h] at this
    have h1 := phi_ne_zero hθ ((P n m k ^ t) x)
    have h2 := phi_ne_zero hθ y
    omega

theorem P_pow_n (x : Fin n) : (P n m k ^ n) x = x := by
  rw [P_pow_eq_iff hline hkm hθ, Nat.add_mod_right, Nat.mod_eq_of_lt]
  have h1 := phi_ne_zero hθ x
  have h2 : phi n θ x ≤ n := phi_le x
  omega

/-- The orbit map: y = P^(tOf x0 y) x0. -/
theorem P_pow_tOf (x0 y : Fin n) :
    (P n m k ^ ((phi n θ y - 1 + n - (phi n θ x0 - 1)) % n)) x0 = y := by
  rw [P_pow_eq_iff hline hkm hθ, Nat.add_mod_mod]
  have h1 := phi_ne_zero hθ x0; have h2 : phi n θ x0 ≤ n := phi_le x0
  have h3 := phi_ne_zero hθ y; have h4 : phi n θ y ≤ n := phi_le y
  rw [show phi n θ x0 - 1 + (phi n θ y - 1 + n - (phi n θ x0 - 1)) = phi n θ y - 1 + n by omega,
    Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]

variable {B : Set (Fin n)} (hB : IsBlockOf (G n m k) B) {x0 : Fin n} (hx0 : x0 ∈ B)
include hB hx0

theorem T_add {t u : ℕ} (ht : (P n m k ^ t) x0 ∈ B) (hu : (P n m k ^ u) x0 ∈ B) :
    (P n m k ^ (t + u)) x0 ∈ B := by
  rw [pow_add, Perm.mul_apply]
  exact block_mem_of hB ((G n m k).pow_mem (P_mem_G n m k) t) hx0 ht hu

theorem T_sub {t u : ℕ} (hut : u ≤ t) (ht : (P n m k ^ t) x0 ∈ B)
    (hu : (P n m k ^ u) x0 ∈ B) : (P n m k ^ (t - u)) x0 ∈ B := by
  have e : (P n m k ^ (t - u)) ((P n m k ^ u) x0) = (P n m k ^ t) x0 := by
    rw [← Perm.mul_apply, ← pow_add, Nat.sub_add_cancel hut]
  exact block_mem_of hB ((G n m k).pow_mem (P_mem_G n m k) (t - u)) hu (by rw [e]; exact ht) hx0

theorem T_mul {d : ℕ} (hd : (P n m k ^ d) x0 ∈ B) (q : ℕ) : (P n m k ^ (d * q)) x0 ∈ B := by
  induction q with
  | zero => rw [Nat.mul_zero, pow_zero, Perm.one_apply]; exact hx0
  | succ q ih => rw [Nat.mul_succ]; exact T_add hline hkm hθ hB hx0 ih hd

theorem T_mod {d t : ℕ} (hd : (P n m k ^ d) x0 ∈ B) (ht : (P n m k ^ t) x0 ∈ B) :
    (P n m k ^ (t % d)) x0 ∈ B := by
  have e : t % d = t - d * (t / d) := by have := Nat.mod_add_div t d; omega
  rw [e]
  exact T_sub hline hkm hθ hB hx0 (by have := Nat.mod_add_div t d; omega) ht
    (T_mul hline hkm hθ hB hx0 hd (t / d))

end blocks

end PrefixRev
