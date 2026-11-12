import Riemann.Deep.PrefixRev.Lean.CycleC
import Riemann.Deep.PrefixRev.Lean.BlocksA
import Riemann.Deep.PrefixRev.Lean.Parity

/-!
# Necessary conditions on the line, part B

* B.2: n even, gcd(m+1, k+1) = 1: the parity classes of the cyclic coordinate phi are two
  blocks (r_n flips the parity, r_k and r_m = r_n r_k P keep it); G ≠ S_n.
* Parity: n ≡ 1 (mod 4) and 4 | k ⇒ G ≤ A_n.
* Theorem A, necessity on the line.
-/

namespace PrefixRev

open Equiv

theorem rm_eq {n m k : ℕ} : rev n m = rev n n * rev n k * P n m k := by
  unfold P
  rw [show rev n n * rev n k * (rev n k * rev n n * rev n m)
      = rev n n * (rev n k * rev n k) * rev n n * rev n m by simp only [mul_assoc],
    rev_mul_self, mul_one, rev_mul_self, one_mul]

/-- B.2: n even on the line with gcd = 1 ⇒ G ≠ S_n. -/
theorem G_ne_top_of_even {n m k θ : ℕ} (hk2 : 2 ≤ k) (hline : n = m + k + 1) (hkm : k < m)
    (hθ : (k + 1) * θ % (n + 1) = 1) (hev : n % 2 = 0) : G n m k ≠ ⊤ := by
  have hle : G n m k ≤ fibreStab (fun x : Fin n => phi n θ x % 2) := by
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
  intro htop
  -- the points x0 = m (phi = n), x1 = P x0 (phi = 1), x2 = P x1 (phi = 2)
  have hn4 : 4 ≤ n := by omega
  let x0 : Fin n := ⟨m, by omega⟩
  have hx0 : phi n θ x0 = n := by
    have h := phi_P hline hkm hθ x0
    have h1 : phi n θ (P n m k x0) = 1 := by
      have : (P n m k x0).val = k := by
        rw [P_val hline hkm x0, ite_eq_right (show ¬ m < m by omega), ite_eq_left rfl]
      unfold phi; rw [this, Nat.mul_comm, hθ]
    have h2 := phi_le (θ := θ) x0
    have h3 := phi_ne_zero hθ x0
    rw [h1] at h
    split_ifs at h <;> omega
  have hx1 : phi n θ (P n m k x0) = 1 := by
    rw [phi_P hline hkm hθ x0, hx0, ite_eq_right (by omega)]
  have hx2 : phi n θ (P n m k (P n m k x0)) = 2 := by
    rw [phi_P hline hkm hθ, hx1, ite_eq_left (by omega)]
  have hs := hle (show swap (P n m k x0) (P n m k (P n m k x0)) ∈ G n m k by rw [htop]; trivial)
  have h01 : x0 ≠ P n m k x0 := fun e => by rw [← e, hx0] at hx1; omega
  have h02 : x0 ≠ P n m k (P n m k x0) := fun e => by rw [← e, hx0] at hx2; omega
  have := (hs x0 (P n m k (P n m k x0))).mp (by show phi n θ x0 % 2 = _ % 2; rw [hx0, hx2]; omega)
  change phi n θ (swap _ _ x0) % 2 = phi n θ (swap _ _ _) % 2 at this
  rw [swap_apply_of_ne_of_ne h01 h02, swap_apply_right, hx0, hx1] at this
  omega

/-- n ≡ 1 (mod 4), 4 | k on the line ⇒ G ≤ A_n, so G ≠ S_n. -/
theorem G_le_alt_of_mod4 {n m k : ℕ} (hline : n = m + k + 1) (_hkm : k < m)
    (h4 : n % 4 = 1 ∧ k % 4 = 0) : ∀ g ∈ G n m k, Perm.sign g = 1 :=
  (G_le_alt_iff (by omega) (by omega)).mpr (by omega)

theorem G_ne_top_of_mod4 {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (h4 : n % 4 = 1 ∧ k % 4 = 0) : G n m k ≠ ⊤ := by
  intro htop
  have := G_le_alt_of_mod4 hline hkm h4 (swap ⟨0, by omega⟩ ⟨1, by omega⟩) (by rw [htop]; trivial)
  rw [Perm.sign_swap (by rw [Ne, Fin.mk.injEq]; omega)] at this
  exact absurd this (by decide)

/-- Theorem A, necessity on the line n = m + k + 1 (2 ≤ k < m):
G = S_n ⇒ n odd, gcd(m+1, k+1) = 1, (k+1) ∤ m, and not (n ≡ 1 mod 4 and 4 | k). -/
theorem line_necessary {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m) (hline : n = m + k + 1)
    (htop : G n m k = ⊤) :
    n % 2 = 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧ ¬ (n % 4 = 1 ∧ k % 4 = 0) := by
  have hg : Nat.Coprime (k + 1) (m + 1) := by
    by_contra h
    exact G_ne_top_of_gcd hline hkm h htop
  refine ⟨?_, Nat.Coprime.symm hg, fun hd => G_ne_top_of_dvd hk hline hkm hd htop,
    fun h4 => G_ne_top_of_mod4 hline hkm h4 htop⟩
  by_contra hodd
  obtain ⟨θ, hθ⟩ := exists_theta hline hg
  exact G_ne_top_of_even hk hline hkm hθ (by omega) htop

end PrefixRev
