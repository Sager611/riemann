import Riemann.Deep.PrefixRev.Lean.PrimB1

/-!
# Lemma 3, part C: the r_k action forces (★_d)
-/

namespace PrefixRev

open Equiv

theorem cancel_shift {d ρ c c0 : ℕ} (hρ : ρ < d) (hc : c < d) (hc0 : c0 < d)
    (h : (ρ + (c0 + d - c)) % d = c0) : ρ = c := by
  rcases mod_cases3 (ρ + (c0 + d - c)) d (by omega) with h1 | h1 | h1 <;> omega

theorem class_neg {n d a c : ℕ} (hdn : d ∣ n) (hdodd : d % 2 = 1) (ha : a ≤ n)
    (h1 : a % d = c) (h2 : (n - a) % d = c) : c = 0 := by
  have hd : 0 < d := by omega
  have hcd : c < d := by rw [← h1]; exact Nat.mod_lt _ hd
  have hsum : ((n - a) % d + a % d) % d = 0 := by
    rw [← Nat.add_mod, show n - a + a = n by omega]
    exact Nat.mod_eq_zero_of_dvd hdn
  rw [h1, h2] at hsum
  rcases mod_cases (c + c) d (by omega) with ⟨_, h⟩ | ⟨_, h⟩ <;> omega

theorem two_d_le {n d : ℕ} (hodd : n % 2 = 1) (hd1 : 1 < d) (hdlt : d < n) (hdn : d ∣ n) :
    2 * d ≤ n := by
  obtain ⟨q, hq⟩ := hdn
  rcases q with _ | _ | q
  · omega
  · rw [Nat.mul_one] at hq; omega
  · rw [hq, Nat.mul_add, Nat.mul_add]; omega

section blocks

variable {n m k θ : ℕ} (hline : n = m + k + 1) (hkm : k < m) (hθ : (k + 1) * θ % (n + 1) = 1)
include hline hkm hθ

theorem psi_shift_mod {d : ℕ} (hdn : d ∣ n) (z : Fin n) (s : ℕ) :
    (phi n θ ((P n m k ^ s) z) - 1) % d = ((phi n θ z - 1) % d + s) % d := by
  rw [psi_P_pow hline hkm hθ, Nat.mod_mod_of_dvd _ hdn, Nat.mod_add_mod]

theorem psi_rk {x : Fin n} (hx : x.val < k) : phi n θ (rev n k x) - 1 = n - (phi n θ x - 1) := by
  have h1 := phi_rk hline hkm hθ x
  rw [ite_eq_left hx] at h1
  have h2 := phi_ne_zero hθ x
  have h3 : phi n θ x ≤ n := phi_le x
  omega

/-- The key step of Lemma 3: in a block system of residue classes mod d (d odd), a class
c ≠ 0 cannot contain a point x of E = {x < k} and a point y outside E. -/
theorem star_key {B : Set (Fin n)} (hB : IsBlockOf (G n m k) B) {x0 : Fin n} {d : ℕ}
    (hd1 : 1 < d) (hdlt : d < n) (hdn : d ∣ n) (hdodd : d % 2 = 1) (hodd : n % 2 = 1)
    (hmem : ∀ y, y ∈ B ↔ (phi n θ y - 1) % d = (phi n θ x0 - 1) % d)
    (x y : Fin n) (hxy : (phi n θ x - 1) % d = (phi n θ y - 1) % d)
    (hc : (phi n θ x - 1) % d ≠ 0) (hx : x.val < k) (hy : ¬ y.val < k) : False := by
  obtain ⟨c, hcdef⟩ : ∃ c, c = (phi n θ x - 1) % d := ⟨_, rfl⟩
  obtain ⟨c0, hc0def⟩ : ∃ c0, c0 = (phi n θ x0 - 1) % d := ⟨_, rfl⟩
  have hcd : c < d := by rw [hcdef]; exact Nat.mod_lt _ (by omega)
  have hc0d : c0 < d := by rw [hc0def]; exact Nat.mod_lt _ (by omega)
  obtain ⟨s, hs⟩ : ∃ s, s = c0 + d - c := ⟨_, rfl⟩
  have hsn : s ≤ n := by have := two_d_le hodd hd1 hdlt hdn; omega
  have hshift : ∀ z : Fin n, (phi n θ z - 1) % d = c → (P n m k ^ s) z ∈ B := by
    intro z hz
    rw [hmem, psi_shift_mod hline hkm hθ hdn, hz, show c + s = c0 + d by omega,
      Nat.add_mod_right, Nat.mod_eq_of_lt hc0d, hc0def]
  have hx' := hshift x hcdef.symm
  have hy' := hshift y (by rw [← hxy, ← hcdef])
  have hg : P n m k ^ s * rev n k * P n m k ^ (n - s) ∈ G n m k :=
    (G n m k).mul_mem ((G n m k).mul_mem ((G n m k).pow_mem (P_mem_G n m k) s) (rk_mem n m k))
      ((G n m k).pow_mem (P_mem_G n m k) (n - s))
  have hback : ∀ z : Fin n, (P n m k ^ (n - s)) ((P n m k ^ s) z) = z := by
    intro z
    rw [← Perm.mul_apply, ← pow_add, Nat.sub_add_cancel hsn, P_pow_n hline hkm hθ]
  have hrky : rev n k y = y := Fin.ext (by rw [rev_val, ite_eq_right (by omega)])
  have hgy : (P n m k ^ s * rev n k * P n m k ^ (n - s)) ((P n m k ^ s) y) ∈ B := by
    rw [Perm.mul_apply, Perm.mul_apply, hback, hrky]; exact hy'
  have hgx := block_mem_of hB hg hy' hgy hx'
  rw [Perm.mul_apply, Perm.mul_apply, hback, hmem, psi_shift_mod hline hkm hθ hdn,
    psi_rk hline hkm hθ hx, ← hc0def, hs] at hgx
  have hρ := cancel_shift (Nat.mod_lt _ (show 0 < d by omega)) hcd hc0d hgx
  have h3 : phi n θ x ≤ n := phi_le x
  have := class_neg hdn hdodd (show phi n θ x - 1 ≤ n by omega) hcdef.symm hρ
  exact hc (by rw [← hcdef]; exact this)

/-- The r_k action forces (★_d) for the modulus d of a nontrivial block. -/
theorem star_of_block (hodd : n % 2 = 1) {B : Set (Fin n)} (hB : IsBlockOf (G n m k) B)
    {x0 : Fin n} {d : ℕ} (hd1 : 1 < d) (hdlt : d < n) (hdn : d ∣ n)
    (hmem : ∀ y, y ∈ B ↔ (phi n θ y - 1) % d = (phi n θ x0 - 1) % d) : StarD n k θ d := by
  have hdodd : d % 2 = 1 := by
    by_contra h
    have : 2 ∣ n := Nat.dvd_trans (Nat.dvd_of_mod_eq_zero (by omega)) hdn
    omega
  intro x y hxy hc
  constructor
  · intro hx; by_contra hy
    exact star_key hline hkm hθ hB hd1 hdlt hdn hdodd hodd hmem x y hxy hc hx hy
  · intro hy; by_contra hx
    exact star_key hline hkm hθ hB hd1 hdlt hdn hdodd hodd hmem y x hxy.symm (by rw [← hxy]; exact hc) hy hx

end blocks

/-- Named hypothesis (OPEN OBLIGATION: Steps 1-4 of the NOTE's Proposition 5, pure
arithmetic): on the odd line with gcd = 1, (★_d) for a divisor 1 < d < n of n forces
(k + 1) ∣ m. -/
def Prop5Core : Prop :=
  ∀ n m k θ d : ℕ, 2 ≤ k → k < m → n = m + k + 1 → n % 2 = 1 → Nat.Coprime (m + 1) (k + 1) →
    (k + 1) * θ % (n + 1) = 1 → 1 < d → d < n → d ∣ n → StarD n k θ d → (k + 1) ∣ m

/-- Lemma 3 + (★_d): Prop5Core implies the primitivity statement Prop5Primitive. -/
theorem prop5_of_core (hC : Prop5Core) : Prop5Primitive := by
  intro n m k hk hkm hline hodd hg hdiv
  obtain ⟨θ, hθ⟩ := exists_theta hline (Nat.Coprime.symm hg)
  refine ⟨fun x y => ?_, fun B hB => ?_⟩
  · obtain ⟨t, ht⟩ := P_cycle hline hkm (Nat.Coprime.symm hg) x y
    exact ⟨P n m k ^ t, (G n m k).pow_mem (P_mem_G n m k) t, ht⟩
  · rcases Set.eq_empty_or_nonempty B with h | ⟨x0, hx0⟩
    · left; rw [h]; exact Set.subsingleton_empty
    · rcases block_classes hline hkm hθ hB hx0 with h | h | ⟨d, hd1, hdn', hdn, hmem⟩
      · exact Or.inl h
      · exact Or.inr h
      · exact absurd (hC n m k θ d hk hkm hline hodd hg hθ hd1 hdn' hdn
          (star_of_block hline hkm hθ hodd hB hd1 hdn' hdn hmem)) hdiv

/-- Theorem A with Prop 5 reduced to its arithmetic core. -/
theorem theoremA_core (hOff : OffLineLZZW) (hC : Prop5Core) (hJ : JonesLemma7) {n m k : ℕ}
    (hk : 2 ≤ k) (hkm : k < m) (hmn : m < n) (hsum : m + k < n) :
    G n m k = ⊤ ↔ (n % 2 = 1 ∧ m + k = n - 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧
      ¬ (n % 4 = 1 ∧ k % 4 = 0)) :=
  theoremA hOff (prop5_of_core hC) hJ hk hkm hmn hsum

end PrefixRev
