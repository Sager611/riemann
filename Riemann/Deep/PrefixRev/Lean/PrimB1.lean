import Riemann.Deep.PrefixRev.Lean.PrimA

/-!
# Lemma 3, part B: block structure and (★_d); Prop 5 reduced to arithmetic

* block_classes (Lemma 3): every block of G is trivial or a residue class of psi mod d,
  1 < d < n, d ∣ n.
* star_of_block: then (★_d) holds (r_k maps the class c ∩ E to the class -c).
* Prop5Core (named, OPEN OBLIGATION: Steps 1-4 of the NOTE's Proposition 5, pure arithmetic):
  (★_d) for some 1 < d < n, d ∣ n ⇒ (k + 1) ∣ m.  prop5_of_core: Prop5Core ⇒ Prop5Primitive.
-/

namespace PrefixRev

open Equiv

theorem dvd_shift_iff {n d a b : ℕ} (hd : d ∣ n) (hb : b ≤ a + n) :
    d ∣ (a + n - b) % n ↔ a % d = b % d := by
  rw [Nat.dvd_mod_iff hd]
  have hn : n % d = 0 := Nat.mod_eq_zero_of_dvd hd
  have e1 : (a + n) % d = a % d := by rw [Nat.add_mod, hn, Nat.add_zero, Nat.mod_mod]
  constructor
  · intro h
    have h2 : (a + n - b) % d = 0 := Nat.mod_eq_zero_of_dvd h
    have : (a + n) % d = b % d := by
      rw [show a + n = (a + n - b) + b by omega, Nat.add_mod, h2, Nat.zero_add, Nat.mod_mod]
    rw [← e1, this]
  · intro h
    apply Nat.dvd_of_mod_eq_zero
    apply Nat.sub_mod_eq_zero_of_mod_eq
    rw [e1, h]

theorem mod_cases3 (v d : ℕ) (h : v < 3 * d) :
    (v < d ∧ v % d = v) ∨ (d ≤ v ∧ v < 2 * d ∧ v % d = v - d) ∨ (2 * d ≤ v ∧ v % d = v - 2 * d) := by
  by_cases h1 : v < d
  · exact Or.inl ⟨h1, Nat.mod_eq_of_lt h1⟩
  · by_cases h2 : v < 2 * d
    · refine Or.inr (Or.inl ⟨by omega, h2, ?_⟩)
      rw [Nat.mod_eq_sub_mod (by omega), Nat.mod_eq_of_lt (by omega)]
    · refine Or.inr (Or.inr ⟨by omega, ?_⟩)
      rw [Nat.mod_eq_sub_mod (by omega), Nat.mod_eq_sub_mod (by omega), Nat.mod_eq_of_lt (by omega)]
      omega

/-- (★_d) in the 0-based cyclic coordinate psi = phi - 1: each class psi ≡ c (mod d) with
c ≠ 0 lies inside E = {x : x < k} or misses it. (NOTE: classes C_c, c ≢ 1, 1-based.) -/
def StarD (n k θ d : ℕ) : Prop :=
  ∀ x y : Fin n, (phi n θ x - 1) % d = (phi n θ y - 1) % d → (phi n θ x - 1) % d ≠ 0 →
    (x.val < k ↔ y.val < k)

section blocks

variable {n m k θ : ℕ} (hline : n = m + k + 1) (hkm : k < m) (hθ : (k + 1) * θ % (n + 1) = 1)
include hline hkm hθ

/-- Lemma 3: a block of G through x0 is a singleton, everything, or a residue class of psi. -/
theorem block_classes {B : Set (Fin n)} (hB : IsBlockOf (G n m k) B) {x0 : Fin n}
    (hx0 : x0 ∈ B) :
    B.Subsingleton ∨ B = Set.univ ∨
      ∃ d, 1 < d ∧ d < n ∧ d ∣ n ∧
        ∀ y, y ∈ B ↔ (phi n θ y - 1) % d = (phi n θ x0 - 1) % d := by
  classical
  by_cases hs : B.Subsingleton
  · exact Or.inl hs
  right
  have hpos : 0 < n := by omega
  have hex : ∃ t, (0 < t ∧ t < n) ∧ (P n m k ^ t) x0 ∈ B := by
    obtain ⟨y, hy, z, hz, hyz⟩ := Set.not_subsingleton_iff.mp hs
    obtain ⟨w, hw, hwx⟩ : ∃ w ∈ B, w ≠ x0 := by
      by_cases h : y = x0
      · exact ⟨z, hz, fun e => hyz (h.trans e.symm)⟩
      · exact ⟨y, hy, h⟩
    refine ⟨(phi n θ w - 1 + n - (phi n θ x0 - 1)) % n, ⟨?_, Nat.mod_lt _ hpos⟩, ?_⟩
    · rcases Nat.eq_zero_or_pos ((phi n θ w - 1 + n - (phi n θ x0 - 1)) % n) with h0 | h0
      · exfalso
        have := P_pow_tOf hline hkm hθ x0 w
        rw [h0, pow_zero, Perm.one_apply] at this
        exact hwx this.symm
      · exact h0
    · rw [P_pow_tOf hline hkm hθ x0 w]; exact hw
  let d := Nat.find hex
  have hd : (0 < d ∧ d < n) ∧ (P n m k ^ d) x0 ∈ B := Nat.find_spec hex
  have hmin : ∀ t, 0 < t → t < n → (P n m k ^ t) x0 ∈ B → d ≤ t :=
    fun t h1 h2 h3 => Nat.find_min' hex ⟨⟨h1, h2⟩, h3⟩
  have hdvd : ∀ t, (P n m k ^ t) x0 ∈ B → d ∣ t := by
    intro t ht
    have h1 := T_mod hline hkm hθ hB hx0 hd.2 ht
    apply Nat.dvd_of_mod_eq_zero
    by_contra h2
    have h3 := hmin (t % d) (by omega) (by have := Nat.mod_lt t hd.1.1; omega) h1
    have := Nat.mod_lt t hd.1.1
    omega
  have hdn : d ∣ n := hdvd n (by rw [P_pow_n hline hkm hθ]; exact hx0)
  have hmem : ∀ y, y ∈ B ↔ (phi n θ y - 1) % d = (phi n θ x0 - 1) % d := by
    intro y
    have h1 := phi_ne_zero hθ x0; have h2 : phi n θ x0 ≤ n := phi_le x0
    rw [← dvd_shift_iff hdn (by omega)]
    constructor
    · intro hy
      apply hdvd
      rw [P_pow_tOf hline hkm hθ x0 y]; exact hy
    · intro h
      obtain ⟨q, hq⟩ := h
      have := T_mul hline hkm hθ hB hx0 hd.2 q
      rw [← hq, P_pow_tOf hline hkm hθ x0 y] at this
      exact this
  rcases Nat.lt_or_ge 1 d with h1 | h1
  · exact Or.inr ⟨d, h1, hd.1.2, hdn, hmem⟩
  · left
    apply Set.eq_univ_of_forall
    intro y
    rw [hmem y, show d = 1 by omega, Nat.mod_one, Nat.mod_one]

end blocks

end PrefixRev
