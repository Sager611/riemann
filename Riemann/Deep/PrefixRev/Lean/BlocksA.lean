import Riemann.Deep.PrefixRev.Lean.Defs

/-!
# Necessary conditions on the line, part A: invariant sets and interval blocks

* g = gcd(k+1, m+1) ≥ 2: the positions x (1-based) with g | x form a G-invariant set
  (each generator maps x to x or to C - x with C ∈ {n+1, m+1, k+1}, all ≡ 0 mod g); G is
  intransitive (Theorem B.1; necessity of g = 1 is the Lemma of 2610.30892, TeX l. 277).
* α = k + 1 divides n: the intervals [iα, (i+1)α) (0-based) are blocks (Lemma 6, B.3).
-/

namespace PrefixRev

open Equiv

/-- A generator that is an involution and acts on the fibres of b through a map F
preserves the partition into fibres. -/
theorem mem_fibreStab_of_invol {n : ℕ} {β : Type} (b : Fin n → β) (r : Perm (Fin n))
    (hr : ∀ x, r (r x) = x) (F : β → β) (hF : ∀ x, b (r x) = F (b x)) : r ∈ fibreStab b := by
  intro x y
  constructor
  · intro h; rw [hF, hF, h]
  · intro h
    have := congrArg F h
    rw [← hF, ← hF, hr, hr] at this
    exact this

theorem mem_fibreFix_of {n : ℕ} {β : Type} (b : Fin n → β) (r : Perm (Fin n))
    (h : ∀ x, b (r x) = b x) : r ∈ fibreFix b := h

/-! ## Invariant set for a common divisor of k + 1 and m + 1 -/

theorem dvd_rev_iff {n i d : ℕ} (hi : i ≤ n) (hd : d ∣ i + 1) (x : Fin n) :
    d ∣ (rev n i x).val + 1 ↔ d ∣ x.val + 1 := by
  rw [rev_val]
  by_cases h : x.val < i
  · rw [ite_eq_left ⟨h, hi⟩, show i - 1 - x.val + 1 = (i + 1) - (x.val + 1) by omega]
    constructor
    · intro h1
      have := Nat.dvd_sub hd h1
      rwa [Nat.sub_sub_self (by omega)] at this
    · intro h1; exact Nat.dvd_sub hd h1
  · rw [ite_eq_right (by omega)]

/-- If d ≥ 2 divides k + 1 and m + 1 (line n = m + k + 1), G is intransitive:
no element maps position 0 to position d - 1. -/
theorem G_intransitive {n m k d : ℕ} (hline : n = m + k + 1) (hkm : k < m) (hd : 2 ≤ d)
    (hdk : d ∣ k + 1) (hdm : d ∣ m + 1) :
    ∀ g ∈ G n m k, ∀ x y : Fin n, x.val = 0 → y.val = d - 1 → g x ≠ y := by
  have hdn : d ∣ n + 1 := by
    rw [hline, show m + k + 1 + 1 = (m + 1) + (k + 1) by omega]; exact Nat.dvd_add hdm hdk
  have hle : G n m k ≤ fibreFix (fun x : Fin n => decide (d ∣ x.val + 1)) := by
    apply G_le <;> apply mem_fibreFix_of <;> intro x <;> simp only [decide_eq_decide]
    · exact dvd_rev_iff le_rfl hdn x
    · exact dvd_rev_iff (by omega) hdm x
    · exact dvd_rev_iff (by omega) hdk x
  intro g hg x y hx hy h
  have := hle hg x
  rw [h] at this
  simp only [decide_eq_decide] at this
  have h2 : d ∣ y.val + 1 := by rw [hy, Nat.sub_add_cancel (by omega)]
  have h3 := this.mp h2
  rw [hx] at h3
  have := Nat.le_of_dvd (by omega) h3
  omega

/-- g = gcd(k+1, m+1) ≠ 1 ⇒ G ≠ S_n. -/
theorem G_ne_top_of_gcd {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (hg : Nat.gcd (k + 1) (m + 1) ≠ 1) : G n m k ≠ ⊤ := by
  have hpos : 0 < Nat.gcd (k + 1) (m + 1) := Nat.gcd_pos_of_pos_left _ (by omega)
  intro htop
  have hdk := Nat.gcd_dvd_left (k + 1) (m + 1)
  have hdm := Nat.gcd_dvd_right (k + 1) (m + 1)
  have hle := Nat.le_of_dvd (by omega) hdk
  have := G_intransitive hline hkm (by omega) hdk hdm
    (swap ⟨0, by omega⟩ ⟨Nat.gcd (k + 1) (m + 1) - 1, by omega⟩) (by rw [htop]; trivial)
    ⟨0, by omega⟩ ⟨Nat.gcd (k + 1) (m + 1) - 1, by omega⟩ rfl rfl
  exact this (swap_apply_left _ _)

/-! ## Interval blocks when α = k + 1 divides n -/

theorem div_reflect {a M x : ℕ} (ha : 0 < a) (hx : x < M * a) :
    (M * a - 1 - x) / a = M - 1 - x / a := by
  have h1 := Nat.div_add_mod x a
  have h2 := Nat.mod_lt x ha
  have h3 : x / a < M := (Nat.div_lt_iff_lt_mul ha).mpr hx
  obtain ⟨t, ht⟩ : ∃ t, M = x / a + 1 + t := ⟨M - x / a - 1, by omega⟩
  have e2 : (x / a + 1 + t) * a = a * (x / a) + a + t * a := by
    rw [Nat.add_mul, Nat.add_mul, Nat.one_mul, Nat.mul_comm (x / a) a]
  have e : M * a - 1 - x = (a - 1 - x % a) + t * a := by rw [ht, e2]; omega
  rw [e, Nat.add_mul_div_right _ _ ha, Nat.div_eq_of_lt (by omega)]
  omega

/-- A prefix reversal r_i with i = M a (a whole number of intervals of length a) acts on
interval indices by i ↦ M - 1 - i on the first M intervals and fixes the others. -/
theorem rev_div {n a M : ℕ} (ha : 0 < a) (hM : M * a ≤ n) (x : Fin n) :
    (rev n (M * a) x).val / a = if x.val / a < M then M - 1 - x.val / a else x.val / a := by
  rw [rev_val]
  by_cases h : x.val < M * a
  · rw [ite_eq_left ⟨h, hM⟩, ite_eq_left ((Nat.div_lt_iff_lt_mul ha).mpr h), div_reflect ha h]
  · rw [ite_eq_right (by omega), ite_eq_right (by rw [Nat.div_lt_iff_lt_mul ha]; omega)]

/-- Lemma 6 / B.3: if α = k + 1 divides n (equivalently divides m) then the intervals of
length α are blocks of G, and G ≠ S_n. -/
theorem G_le_intervals {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) (hdiv : (k + 1) ∣ n) :
    G n m k ≤ fibreStab (fun x : Fin n => x.val / (k + 1)) := by
  obtain ⟨q, hq⟩ := hdiv
  have hq1 : 2 ≤ q := by
    rcases q with _ | _ | q
    · rw [Nat.mul_zero] at hq; omega
    · rw [Nat.mul_one] at hq; omega
    · omega
  apply G_le
  · apply mem_fibreStab_of_invol _ _ (rev_rev n n) (fun i => if i < q then q - 1 - i else i)
    intro x
    have := rev_div (n := n) (a := k + 1) (M := q) (by omega) (by rw [hq, Nat.mul_comm]) x
    rwa [show q * (k + 1) = n by rw [hq, Nat.mul_comm]] at this
  · apply mem_fibreStab_of_invol _ _ (rev_rev n m)
      (fun i => if i < q - 1 then q - 1 - 1 - i else i)
    intro x
    have hm : m = (q - 1) * (k + 1) := by
      rw [Nat.sub_mul, Nat.one_mul, Nat.mul_comm, ← hq]; omega
    have := rev_div (n := n) (a := k + 1) (M := q - 1) (by omega) (by rw [← hm]; omega) x
    rwa [← hm] at this
  · apply mem_fibreStab_of_invol _ _ (rev_rev n k) id
    intro x
    rw [rev_val, id]
    split_ifs with h
    · rw [Nat.div_eq_of_lt (by omega), Nat.div_eq_of_lt (by omega)]
    · rfl

theorem G_ne_top_of_dvd {n m k : ℕ} (hk : 2 ≤ k) (hline : n = m + k + 1) (hkm : k < m)
    (hdiv : (k + 1) ∣ m) : G n m k ≠ ⊤ := by
  have hdn : (k + 1) ∣ n := by rw [hline, Nat.add_assoc]; exact Nat.dvd_add hdiv dvd_rfl
  intro htop
  have hle := G_le_intervals hline hkm hdn
  have hs := hle (show swap (⟨0, by omega⟩ : Fin n) ⟨k + 1, by omega⟩ ∈ G n m k by
    rw [htop]; trivial)
  have := (hs ⟨0, by omega⟩ ⟨1, by omega⟩).mp
    (by simp only; rw [Nat.div_eq_of_lt (by omega), Nat.div_eq_of_lt (by omega)])
  rw [swap_apply_left, swap_apply_of_ne_of_ne (by rw [Ne, Fin.mk.injEq]; omega)
    (by rw [Ne, Fin.mk.injEq]; omega)] at this
  simp only at this
  rw [Nat.div_self (by omega), Nat.div_eq_of_lt (by omega)] at this
  omega

end PrefixRev
