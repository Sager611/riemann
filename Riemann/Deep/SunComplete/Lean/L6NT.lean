import Riemann.Deep.SunComplete.Lean.Sieve
import Mathlib.Data.Nat.Totient
import Mathlib.Data.Nat.PrimeFin
import Mathlib.Data.Finset.Max

/-!
# sunlean · L6NT: number theory for Lemma 6 (no sieve)

For `m ≠ 0` with `r = ω(m)` prime divisors: `φ(m) = m·∏(1 − 1/p)` (Mathlib's Euler product),
`(r + 1)·∏(1 − 1/p) ≥ 1` and `(r + 1)! ≤ ∏ p ≤ m` (distinct integers `≥ 2`).
-/

namespace SunLean

theorem primeDivisors_eq' (m : ℕ) (hm : m ≠ 0) : primeDivisors m = m.primeFactors := by
  ext p
  simp only [primeDivisors, Finset.mem_filter, Finset.mem_range, Nat.mem_primeFactors]
  constructor
  · rintro ⟨_, hp, hd⟩
    exact ⟨hp, hd, hm⟩
  · rintro ⟨hp, hd, -⟩
    exact ⟨Nat.lt_succ_of_le (Nat.le_of_dvd (Nat.pos_of_ne_zero hm) hd), hp, hd⟩

theorem two_le_of_mem_primeDivisors {m p : ℕ} (h : p ∈ primeDivisors m) : 2 ≤ p := by
  unfold primeDivisors at h
  exact (Finset.mem_filter.mp h).2.1.two_le

/-- Euler's product: `φ(m) = m·∏_{p | m} (1 − 1/p)`. -/
theorem phi_eq_prod (m : ℕ) (hm : m ≠ 0) :
    (phi m : ℚ) = m * ∏ p ∈ primeDivisors m, (1 - 1 / (p : ℚ)) := by
  rw [primeDivisors_eq' m hm]
  have h : phi m = Nat.totient m := rfl
  rw [h, Nat.totient_eq_mul_prod_factors]
  simp only [one_div]

/-- `∏_{p | m} p ≤ m`. -/
theorem prod_primeDivisors_le (m : ℕ) (hm : m ≠ 0) : ∏ p ∈ primeDivisors m, p ≤ m := by
  rw [primeDivisors_eq' m hm]
  exact Nat.le_of_dvd (Nat.pos_of_ne_zero hm) (Nat.prod_primeFactors_dvd m)

/-- For distinct integers `≥ 2`: `(|S| + 1)·∏(1 − 1/p) ≥ 1`. -/
theorem theta_card_bound (S : Finset ℕ) (hS : ∀ p ∈ S, 2 ≤ p) :
    1 ≤ ((S.card : ℚ) + 1) * ∏ p ∈ S, (1 - 1 / (p : ℚ)) := by
  induction S using Finset.induction_on_max with
  | empty => simp
  | insert a s hlt ih =>
    have ha : a ∉ s := fun h => lt_irrefl a (hlt a h)
    have hs2 : ∀ p ∈ s, 2 ≤ p := fun p hp => hS p (Finset.mem_insert_of_mem hp)
    have ha2 : 2 ≤ a := hS a (Finset.mem_insert_self a s)
    have hsub : s ⊆ Finset.Ico 2 a := fun x hx => Finset.mem_Ico.mpr ⟨hs2 x hx, hlt x hx⟩
    have hcard : s.card + 2 ≤ a := by
      have := Finset.card_le_card hsub
      rw [Nat.card_Ico] at this
      omega
    rw [Finset.prod_insert ha, Finset.card_insert_of_notMem ha]
    have ih' := ih hs2
    have hθ0 : 0 ≤ ∏ p ∈ s, (1 - 1 / (p : ℚ)) := Finset.prod_nonneg (fun p hp => by
      have h2 : (2 : ℚ) ≤ p := by exact_mod_cast hs2 p hp
      have h3 : 1 / (p : ℚ) ≤ 1 := by rw [div_le_one (by linarith)]; linarith
      linarith)
    have haQ : ((s.card : ℚ) + 2) ≤ a := by exact_mod_cast hcard
    have ha0 : (0 : ℚ) < a := by linarith
    have key : ((s.card : ℚ) + 1) ≤ ((s.card : ℚ) + 2) * (1 - 1 / (a : ℚ)) := by
      have h4 : ((s.card : ℚ) + 2) / a ≤ 1 := by rw [div_le_one ha0]; exact haQ
      have e : ((s.card : ℚ) + 2) * (1 - 1 / (a : ℚ)) = ((s.card : ℚ) + 2) - ((s.card : ℚ) + 2) / a := by
        ring
      linarith
    push_cast
    calc (1 : ℚ) ≤ ((s.card : ℚ) + 1) * ∏ p ∈ s, (1 - 1 / (p : ℚ)) := ih'
      _ ≤ (((s.card : ℚ) + 2) * (1 - 1 / (a : ℚ))) * ∏ p ∈ s, (1 - 1 / (p : ℚ)) :=
          mul_le_mul_of_nonneg_right key hθ0
      _ = ((s.card : ℚ) + 1 + 1) * ((1 - 1 / (a : ℚ)) * ∏ p ∈ s, (1 - 1 / (p : ℚ))) := by ring

/-- For distinct integers `≥ 2`: `(|S| + 1)! ≤ ∏ p`. -/
theorem fact_le_prod (S : Finset ℕ) (hS : ∀ p ∈ S, 2 ≤ p) : (S.card + 1).factorial ≤ ∏ p ∈ S, p := by
  induction S using Finset.induction_on_max with
  | empty => simp
  | insert a s hlt ih =>
    have ha : a ∉ s := fun h => lt_irrefl a (hlt a h)
    have hs2 : ∀ p ∈ s, 2 ≤ p := fun p hp => hS p (Finset.mem_insert_of_mem hp)
    have ha2 : 2 ≤ a := hS a (Finset.mem_insert_self a s)
    have hsub : s ⊆ Finset.Ico 2 a := fun x hx => Finset.mem_Ico.mpr ⟨hs2 x hx, hlt x hx⟩
    have hcard : s.card + 1 + 1 ≤ a := by
      have := Finset.card_le_card hsub
      rw [Nat.card_Ico] at this
      omega
    rw [Finset.prod_insert ha, Finset.card_insert_of_notMem ha, Nat.factorial_succ]
    exact Nat.mul_le_mul hcard (ih hs2)

end SunLean
