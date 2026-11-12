import Riemann.Deep.SunComplete.Lean.L6

/-!
# sunlean · L6b: a sharper prime-product bound narrows Lemma 6's middle range

`k` distinct primes have product `≥ 2·3·5·⋯·(2k − 1)` (at most one is even), so the no-sieve
case covers every `q` with `ω(q − 1) ≥ 12` and every `q ≥ 196·12²·4^11 = 118,380,036,096`.
-/

namespace SunLean

theorem prime_of_mem_primeDivisors {m p : ℕ} (h : p ∈ primeDivisors m) : p.Prime := by
  unfold primeDivisors at h
  exact (Finset.mem_filter.mp h).2.1

/-- `oddProd k = 2·3·5·⋯·(2k − 1)`. -/
def oddProd : ℕ → ℕ
  | 0 => 1
  | k + 1 => (if k = 0 then 2 else 2 * k + 1) * oddProd k

/-- A prime above `|s|` smaller primes is `≥ 2|s| + 1` (`≥ 2` if `s = ∅`). -/
theorem max_prime_bound (s : Finset ℕ) (a : ℕ) (hs : ∀ p ∈ s, p.Prime ∧ p < a) (ha : a.Prime) :
    (if s.card = 0 then 2 else 2 * s.card + 1) ≤ a := by
  split_ifs with h
  · exact ha.two_le
  · have hinj : (s : Set ℕ).InjOn (fun p => (p + 1) / 2) := by
      intro p hp p' hp' heq
      have hp1 := hs p (Finset.mem_coe.mp hp)
      have hp2 := hs p' (Finset.mem_coe.mp hp')
      have h1 := hp1.1.eq_two_or_odd
      have h2 := hp2.1.eq_two_or_odd
      have h3 := hp1.1.two_le
      have h4 := hp2.1.two_le
      simp only at heq
      omega
    have ha2 := ha.eq_two_or_odd
    have hmap : Set.MapsTo (fun p => (p + 1) / 2) (s : Set ℕ) (Finset.Icc 1 ((a - 1) / 2) : Set ℕ) := by
      intro p hp
      have hp1 := hs p (Finset.mem_coe.mp hp)
      have h3 := hp1.1.two_le
      have h5 := hp1.2
      show (p + 1) / 2 ∈ (Finset.Icc 1 ((a - 1) / 2) : Set ℕ)
      rw [Finset.mem_coe, Finset.mem_Icc]
      omega
    have hcard := Finset.card_le_card_of_injOn (fun p => (p + 1) / 2) hmap hinj
    rw [Nat.card_Icc] at hcard
    omega

/-- `k` distinct primes have product at least `oddProd k`. -/
theorem oddProd_le_prod (S : Finset ℕ) (hS : ∀ p ∈ S, p.Prime) : oddProd S.card ≤ ∏ p ∈ S, p := by
  induction S using Finset.induction_on_max with
  | empty => simp [oddProd]
  | insert a s hlt ih =>
    have ha : a ∉ s := fun h => lt_irrefl a (hlt a h)
    have hs : ∀ p ∈ s, p.Prime := fun p hp => hS p (Finset.mem_insert_of_mem hp)
    rw [Finset.prod_insert ha, Finset.card_insert_of_notMem ha]
    have hb := max_prime_bound s a (fun p hp => ⟨hs p hp, hlt p hp⟩) (hS a (Finset.mem_insert_self a s))
    show (if s.card = 0 then 2 else 2 * s.card + 1) * oddProd s.card ≤ a * ∏ p ∈ s, p
    exact Nat.mul_le_mul hb (ih hs)

/-- `196(r + 1)²4^r ≤ oddProd r` for `r ≥ 12`. -/
theorem bound_odd (r : ℕ) (hr : 12 ≤ r) : 196 * (r + 1) ^ 2 * 4 ^ r ≤ oddProd r := by
  induction r, hr using Nat.le_induction with
  | base => decide
  | succ r hr ih =>
    have hrr := Nat.mul_le_mul hr hr
    have h1 : 4 * (r + 2) ^ 2 ≤ (2 * r + 1) * (r + 1) ^ 2 := by nlinarith
    have hr0 : r ≠ 0 := by omega
    show 196 * (r + 1 + 1) ^ 2 * 4 ^ (r + 1) ≤ (if r = 0 then 2 else 2 * r + 1) * oddProd r
    rw [ite_eq_right_iff.mpr (fun h => absurd h hr0)]
    calc 196 * (r + 1 + 1) ^ 2 * 4 ^ (r + 1) = 196 * 4 ^ r * (4 * (r + 2) ^ 2) := by ring
      _ ≤ 196 * 4 ^ r * ((2 * r + 1) * (r + 1) ^ 2) := Nat.mul_le_mul_left _ h1
      _ = (2 * r + 1) * (196 * (r + 1) ^ 2 * 4 ^ r) := by ring
      _ ≤ (2 * r + 1) * oddProd r := Nat.mul_le_mul_left _ ih

/-- **Lemma 6 (case (A)), proved part, sharper**: every `q ≥ 118,380,036,096` and every `q ≥ 3`
with `ω(q − 1) ≥ 12` has sieve data with (CS) (no sieve needed). -/
theorem lemma6A_big2 (q : ℕ) (hq : 3 ≤ q)
    (h : 118380036096 ≤ q ∨ 12 ≤ (primeDivisors (q - 1)).card) :
    ∃ (d : SieveData q) (s : ℚ), CSA q d s := by
  refine ⟨noSieve q, _, csa_noSieve q hq ?_⟩
  by_cases hr : 12 ≤ (primeDivisors (q - 1)).card
  · have h1 := bound_odd _ hr
    have h2 := oddProd_le_prod (primeDivisors (q - 1)) (fun p hp => prime_of_mem_primeDivisors hp)
    have h3 := prod_primeDivisors_le (q - 1) (by omega)
    exact le_trans h1 (le_trans h2 (le_trans h3 (Nat.sub_le q 1)))
  · rcases h with h | h
    · push Not at hr
      have h1 : ((primeDivisors (q - 1)).card + 1) ^ 2 ≤ 12 ^ 2 := Nat.pow_le_pow_left (by omega) 2
      have h2 : 4 ^ (primeDivisors (q - 1)).card ≤ 4 ^ 11 := Nat.pow_le_pow_right (by norm_num) (by omega)
      calc 196 * ((primeDivisors (q - 1)).card + 1) ^ 2 * 4 ^ (primeDivisors (q - 1)).card
          ≤ 196 * 12 ^ 2 * 4 ^ 11 := Nat.mul_le_mul (Nat.mul_le_mul_left _ h1) h2
        _ = 118380036096 := by norm_num
        _ ≤ q := h
    · exact absurd h hr

/-- **What is left of Lemma 6** (sharper), NAMED HYPOTHESIS: (CS) for the prime powers
`30140 ≤ q < 118,380,036,096` with `ω(q − 1) ≤ 11` and `φ(q − 1) < 7·2^{ω(q−1)}·(⌊√q⌋ + 1)`. -/
def Lemma6Amid2 : Prop :=
  ∀ q : ℕ, PrimePow q → 30140 ≤ q → q < 118380036096 → (primeDivisors (q - 1)).card ≤ 11 →
    phi (q - 1) < 7 * 2 ^ (primeDivisors (q - 1)).card * (Nat.sqrt q + 1) →
    ∃ (d : SieveData q) (s : ℚ), CSA q d s

/-- Lemma 6 from its sharper middle range. -/
theorem lemma6A_of_mid2 (h : Lemma6Amid2) : Lemma6A := by
  intro q hq h1
  by_cases hbig : 118380036096 ≤ q
  · exact lemma6A_big2 q (by omega) (Or.inl hbig)
  · by_cases hr : 12 ≤ (primeDivisors (q - 1)).card
    · exact lemma6A_big2 q (by omega) (Or.inr hr)
    · by_cases hφ : 7 * 2 ^ (primeDivisors (q - 1)).card * (Nat.sqrt q + 1) ≤ phi (q - 1)
      · exact ⟨noSieve q, _, csa_noSieve_phi q (by omega) hφ⟩
      · exact h q hq h1 (by omega) (by omega) (by omega)

end SunLean

#print axioms SunLean.lemma6A_big2
#print axioms SunLean.lemma6A_of_mid2
