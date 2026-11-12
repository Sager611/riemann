import Riemann.Deep.OddRobin.Lean.GapD
import Riemann.Deep.OddRobin.Lean.GapC

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

/-- THE GAP THEOREM (PROOF.md §1): every odd n with N < n < c0 satisfies (OR). -/
theorem gap_thm (hL : LogLike L) {E g : ℝ} (hE : 0 < E) (hLE : L E = g)
    (hg : (5772145 / 10000000 : ℝ) ≤ g) {n : ℕ} (hodd : Odd n) (h1 : NN < n) (h2 : n < c0) :
    ORf L E n := by
  by_contra hv
  have hn0 : n ≠ 0 := by
    rintro rfl
    exact absurd h1 (Nat.not_lt_zero _)
  have hsum := total_cost hL hE hLE hg hn0 h1 h2 hv
  have hprime : ∀ p ∈ n.primeFactors ∪ P61, p.Prime ∧ p ≠ 2 := by
    intro p hp
    rcases Finset.mem_union.1 hp with hp | hp
    · refine ⟨Nat.prime_of_mem_primeFactors hp, ?_⟩
      rintro rfl
      have h2d := Nat.dvd_of_mem_primeFactors hp
      exact (Nat.not_even_iff_odd.2 hodd) (even_iff_two_dvd.2 h2d)
    · refine ⟨P61_prime p hp, ?_⟩
      rintro rfl
      simp [P61] at hp
  have hnn : ∀ p ∈ n.primeFactors ∪ P61,
      0 ≤ G L p (NN.factorization p) - G L p (n.factorization p) :=
    fun p hp => (cost_prop hL (hprime p hp).1 (hprime p hp).2 _).1
  have hsingle : ∀ p ∈ n.primeFactors ∪ P61,
      G L p (NN.factorization p) - G L p (n.factorization p) < 13 / 10000 :=
    fun p hp => lt_of_le_of_lt (Finset.single_le_sum hnn hp) hsum
  have hpair : ∀ p ∈ n.primeFactors ∪ P61, ∀ q ∈ n.primeFactors ∪ P61, p ≠ q →
      (G L p (NN.factorization p) - G L p (n.factorization p)) +
        (G L q (NN.factorization q) - G L q (n.factorization q)) < 13 / 10000 := by
    intro p hp q hq hpq
    have hsub : ({p, q} : Finset ℕ) ⊆ n.primeFactors ∪ P61 := by
      intro x hx
      simp only [Finset.mem_insert, Finset.mem_singleton] at hx
      rcases hx with rfl | rfl <;> assumption
    have := Finset.sum_le_sum_of_subset_of_nonneg hsub (fun x hx _ => hnn x hx)
    rw [Finset.sum_pair hpq] at this
    linarith
  have hall : ∀ p, allowed p (n.factorization p) := by
    intro p
    by_cases hpS : p ∈ n.primeFactors ∪ P61
    · exact (cost_prop hL (hprime p hpS).1 (hprime p hpS).2 _).2 (hsingle p hpS)
    · have h1' : p ∉ n.primeFactors := fun h => hpS (Finset.mem_union_left _ h)
      have h2' : p ∉ P61 := fun h => hpS (Finset.mem_union_right _ h)
      have e1 : n.factorization p = 0 := by
        rw [← Finsupp.notMem_support_iff, Nat.support_factorization]; exact h1'
      have e2 : NN.factorization p = 0 := by rw [NN_fact, if_neg h2']
      exact Or.inl (e1.trans e2.symm)
  -- the exponents
  have nf3 : NN.factorization 3 = 4 := by rw [NN_fact]; simp [P61, αN]
  have nf5 : NN.factorization 5 = 3 := by rw [NN_fact]; simp [P61, αN]
  have nf11 : NN.factorization 11 = 1 := by rw [NN_fact]; simp [P61, αN]
  have nf67 : NN.factorization 67 = 0 := by rw [NN_fact]; simp [P61]
  have nf71 : NN.factorization 71 = 0 := by rw [NN_fact]; simp [P61]
  have e3 : n.factorization 3 = 4 ∨ n.factorization 3 = 5 := by
    rcases hall 3 with h | h | h | h | h | h
    · exact Or.inl (h.trans nf3)
    · exact Or.inr h.2
    all_goals (exfalso; omega)
  have e5 : n.factorization 5 = 3 ∨ n.factorization 5 = 2 := by
    rcases hall 5 with h | h | h | h | h | h
    · exact Or.inl (h.trans nf5)
    · exfalso; omega
    · exact Or.inr h.2
    all_goals (exfalso; omega)
  have e11 : n.factorization 11 = 1 ∨ n.factorization 11 = 2 := by
    rcases hall 11 with h | h | h | h | h | h
    · exact Or.inl (h.trans nf11)
    · exfalso; omega
    · exfalso; omega
    · exact Or.inr h.2
    all_goals (exfalso; omega)
  have e67 : n.factorization 67 = 0 ∨ n.factorization 67 = 1 := by
    rcases hall 67 with h | h | h | h | h | h
    · exact Or.inl (h.trans nf67)
    · exfalso; omega
    · exfalso; omega
    · exfalso; omega
    · exact Or.inr h.2
    · exfalso; omega
  have e71 : n.factorization 71 = 0 ∨ n.factorization 71 = 1 := by
    rcases hall 71 with h | h | h | h | h | h
    · exact Or.inl (h.trans nf71)
    · exfalso; omega
    · exfalso; omega
    · exfalso; omega
    · exfalso; omega
    · exact Or.inr h.2
  have hfix : ∀ p ∈ Pfix, n.factorization p = αN p := by
    intro p hp
    have hp61 : p ∈ P61 := Pfix_sub_P61 hp
    rcases hall p with h | h | h | h | h | h
    · rw [h, NN_fact, if_pos hp61]
    all_goals (obtain ⟨rfl, _⟩ := h; simp [Pfix] at hp)
  -- pair exclusions (any two deviations cost more than 13/10000)
  have m3S : 3 ∈ n.primeFactors ∪ P61 := Finset.mem_union_right _ (by simp [P61])
  have m5S : 5 ∈ n.primeFactors ∪ P61 := Finset.mem_union_right _ (by simp [P61])
  have m11S : 11 ∈ n.primeFactors ∪ P61 := Finset.mem_union_right _ (by simp [P61])
  have m71S : n.factorization 71 = 1 → 71 ∈ n.primeFactors ∪ P61 := by
    intro h
    refine Finset.mem_union_left _ ?_
    rw [← Nat.support_factorization, Finsupp.mem_support_iff, h]
    norm_num
  have c3 : n.factorization 3 = 5 →
      11 / 10000 < G L 3 (NN.factorization 3) - G L 3 (n.factorization 3) := by
    intro h; rw [h, nf3]; exact low3 hL
  have c5 : n.factorization 5 = 2 →
      75 / 100000 < G L 5 (NN.factorization 5) - G L 5 (n.factorization 5) := by
    intro h; rw [h, nf5]; exact low5 hL
  have c11 : n.factorization 11 = 2 →
      89 / 100000 < G L 11 (NN.factorization 11) - G L 11 (n.factorization 11) := by
    intro h; rw [h, nf11]; exact low11 hL
  have c71 : n.factorization 71 = 1 →
      1 / 1000 < G L 71 (NN.factorization 71) - G L 71 (n.factorization 71) := by
    intro h; rw [h, nf71]; exact low71 hL
  have x35 : ¬ (n.factorization 3 = 5 ∧ n.factorization 5 = 2) := fun ⟨ha, hb⟩ => by
    have := hpair 3 m3S 5 m5S (by norm_num); linarith [c3 ha, c5 hb]
  have x311 : ¬ (n.factorization 3 = 5 ∧ n.factorization 11 = 2) := fun ⟨ha, hb⟩ => by
    have := hpair 3 m3S 11 m11S (by norm_num); linarith [c3 ha, c11 hb]
  have x371 : ¬ (n.factorization 3 = 5 ∧ n.factorization 71 = 1) := fun ⟨ha, hb⟩ => by
    have := hpair 3 m3S 71 (m71S hb) (by norm_num); linarith [c3 ha, c71 hb]
  have x511 : ¬ (n.factorization 5 = 2 ∧ n.factorization 11 = 2) := fun ⟨ha, hb⟩ => by
    have := hpair 5 m5S 11 m11S (by norm_num); linarith [c5 ha, c11 hb]
  have x571 : ¬ (n.factorization 5 = 2 ∧ n.factorization 71 = 1) := fun ⟨ha, hb⟩ => by
    have := hpair 5 m5S 71 (m71S hb) (by norm_num); linarith [c5 ha, c71 hb]
  have x1171 : ¬ (n.factorization 11 = 2 ∧ n.factorization 71 = 1) := fun ⟨ha, hb⟩ => by
    have := hpair 11 m11S 71 (m71S hb) (by norm_num); linarith [c11 ha, c71 hb]
  -- reconstruct n
  have hpf : n.primeFactors ⊆ P71 := by
    intro p hp
    have hne : n.factorization p ≠ 0 := by
      rw [← Finsupp.mem_support_iff, Nat.support_factorization]; exact hp
    rcases hall p with h | h | h | h | h | h
    · rw [h, NN_fact] at hne
      by_cases hm : p ∈ P61
      · exact P61_sub_P71 hm
      · simp [hm] at hne
    all_goals (obtain ⟨rfl, _⟩ := h; decide)
  have hn_eq : n = ∏ p ∈ P71, p ^ n.factorization p := by
    have hpe := Nat.prod_factorization_pow_eq_self hn0
    simp only [Finsupp.prod, Nat.support_factorization] at hpe
    have hsub : ∀ p ∈ P71, p ∉ n.primeFactors → p ^ n.factorization p = 1 := by
      intro p _ hp
      have : n.factorization p = 0 := by
        rw [← Finsupp.notMem_support_iff, Nat.support_factorization]; exact hp
      rw [this, pow_zero]
    rw [Finset.prod_subset hpf hsub] at hpe
    exact hpe.symm
  have hfixprod : ∏ p ∈ Pfix, p ^ n.factorization p = KK := by
    rw [Finset.prod_congr rfl (fun p hp => by rw [hfix p hp]), Pfix_prod]
  have hn_eq2 : n = 3 ^ n.factorization 3 * (5 ^ n.factorization 5 * (11 ^ n.factorization 11 *
      (67 ^ n.factorization 67 * (71 ^ n.factorization 71 * KK)))) := by
    refine hn_eq.trans ?_
    rw [P71, Finset.prod_insert (by decide), Finset.prod_insert (by decide),
      Finset.prod_insert (by decide), Finset.prod_insert (by decide),
      Finset.prod_insert (by decide), hfixprod]
  have hcases := final_cases _ _ _ _ _ e3 e5 e11 e67 e71 x35 x311 x371 x511 x571 x1171
    (lt_of_lt_of_eq h1 hn_eq2) (lt_of_eq_of_lt hn_eq2.symm h2)
  rw [← hn_eq2] at hcases
  rcases hcases with h | h | h
  · exact hv (h ▸ OR_m3 hL hE hLE hg)
  · exact hv (h ▸ OR_m11 hL hE hLE hg)
  · exact hv (h ▸ OR_m67 hL hE hLE hg)


end OddRobin
