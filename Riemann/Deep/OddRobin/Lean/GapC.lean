import Riemann.Deep.OddRobin.Lean.GapP
import Riemann.Deep.OddRobin.Lean.GapF
import Mathlib.Tactic.IntervalCases

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

/-- The exponents allowed for a violator (PROOF.md §1 (1f)). -/
def allowed (p a : ℕ) : Prop :=
  a = NN.factorization p ∨ (p = 3 ∧ a = 5) ∨ (p = 5 ∧ a = 2) ∨ (p = 11 ∧ a = 2) ∨
    (p = 67 ∧ a = 1) ∨ (p = 71 ∧ a = 1)

/-- The other prime divisors of N. -/
def Pfix : Finset ℕ := {7, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61}

/-- The odd primes below 73. -/
def P71 : Finset ℕ := insert 3 (insert 5 (insert 11 (insert 67 (insert 71 Pfix))))

lemma Pfix_sub_P61 : Pfix ⊆ P61 := by decide

lemma P61_sub_P71 : P61 ⊆ P71 := by decide

lemma Pfix_prod : ∏ p ∈ Pfix, p ^ αN p = KK := by
  simp [Pfix, αN, KK]

lemma mid_case (hL : LogLike L) {p : ℕ} (hp : p.Prime) (h13 : 13 ≤ p) (h61 : p ≤ 61)
    (hpP : p ∈ P61) (a : ℕ) :
    0 ≤ G L p (NN.factorization p) - G L p a ∧
      (G L p (NN.factorization p) - G L p a < 13 / 10000 → allowed p a) := by
  have hf : NN.factorization p = 1 := by
    rw [NN_fact, if_pos hpP]
    unfold αN
    split_ifs <;> omega
  rw [hf]
  obtain ⟨h1, h2⟩ := pmid hL hp h13 h61 a
  exact ⟨h1, fun hc => Or.inl ((h2 hc).trans hf.symm)⟩

lemma prime_lt_73 {p : ℕ} (hp : p.Prime) (h2 : p ≠ 2) (h : p < 73) : p ∈ P71 := by
  interval_cases p <;> first | (norm_num at hp; done) | (exact absurd rfl h2) | decide

lemma cost_prop (hL : LogLike L) {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) (a : ℕ) :
    0 ≤ G L p (NN.factorization p) - G L p a ∧
      (G L p (NN.factorization p) - G L p a < 13 / 10000 → allowed p a) := by
  rcases lt_or_ge p 73 with h | h
  · have hmem := prime_lt_73 hp hp2 h
    simp only [P71, Pfix, Finset.mem_insert, Finset.mem_singleton] at hmem
    rcases hmem with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
      rfl | rfl | rfl | rfl | rfl | rfl
    · have hf : NN.factorization 3 = 4 := by rw [NN_fact]; simp [P61, αN]
      rw [hf]
      obtain ⟨h1, h2⟩ := p3 hL a
      refine ⟨h1, fun hc => ?_⟩
      rcases h2 hc with h' | h'
      · exact Or.inl (h'.trans hf.symm)
      · exact Or.inr (Or.inl ⟨rfl, h'⟩)
    · have hf : NN.factorization 5 = 3 := by rw [NN_fact]; simp [P61, αN]
      rw [hf]
      obtain ⟨h1, h2⟩ := p5 hL a
      refine ⟨h1, fun hc => ?_⟩
      rcases h2 hc with h' | h'
      · exact Or.inr (Or.inr (Or.inl ⟨rfl, h'⟩))
      · exact Or.inl (h'.trans hf.symm)
    · have hf : NN.factorization 11 = 1 := by rw [NN_fact]; simp [P61, αN]
      rw [hf]
      obtain ⟨h1, h2⟩ := p11 hL a
      refine ⟨h1, fun hc => ?_⟩
      rcases h2 hc with h' | h'
      · exact Or.inl (h'.trans hf.symm)
      · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨rfl, h'⟩)))
    · have hf : NN.factorization 67 = 0 := by rw [NN_fact]; simp [P61]
      rw [hf]
      obtain ⟨h1, h2⟩ := p67 hL a
      refine ⟨h1, fun hc => ?_⟩
      rcases h2 hc with h' | h'
      · exact Or.inl (h'.trans hf.symm)
      · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inl ⟨rfl, h'⟩))))
    · have hf : NN.factorization 71 = 0 := by rw [NN_fact]; simp [P61]
      rw [hf]
      obtain ⟨h1, h2⟩ := p71 hL a
      refine ⟨h1, fun hc => ?_⟩
      rcases h2 hc with h' | h'
      · exact Or.inl (h'.trans hf.symm)
      · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr ⟨rfl, h'⟩))))
    · have hf : NN.factorization 7 = 2 := by rw [NN_fact]; simp [P61, αN]
      rw [hf]
      obtain ⟨h1, h2⟩ := p7 hL a
      exact ⟨h1, fun hc => Or.inl ((h2 hc).trans hf.symm)⟩
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
    · exact mid_case hL hp (by norm_num) (by norm_num) (by simp [P61]) a
  · have hf : NN.factorization p = 0 := by
      rw [NN_fact, if_neg]
      intro hm
      simp only [P61, Finset.mem_insert, Finset.mem_singleton] at hm
      omega
    rw [hf]
    obtain ⟨h1, h2⟩ := pbig hL hp h a
    exact ⟨h1, fun hc => Or.inl ((h2 hc).trans hf.symm)⟩


end OddRobin
