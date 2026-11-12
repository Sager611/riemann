import Riemann.Deep.OddRobin.Lean.LogLike
import Riemann.Deep.OddRobin.Lean.Defs
import Mathlib.Data.Nat.Factorization.Basic
import Mathlib.Tactic.NormNum.Prime

/-!
# The eps-cost machinery (PROOF.md §1 (1a)–(1b)), for an abstract logarithm L

eps = L(68/67)/L 67; G p a = log g_p(a) = L sigma(p^a) − a(1+eps) L p;
inc p a = G p (a+1) − G p a = L(1 + 1/(p sigma(p^a))) − eps L p, antitone in a (concavity).
-/

namespace OddRobin

open Finset ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

/-- eps* = log(68/67)/log 67 (PROOF.md §1). -/
noncomputable def eps (L : ℝ → ℝ) : ℝ := L (68 / 67) / L 67

/-- log g_p(a), g_p(a) = sigma(p^a)/p^(a(1+eps)). -/
noncomputable def G (L : ℝ → ℝ) (p a : ℕ) : ℝ := L (σ 1 (p ^ a) : ℝ) - a * (1 + eps L) * L p

/-- The increment log g_p(a+1) − log g_p(a). -/
noncomputable def inc (L : ℝ → ℝ) (p a : ℕ) : ℝ :=
  L (1 + 1 / ((p : ℝ) * (σ 1 (p ^ a) : ℝ))) - eps L * L p

lemma sigma_pp_pos (p a : ℕ) (hp : p.Prime) : 0 < σ 1 (p ^ a) :=
  sigma_pos_iff.2 (pow_pos hp.pos a)

lemma sigma_succ {p : ℕ} (hp : p.Prime) (a : ℕ) : σ 1 (p ^ (a + 1)) = p * σ 1 (p ^ a) + 1 := by
  rw [sigma_one_apply_prime_pow hp, sigma_one_apply_prime_pow hp, Finset.sum_range_succ',
    Finset.mul_sum, pow_zero]
  congr 1
  exact Finset.sum_congr rfl (fun k _ => by ring)

lemma sigma_mono {p : ℕ} (hp : p.Prime) {a b : ℕ} (h : a ≤ b) : σ 1 (p ^ a) ≤ σ 1 (p ^ b) := by
  rw [sigma_one_apply_prime_pow hp, sigma_one_apply_prime_pow hp]
  exact Finset.sum_le_sum_of_subset (Finset.range_subset_range.2 (by omega))

lemma G_zero (hL : LogLike L) (p : ℕ) : G L p 0 = 0 := by
  simp [G, isMultiplicative_sigma.map_one, hL.map_one]

lemma G_succ (hL : LogLike L) {p : ℕ} (hp : p.Prime) (a : ℕ) :
    G L p (a + 1) = G L p a + inc L p a := by
  have hp0 : (0 : ℝ) < p := by exact_mod_cast hp.pos
  have hs0 : (0 : ℝ) < (σ 1 (p ^ a) : ℝ) := by exact_mod_cast sigma_pp_pos p a hp
  have hps : (0 : ℝ) < (p : ℝ) * (σ 1 (p ^ a) : ℝ) := mul_pos hp0 hs0
  have key : ((σ 1 (p ^ (a + 1)) : ℕ) : ℝ)
      = (p : ℝ) * (σ 1 (p ^ a) : ℝ) * (1 + 1 / ((p : ℝ) * (σ 1 (p ^ a) : ℝ))) := by
    rw [sigma_succ hp a, mul_add, mul_one, mul_one_div_cancel hps.ne']
    push_cast
    ring
  have h1 : (0 : ℝ) < 1 + 1 / ((p : ℝ) * (σ 1 (p ^ a) : ℝ)) := by positivity
  unfold G inc
  rw [key, hL.map_mul _ _ hps h1, hL.map_mul _ _ hp0 hs0]
  push_cast
  ring

lemma inc_anti (hL : LogLike L) {p : ℕ} (hp : p.Prime) {a b : ℕ} (h : a ≤ b) :
    inc L p b ≤ inc L p a := by
  unfold inc
  have hp0 : (0 : ℝ) < p := by exact_mod_cast hp.pos
  have hsa : (0 : ℝ) < (σ 1 (p ^ a) : ℝ) := by exact_mod_cast sigma_pp_pos p a hp
  have hsb : (0 : ℝ) < (σ 1 (p ^ b) : ℝ) := by exact_mod_cast sigma_pp_pos p b hp
  have hab : ((σ 1 (p ^ a) : ℕ) : ℝ) ≤ (σ 1 (p ^ b) : ℝ) := by exact_mod_cast sigma_mono hp h
  have h1 : 1 / ((p : ℝ) * (σ 1 (p ^ b) : ℝ)) ≤ 1 / ((p : ℝ) * (σ 1 (p ^ a) : ℝ)) :=
    one_div_le_one_div_of_le (by positivity) (by gcongr)
  have h2 : (0 : ℝ) < 1 + 1 / ((p : ℝ) * (σ 1 (p ^ b) : ℝ)) := by positivity
  have := hL.mono _ _ h2 (by linarith : 1 + 1 / ((p : ℝ) * (σ 1 (p ^ b) : ℝ))
    ≤ 1 + 1 / ((p : ℝ) * (σ 1 (p ^ a) : ℝ)))
  linarith

lemma G_le_up_aux (hL : LogLike L) {p b : ℕ} (hp : p.Prime) (hb : inc L p b ≤ 0) (d : ℕ) :
    G L p (b + d) ≤ G L p b := by
  induction d with
  | zero => simp
  | succ d ih =>
    rw [← add_assoc, G_succ hL hp]
    have := inc_anti hL hp (Nat.le_add_right b d)
    linarith

/-- Concavity, upper side: if inc p b ≤ 0 then G p a ≤ G p b for a ≥ b. -/
lemma G_le_up (hL : LogLike L) {p b : ℕ} (hp : p.Prime) (hb : inc L p b ≤ 0) {a : ℕ}
    (hab : b ≤ a) : G L p a ≤ G L p b := by
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hab
  exact G_le_up_aux hL hp hb d

lemma G_le_down_aux (hL : LogLike L) {p b : ℕ} (hp : p.Prime) (hb : 0 ≤ inc L p b) (d : ℕ) :
    ∀ a : ℕ, a + d = b + 1 → G L p a ≤ G L p (b + 1) := by
  induction d with
  | zero =>
    intro a h
    have h' : a = b + 1 := by simpa using h
    rw [h']
  | succ d ih =>
    intro a h
    have ha : a ≤ b := by omega
    have h1 : G L p a ≤ G L p (a + 1) := by
      rw [G_succ hL hp]
      have := inc_anti hL hp ha
      linarith
    exact h1.trans (ih (a + 1) (by omega))

/-- Concavity, lower side: if 0 ≤ inc p b then G p a ≤ G p (b + 1) for a ≤ b + 1. -/
lemma G_le_down (hL : LogLike L) {p b : ℕ} (hp : p.Prime) (hb : 0 ≤ inc L p b) {a : ℕ}
    (hab : a ≤ b + 1) : G L p a ≤ G L p (b + 1) :=
  G_le_down_aux hL hp hb (b + 1 - a) a (by omega)

lemma inc_val (p a s : ℕ) (hs : σ 1 (p ^ a) = s) :
    inc L p a = L (1 + 1 / ((p : ℝ) * (s : ℝ))) - eps L * L p := by
  rw [inc, hs]

/-- (1b): the sum of log g_p over the factorization. -/
lemma decomp (hL : LogLike L) {m : ℕ} (hm : m ≠ 0) {S : Finset ℕ} (hS : m.primeFactors ⊆ S) :
    ∑ p ∈ S, G L p (m.factorization p) = L (σ 1 m : ℝ) - (1 + eps L) * L m := by
  have h1 : ∑ p ∈ S, G L p (m.factorization p)
      = ∑ p ∈ m.primeFactors, G L p (m.factorization p) := by
    refine (Finset.sum_subset hS (fun p _ hp => ?_)).symm
    have : m.factorization p = 0 := by
      rw [← Finsupp.notMem_support_iff, Nat.support_factorization]; exact hp
    rw [this, G_zero hL]
  have h2 : (σ 1 m : ℝ) = ∏ p ∈ m.primeFactors, (σ 1 (p ^ m.factorization p) : ℝ) := by
    rw [IsMultiplicative.multiplicative_factorization (σ 1) isMultiplicative_sigma hm, Finsupp.prod,
      Nat.support_factorization, Nat.cast_prod]
  have hpf := Nat.prod_factorization_pow_eq_self hm
  simp only [Finsupp.prod, Nat.support_factorization] at hpf
  have h3 : (m : ℝ) = ∏ p ∈ m.primeFactors, (p : ℝ) ^ m.factorization p := by
    conv_lhs => rw [← hpf]
    push_cast
    rfl
  have h4 : ∀ p ∈ m.primeFactors, (0 : ℝ) < (σ 1 (p ^ m.factorization p) : ℝ) := by
    intro p hp
    exact_mod_cast sigma_pp_pos p _ (Nat.prime_of_mem_primeFactors hp)
  have h5 : ∀ p ∈ m.primeFactors, (0 : ℝ) < (p : ℝ) ^ m.factorization p := by
    intro p hp
    have : (0 : ℝ) < p := by exact_mod_cast (Nat.prime_of_mem_primeFactors hp).pos
    positivity
  rw [h1, h2, hL.map_prod _ _ h4, h3, hL.map_prod _ _ h5, Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun p hp => ?_)
  have : (0 : ℝ) < p := by exact_mod_cast (Nat.prime_of_mem_primeFactors hp).pos
  rw [hL.map_pow this]
  unfold G
  ring

/-! ## N as a product of prime powers -/

/-- The odd primes ≤ 61, i.e. the prime divisors of N. -/
def P61 : Finset ℕ := {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61}

/-- Exponents of N. -/
def αN (p : ℕ) : ℕ := if p = 3 then 4 else if p = 5 then 3 else if p = 7 then 2 else 1

lemma P61_prime : ∀ p ∈ P61, p.Prime := by
  intro p hp
  simp only [P61, Finset.mem_insert, Finset.mem_singleton] at hp
  rcases hp with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
    rfl | rfl <;> norm_num

lemma NN_prod : NN = ∏ p ∈ P61, p ^ αN p := by
  simp [NN, P61, αN]

lemma factorization_prod_pow {s : Finset ℕ} (hs : ∀ p ∈ s, p.Prime) (e : ℕ → ℕ) (q : ℕ) :
    (∏ p ∈ s, p ^ e p).factorization q = if q ∈ s then e q else 0 := by
  rw [Nat.factorization_prod (fun p hp => pow_ne_zero _ (hs p hp).ne_zero), Finsupp.finsetSum_apply,
    Finset.sum_congr rfl (fun p hp => by rw [(hs p hp).factorization_pow])]
  simp only [Finsupp.single_apply]
  exact Finset.sum_ite_eq' s q e

lemma NN_fact (q : ℕ) : NN.factorization q = if q ∈ P61 then αN q else 0 := by
  rw [NN_prod]; exact factorization_prod_pow P61_prime αN q

lemma NN_ne_zero : NN ≠ 0 := by
  norm_num [NN]

lemma NN_pf_sub : NN.primeFactors ⊆ P61 := by
  intro p hp
  have : NN.factorization p ≠ 0 := by
    rw [← Finsupp.mem_support_iff, Nat.support_factorization]; exact hp
  rw [NN_fact] at this
  by_contra h
  simp [h] at this

/-! ## (OR) in log form -/

lemma ORf_iff_log (hL : LogLike L) {E : ℝ} (hE : 0 < E) {n : ℕ} (hn : 0 < n)
    (hll : 0 < L (L n)) :
    ORf L E n ↔ L ((σ 1 n : ℝ) / n) < L E - L 2 + L (L (L n)) := by
  unfold ORf
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have hs : (0 : ℝ) < (σ 1 n : ℝ) := by exact_mod_cast sigma_pos_iff.2 hn
  have hq : (0 : ℝ) < (σ 1 n : ℝ) / n := div_pos hs hn'
  have hr : (0 : ℝ) < E / 2 * L (L n) := by positivity
  rw [← hL.lt_iff hq hr, hL.map_mul _ _ (by positivity) hll, hL.map_div hE two_pos]

end OddRobin

