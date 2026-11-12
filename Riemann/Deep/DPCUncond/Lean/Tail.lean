import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.Order.Ring.Pow
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Positivity

/-!
# dpcuncond Step 6f, part 1 (new step (d)): abstract exp/log and the polynomial-times-geometric tail

ExpLog ex lg: an abstract exponential and logarithm with exactly the four Mathlib facts used (Real.exp_add,
Real.exp_pos, Real.exp_strictMono, Real.exp_log); Mathlib's Real.log is not importable under the 290 MB cap
(OddRobin/LEAN.md: Log.Basic alone 326 MB). polygeom: (a + 2)^C t^a ≤ K s^a with s < 1 (Bernoulli twice).
tail_small: Σ_{A<a≤M} (a + 2)^C q ρ^a ≤ τ for all A ≥ N and all M, for abstract q ≥ 0, 0 ≤ ρ < 1 (the
summability fact of 6f, proved, uniform in the truncation M).
-/

namespace DPCLean

open Finset

/-- Abstract exp and log: the four facts of Mathlib's Real.exp, Real.log used here. -/
structure ExpLog (ex lg : ℝ → ℝ) : Prop where
  add : ∀ x y : ℝ, ex (x + y) = ex x * ex y
  pos : ∀ x : ℝ, 0 < ex x
  mono : StrictMono ex
  ex_lg : ∀ x : ℝ, 0 < x → ex (lg x) = x

variable {ex lg : ℝ → ℝ}

lemma ExpLog.zero (h : ExpLog ex lg) : ex 0 = 1 := by
  have h1 := h.add 0 0
  rw [add_zero] at h1
  exact mul_left_cancel₀ (h.pos 0).ne' (by rw [mul_one]; exact h1.symm)

lemma ExpLog.nat_mul (h : ExpLog ex lg) (n : ℕ) (x : ℝ) : ex (n * x) = ex x ^ n := by
  induction n with
  | zero => rw [Nat.cast_zero, zero_mul, h.zero, pow_zero]
  | succ n ih => rw [Nat.cast_succ, add_mul, one_mul, h.add, ih, pow_succ]

lemma ExpLog.lg_mul (h : ExpLog ex lg) {x y : ℝ} (hx : 0 < x) (hy : 0 < y) :
    lg (x * y) = lg x + lg y :=
  h.mono.injective (by rw [h.add, h.ex_lg _ hx, h.ex_lg _ hy, h.ex_lg _ (mul_pos hx hy)])

lemma ExpLog.lg_one (h : ExpLog ex lg) : lg 1 = 0 :=
  h.mono.injective (by rw [h.ex_lg _ one_pos, h.zero])

lemma ExpLog.lg_inv (h : ExpLog ex lg) {x : ℝ} (hx : 0 < x) : lg x = -lg (1 / x) := by
  have h1 := h.lg_mul hx (one_div_pos.2 hx)
  have e : x * (1 / x) = 1 := by rw [mul_one_div, div_self hx.ne']
  rw [e, h.lg_one] at h1
  linarith

lemma ExpLog.le_lg (h : ExpLog ex lg) {M x : ℝ} (hx : ex M ≤ x) : M ≤ lg x := by
  rcases le_or_gt M (lg x) with h1 | h1
  · exact h1
  · have hx0 : 0 < x := lt_of_lt_of_le (h.pos M) hx
    have h2 := h.mono h1
    rw [h.ex_lg _ hx0] at h2
    linarith

lemma ExpLog.lt_one (h : ExpLog ex lg) {x : ℝ} (hx : x < 0) : ex x < 1 := by
  have h1 := h.mono hx
  rwa [h.zero] at h1

/-- Bernoulli: s, u ≥ 0 and s(1 + u) = 1 give s^a (1 + a u) ≤ 1. -/
lemma pow_mul_le_one {s u : ℝ} (hs : 0 ≤ s) (hu : 0 ≤ u) (hsu : s * (1 + u) = 1) (a : ℕ) :
    s ^ a * (1 + a * u) ≤ 1 := by
  have hb : 1 + a * u ≤ (1 + u) ^ a := one_add_mul_le_pow (by linarith) a
  calc s ^ a * (1 + a * u) ≤ s ^ a * (1 + u) ^ a := mul_le_mul_of_nonneg_left hb (pow_nonneg hs a)
    _ = (s * (1 + u)) ^ a := (mul_pow s (1 + u) a).symm
    _ = 1 := by rw [hsu, one_pow]

/-- Polynomial times geometric: for 0 ≤ t < 1 there are K ≥ 0, 0 < s < 1, u > 0 with s(1 + u) = 1 and
(a + 2)^C t^a ≤ K s^a for every a ∈ ℕ (s = 1 − (1 − t)/(C + 2), t ≤ s^{C+2} by Bernoulli). -/
lemma polygeom (C : ℕ) {t : ℝ} (ht0 : 0 ≤ t) (ht1 : t < 1) :
    ∃ K s u : ℝ, 0 ≤ K ∧ 0 < s ∧ s < 1 ∧ 0 < u ∧ s * (1 + u) = 1 ∧
      ∀ a : ℕ, ((a : ℝ) + 2) ^ C * t ^ a ≤ K * s ^ a := by
  have hC : (0 : ℝ) < (C : ℝ) + 2 := by positivity
  have hCne : (C : ℝ) + 2 ≠ 0 := hC.ne'
  have hd0 : 0 < (1 - t) / ((C : ℝ) + 2) := div_pos (by linarith) hC
  have hd1 : (1 - t) / ((C : ℝ) + 2) ≤ 1 / 2 := by
    rw [div_le_iff₀ hC]
    have hC0 : (0 : ℝ) ≤ C := Nat.cast_nonneg C
    nlinarith
  obtain ⟨s, hs⟩ : ∃ s : ℝ, s = 1 - (1 - t) / ((C : ℝ) + 2) := ⟨_, rfl⟩
  have hs0 : 0 < s := by rw [hs]; linarith
  have hs1 : s < 1 := by rw [hs]; linarith
  obtain ⟨u, hu⟩ : ∃ u : ℝ, u = (1 - s) / s := ⟨_, rfl⟩
  have hu0 : 0 < u := by rw [hu]; exact div_pos (by linarith) hs0
  have hsu : s * (1 + u) = 1 := by
    calc s * (1 + u) = s + (1 - s) / s * s := by rw [hu]; ring
      _ = s + (1 - s) := by rw [div_mul_cancel₀ _ hs0.ne']
      _ = 1 := by ring
  have hts : t ≤ s ^ (C + 2) := by
    have hb := one_add_mul_le_pow (show (-2 : ℝ) ≤ -((1 - t) / ((C : ℝ) + 2)) by linarith) (C + 2)
    have e2 : ((C + 2 : ℕ) : ℝ) = (C : ℝ) + 2 := by norm_num
    have e1 : (1 : ℝ) + ((C + 2 : ℕ) : ℝ) * -((1 - t) / ((C : ℝ) + 2)) = t := by
      rw [e2, mul_neg, mul_comm, div_mul_cancel₀ _ hCne]
      ring
    have e3 : (1 : ℝ) + -((1 - t) / ((C : ℝ) + 2)) = s := by rw [hs]; ring
    rw [e1, e3] at hb
    exact hb
  obtain ⟨K₀, hK₀⟩ : ∃ K₀ : ℝ, K₀ = max 2 (1 / u) := ⟨_, rfl⟩
  have hK₀2 : 2 ≤ K₀ := by rw [hK₀]; exact le_max_left _ _
  have hK₀u : 1 ≤ K₀ * u := by
    have h1 : 1 / u ≤ K₀ := by rw [hK₀]; exact le_max_right _ _
    rwa [div_le_iff₀ hu0] at h1
  refine ⟨K₀ ^ C, s, u, pow_nonneg (by linarith) C, hs0, hs1, hu0, hsu, fun a => ?_⟩
  have hsa0 : 0 ≤ s ^ a := pow_nonneg hs0.le a
  have hv : s ^ a * (1 + a * u) ≤ 1 := pow_mul_le_one hs0.le hu0.le hsu a
  have ha0 : (0 : ℝ) ≤ a := Nat.cast_nonneg a
  have hlin : (a : ℝ) + 2 ≤ K₀ * (1 + a * u) := by
    nlinarith [mul_nonneg ha0 (sub_nonneg.2 hK₀u)]
  have hA : ((a : ℝ) + 2) * s ^ a ≤ K₀ := by
    calc ((a : ℝ) + 2) * s ^ a ≤ K₀ * (1 + a * u) * s ^ a := mul_le_mul_of_nonneg_right hlin hsa0
      _ = K₀ * (s ^ a * (1 + a * u)) := by ring
      _ ≤ K₀ * 1 := mul_le_mul_of_nonneg_left hv (by linarith)
      _ = K₀ := mul_one _
  have h1 : t ≤ s ^ (C + 1) := by
    calc t ≤ s ^ (C + 2) := hts
      _ = s ^ (C + 1) * s := pow_succ s (C + 1)
      _ ≤ s ^ (C + 1) * 1 := mul_le_mul_of_nonneg_left hs1.le (pow_nonneg hs0.le _)
      _ = s ^ (C + 1) := mul_one _
  have hta : t ^ a ≤ (s ^ a) ^ C * s ^ a := by
    calc t ^ a ≤ (s ^ (C + 1)) ^ a := pow_le_pow_left₀ ht0 h1 a
      _ = (s ^ a) ^ C * s ^ a := by rw [pow_right_comm s (C + 1) a, pow_succ]
  calc ((a : ℝ) + 2) ^ C * t ^ a ≤ ((a : ℝ) + 2) ^ C * ((s ^ a) ^ C * s ^ a) :=
        mul_le_mul_of_nonneg_left hta (by positivity)
    _ = (((a : ℝ) + 2) * s ^ a) ^ C * s ^ a := by rw [mul_pow]; ring
    _ ≤ K₀ ^ C * s ^ a := mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (by positivity) hA C) hsa0

/-- (1 − s) Σ_{A<a≤M} s^a = s^{A+1} − s^{M+1}. -/
lemma geom_Ioc (s : ℝ) {A M : ℕ} (hAM : A ≤ M) :
    (1 - s) * ∑ a ∈ Ioc A M, s ^ a = s ^ (A + 1) - s ^ (M + 1) := by
  induction M, hAM using Nat.le_induction with
  | base => simp
  | succ M hAM ih => rw [sum_Ioc_succ_top hAM, mul_add, ih]; ring

lemma geom_tail_le {s : ℝ} (hs0 : 0 ≤ s) (hs1 : s < 1) (A M : ℕ) :
    ∑ a ∈ Ioc A M, s ^ a ≤ s ^ (A + 1) / (1 - s) := by
  rw [le_div_iff₀ (by linarith)]
  rcases le_or_gt A M with h | h
  · have h1 := geom_Ioc s h
    have h2 : 0 ≤ s ^ (M + 1) := pow_nonneg hs0 _
    nlinarith
  · rw [Ioc_eq_empty_of_le h.le, sum_empty, zero_mul]
    exact pow_nonneg hs0 _

/-- **The tail of 6f**: for q ≥ 0 and 0 ≤ ρ < 1, Σ_{A<a≤M} (a + 2)^C q ρ^a ≤ τ for all A ≥ N, uniformly in M. -/
theorem tail_small (C : ℕ) {q ρ : ℝ} (hq : 0 ≤ q) (hρ0 : 0 ≤ ρ) (hρ1 : ρ < 1) :
    ∀ τ : ℝ, 0 < τ → ∃ N : ℕ, ∀ A M : ℕ, N ≤ A →
      ∑ a ∈ Ioc A M, ((a : ℝ) + 2) ^ C * q * ρ ^ a ≤ τ := by
  intro τ hτ
  obtain ⟨K, s, u, hK, hs0, hs1, hu, hsu, hb⟩ := polygeom C hρ0 hρ1
  obtain ⟨N, hN⟩ := exists_nat_gt (q * K / ((1 - s) * u * τ))
  refine ⟨N, fun A M hNA => ?_⟩
  have h1s : 0 < 1 - s := by linarith
  have hden : 0 < (1 - s) * u * τ := by positivity
  have h1 : ∑ a ∈ Ioc A M, ((a : ℝ) + 2) ^ C * q * ρ ^ a ≤ ∑ a ∈ Ioc A M, q * K * s ^ a := by
    apply sum_le_sum
    intro a _
    calc ((a : ℝ) + 2) ^ C * q * ρ ^ a = q * (((a : ℝ) + 2) ^ C * ρ ^ a) := by ring
      _ ≤ q * (K * s ^ a) := mul_le_mul_of_nonneg_left (hb a) hq
      _ = q * K * s ^ a := by ring
  have h2 : ∑ a ∈ Ioc A M, q * K * s ^ a = q * K * ∑ a ∈ Ioc A M, s ^ a := by rw [mul_sum]
  have h3 := geom_tail_le hs0.le hs1 A M
  have h4 := pow_mul_le_one hs0.le hu.le hsu (A + 1)
  have hX0 : 0 ≤ s ^ (A + 1) := pow_nonneg hs0.le _
  have hA1 : (N : ℝ) ≤ ((A + 1 : ℕ) : ℝ) := by exact_mod_cast (by omega : N ≤ A + 1)
  have hqK : q * K < ((A + 1 : ℕ) : ℝ) * ((1 - s) * u * τ) := by
    rw [div_lt_iff₀ hden] at hN
    exact lt_of_lt_of_le hN (mul_le_mul_of_nonneg_right hA1 hden.le)
  have hXu : s ^ (A + 1) * (((A + 1 : ℕ) : ℝ) * u) ≤ 1 := by nlinarith
  have h5 : q * K * s ^ (A + 1) ≤ (1 - s) * τ := by
    calc q * K * s ^ (A + 1) ≤ (((A + 1 : ℕ) : ℝ) * ((1 - s) * u * τ)) * s ^ (A + 1) :=
          mul_le_mul_of_nonneg_right hqK.le hX0
      _ = ((1 - s) * τ) * (s ^ (A + 1) * (((A + 1 : ℕ) : ℝ) * u)) := by ring
      _ ≤ ((1 - s) * τ) * 1 := mul_le_mul_of_nonneg_left hXu (by positivity)
      _ = (1 - s) * τ := mul_one _
  have h6 : q * K * ∑ a ∈ Ioc A M, s ^ a ≤ q * K * (s ^ (A + 1) / (1 - s)) :=
    mul_le_mul_of_nonneg_left h3 (mul_nonneg hq hK)
  have h7 : q * K * (s ^ (A + 1) / (1 - s)) ≤ τ := by
    rw [mul_div_assoc', div_le_iff₀ h1s]
    linarith
  calc ∑ a ∈ Ioc A M, ((a : ℝ) + 2) ^ C * q * ρ ^ a ≤ ∑ a ∈ Ioc A M, q * K * s ^ a := h1
    _ = q * K * ∑ a ∈ Ioc A M, s ^ a := h2
    _ ≤ q * K * (s ^ (A + 1) / (1 - s)) := h6
    _ ≤ τ := h7

end DPCLean

