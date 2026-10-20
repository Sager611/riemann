import Riemann.Deep.Hasse.Count
import Mathlib.Analysis.Real.Sqrt
import Mathlib.Tactic.Linarith

/-!
# Hasse's theorem in odd characteristic (Manin's method, upper-bound form)

For a finite field F of odd characteristic with q elements and any monic cubic
f = X³ + c₂X² + c₁X + c₀, the character sum S = ∑ₓ χ(f(x)) satisfies S² ≤ 4q.
Equivalently the curve y² = f(x) has q + 1 + S projective points, so
|#E(F) − (q + 1)| ≤ 2√q.

Proof. On Manin's twist f(t) y² = f(x) over F(t) let hₙ = ht(P₁ + n P₀) (degree of the
x-coordinate). Then hₙ₊₁ + hₙ₋₁ ≤ 2hₙ + 2 (ht_pgram), h₀ = q (ht_P1) and
h₁ ≤ q + 1 − S, h₋₁ ≤ q + 1 + S (ht_P1_add, ht_P1_sub). Hence 0 ≤ hₙ ≤ n² − Sn + q for all n,
and q odd forces S² ≤ 4q.
-/

open Polynomial

namespace Hasse

open scoped Classical

/-- A sequence with second differences at most 2, pinned at 0 and ±1, lies below the parabola. -/
lemma le_parabola (h : ℤ → ℤ) (S q : ℤ) (hrec : ∀ n, h (n + 1) + h (n - 1) ≤ 2 * h n + 2)
    (h0 : h 0 = q) (h1 : h 1 ≤ q + 1 - S) (hm1 : h (-1) ≤ q + 1 + S) (n : ℤ) :
    h n ≤ n ^ 2 - S * n + q := by
  set e : ℤ → ℤ := fun n ↦ h n - (n ^ 2 - S * n + q) with he
  have hre : ∀ n, e (n + 1) + e (n - 1) ≤ 2 * e n := by
    intro n
    have := hrec n
    have hg : (n + 1) ^ 2 - S * (n + 1) + q + ((n - 1) ^ 2 - S * (n - 1) + q) =
        2 * (n ^ 2 - S * n + q) + 2 := by ring
    simp only [he]
    linarith
  have e0 : e 0 = 0 := by simp [he, h0]
  have e1 : e 1 ≤ 0 := by simp only [he]; linarith
  have em1 : e (-1) ≤ 0 := by simp only [he]; linarith
  have fwd : ∀ k : ℕ, e ((k : ℤ) + 1) ≤ e k ∧ e k ≤ 0 := by
    intro k
    induction k with
    | zero => simp only [Nat.cast_zero, zero_add]; exact ⟨by linarith, by linarith⟩
    | succ k ih =>
      have := hre ((k : ℤ) + 1)
      simp only [add_sub_cancel_right] at this
      push_cast
      constructor <;> linarith [ih.1, ih.2]
  have bwd : ∀ k : ℕ, e (-(k : ℤ) - 1) ≤ e (-(k : ℤ)) ∧ e (-(k : ℤ)) ≤ 0 := by
    intro k
    induction k with
    | zero => simp only [Nat.cast_zero, neg_zero, zero_sub]; exact ⟨by linarith, by linarith⟩
    | succ k ih =>
      have := hre (-(k : ℤ) - 1)
      have e1' : -(k : ℤ) - 1 + 1 = -(k : ℤ) := by ring
      have e2' : -(k : ℤ) - 1 - 1 = -((k : ℤ) + 1) - 1 := by ring
      rw [e1', e2'] at this
      push_cast
      have e3' : -((k : ℤ) + 1) = -(k : ℤ) - 1 := by ring
      rw [e3']
      rw [e3'] at this
      constructor <;> linarith [ih.1, ih.2]
  have key : e n ≤ 0 := by
    rcases le_or_gt 0 n with hn | hn
    · obtain ⟨k, rfl⟩ := Int.eq_ofNat_of_zero_le hn
      exact (fwd k).2
    · obtain ⟨k, hk⟩ := Int.exists_eq_neg_ofNat hn.le
      rw [hk]; exact (bwd k).2
  simp only [he] at key
  linarith

/-- An integer quadratic n² − Sn + q ≥ 0 on all integers, with q odd, has S² ≤ 4q. -/
lemma sq_le_of_parabola_nonneg (S q : ℤ) (hq : q % 2 = 1) (h : ∀ n : ℤ, 0 ≤ n ^ 2 - S * n + q) :
    S ^ 2 ≤ 4 * q := by
  obtain ⟨n, r, hr, hS⟩ : ∃ n r : ℤ, (r = 0 ∨ r = 1) ∧ S = 2 * n + r :=
    ⟨S / 2, S % 2, by omega, by omega⟩
  have hn := h n
  subst hS
  rcases hr with rfl | rfl
  · nlinarith
  · obtain ⟨j, hj⟩ : ∃ j, n * (n + 1) = 2 * j := by
      rcases Int.even_mul_succ_self n with ⟨j, hj⟩; exact ⟨j, by omega⟩
    have hN : n ^ 2 - (2 * n + 1) * n + q = q - n * (n + 1) := by ring
    rw [hN] at hn
    have hsq : (2 * n + 1) ^ 2 = 4 * (n * (n + 1)) + 1 := by ring
    rw [hsq]
    omega

variable {F : Type*} [Field F] [Fintype F]

/-- Hasse's bound for the character sum S = ∑ₓ χ(x³ + c₂x² + c₁x + c₀): S² ≤ 4q. -/
theorem sum_quadraticChar_sq_le (hF : ringChar F ≠ 2) (c₂ c₁ c₀ : F) :
    (∑ a : F, quadraticChar F ((cub c₂ c₁ c₀).eval a)) ^ 2 ≤ 4 * Fintype.card F := by
  set S := ∑ a : F, quadraticChar F ((cub c₂ c₁ c₀).eval a) with hS
  set h2 := Ring.two_ne_zero hF
  set P₀ := P0 c₂ c₁ c₀ h2
  set P₁ := P1 c₂ c₁ c₀ hF
  set hs : ℤ → ℤ := fun n ↦ (ht c₂ c₁ c₀ (P₁ + n • P₀) : ℤ) with hhs
  have hrec : ∀ n, hs (n + 1) + hs (n - 1) ≤ 2 * hs n + 2 := by
    intro n
    have := ht_pgram h2 (P₁ + n • P₀)
    have e1 : P₁ + n • P₀ + P₀ = P₁ + (n + 1) • P₀ := by rw [add_zsmul, one_zsmul]; abel
    have e2 : P₁ + n • P₀ - P₀ = P₁ + (n - 1) • P₀ := by rw [sub_zsmul, one_zsmul]; abel
    rw [e1, e2] at this
    simp only [hhs]
    exact_mod_cast this
  have h0 : hs 0 = Fintype.card F := by
    have e0 : P₁ + (0 : ℤ) • P₀ = P₁ := by rw [zero_zsmul, add_zero]
    simp only [hhs, e0]
    exact_mod_cast ht_P1 hF
  have h1 : hs 1 ≤ Fintype.card F + 1 - S := by
    simp only [hhs, one_smul]
    exact ht_P1_add hF
  have hm1 : hs (-1) ≤ Fintype.card F + 1 + S := by
    simp only [hhs, neg_smul, one_smul, ← sub_eq_add_neg]
    exact ht_P1_sub hF
  have hpar := le_parabola hs S (Fintype.card F) hrec h0 h1 hm1
  have hq : ((Fintype.card F : ℕ) : ℤ) % 2 = 1 := by
    have := FiniteField.odd_card_of_char_ne_two hF; omega
  refine sq_le_of_parabola_nonneg S _ hq (fun n ↦ ?_)
  have := hpar n
  have : (0 : ℤ) ≤ hs n := by simp only [hhs]; positivity
  linarith

/-- Point count: #{(x, y) | y² = f(x)} = q + S. -/
theorem card_affine (hF : ringChar F ≠ 2) (c₂ c₁ c₀ : F) :
    ((Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = (cub c₂ c₁ c₀).eval p.1)).card : ℤ) =
      Fintype.card F + ∑ a : F, quadraticChar F ((cub c₂ c₁ c₀).eval a) := by
  have hfib : ∀ x : F, ((Finset.univ.filter (fun y : F ↦ y ^ 2 = (cub c₂ c₁ c₀).eval x)).card : ℤ) =
      quadraticChar F ((cub c₂ c₁ c₀).eval x) + 1 := by
    intro x
    rw [← quadraticChar_card_sqrts hF]
    congr 2
    ext y
    simp
  have hsplit : (Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = (cub c₂ c₁ c₀).eval p.1)).card =
      ∑ x : F, (Finset.univ.filter (fun y : F ↦ y ^ 2 = (cub c₂ c₁ c₀).eval x)).card := by
    rw [Finset.card_filter, Fintype.sum_prod_type]
    simp only [Finset.card_filter]
  rw [hsplit]
  push_cast
  simp_rw [hfib]
  rw [Finset.sum_add_distrib]
  simp [add_comm]

/-- Hasse's theorem for y² = x³ + c₂x² + c₁x + c₀ over a finite field of odd characteristic
(no nonsingularity assumption is needed): with N = #{affine solutions} + 1,
|N − (q + 1)| ≤ 2√q. -/
theorem hasse (hF : ringChar F ≠ 2) (c₂ c₁ c₀ : F) :
    |((Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 + c₂ * p.1 ^ 2 + c₁ * p.1 + c₀)).card
        + 1 : ℝ) - (Fintype.card F + 1)| ≤ 2 * Real.sqrt (Fintype.card F) := by
  have hc := card_affine hF c₂ c₁ c₀
  have hb := sum_quadraticChar_sq_le hF c₂ c₁ c₀
  have hev : ∀ x : F, (cub c₂ c₁ c₀).eval x = x ^ 3 + c₂ * x ^ 2 + c₁ * x + c₀ := by
    intro x; simp [cub]
  simp_rw [hev] at hc hb
  set S := ∑ a : F, quadraticChar F (a ^ 3 + c₂ * a ^ 2 + c₁ * a + c₀)
  have hcR : ((Finset.univ.filter (fun p : F × F ↦ p.2 ^ 2 = p.1 ^ 3 + c₂ * p.1 ^ 2 + c₁ * p.1 +
      c₀)).card : ℝ) = Fintype.card F + S := by exact_mod_cast hc
  rw [hcR]
  have hbR : (S : ℝ) ^ 2 ≤ 4 * Fintype.card F := by exact_mod_cast hb
  have : ((Fintype.card F : ℝ) + S + 1 - (Fintype.card F + 1)) = S := by ring
  rw [this]
  calc |(S : ℝ)| ≤ Real.sqrt (4 * Fintype.card F) := Real.abs_le_sqrt hbR
    _ = 2 * Real.sqrt (Fintype.card F) := by
      rw [Real.sqrt_mul (by norm_num), show (4 : ℝ) = 2 ^ 2 by norm_num,
        Real.sqrt_sq (by norm_num)]

end Hasse
