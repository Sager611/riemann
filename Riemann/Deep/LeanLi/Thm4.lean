module

public import Riemann.Deep.LeanLi.Thm4b

/-!
# LeanLi: Theorem 4 (RH): `Σ_n φ(n/N) R(n) = (7/4) Σ_n φ(n/N) + O(log N)` for `φ ∈ C_c^∞(0, ∞)`

Lemma A at a height `T₀` below the first zero and `X = N⁸`; the smoothed oscillating terms are
`O(1)` by the decay of `Σ φ(n/N) n^j e^{inθ}`; the zero term `Σ φ(n/N) 2n a(n) = ∫ S G` is `O(log N)` after
one integration by parts against `S₁ = O(log)`; the constant is `2 + 2∫_0^∞ r = 2 - 1/4 = 7/4`.
-/

@[expose] public section

open Real Filter Set Finset MeasureTheory Topology

namespace LeanLi

/-- Exact form of `R(n) - 7/4` at a height `T₀` below the first zero. -/
lemma R_sub_eq (L : LitCore) (hRH : RiemannHypothesis) {T₀ X : ℝ} (hT₀ : 0 < T₀) (hTX : T₀ ≤ X)
    (hz : zset T₀ = ∅) (n : ℕ) :
    R n - 7 / 4 = 2 * n * aCoef T₀ X n + (-2 * Afun T₀ * Real.cos (n * liTheta T₀)) +
      2 * ((∫ t in (0 : ℝ)..X, rfun t) + 1 / 8) + (-2 * ∫ t in T₀..X, Real.cos (n * liTheta t) * rfun t) +
      (Iclosed n - liTrend n) + (-∫ t in Ioi X, hI n t) +
      1 / Real.pi * (∫ t in (0 : ℝ)..T₀, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))) +
      gfun n X * ((zeroCount X : ℝ) - Afun X) + ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k) := by
  have hX : 0 < X := lt_of_lt_of_le hT₀ hTX
  have hd := decomp L hRH n hT₀ hTX
  have hW : Wv T₀ n = 0 := by unfold Wv; rw [ZK_eq_empty hz, Finset.sum_empty]
  have hN0 := zeroCount_eq_zero hz
  have hc : Continuous fun t => Real.cos (n * liTheta t) := by have := continuous_liTheta; fun_prop
  have hr1 := intervalIntegrable_rfun hT₀.le hX.le
  have hr2 : IntervalIntegrable (fun t => Real.cos (n * liTheta t) * rfun t) volume T₀ X :=
    hr1.continuousOn_mul hc.continuousOn
  have hg : ∫ t in T₀..X, gfun n t * rfun t =
      2 * (∫ t in T₀..X, rfun t) - 2 * ∫ t in T₀..X, Real.cos (n * liTheta t) * rfun t := by
    rw [← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_const_mul,
      ← intervalIntegral.integral_sub (hr1.const_mul 2) (hr2.const_mul 2)]
    refine intervalIntegral.integral_congr fun t _ => ?_
    simp only [gfun]; ring
  have hadd := intervalIntegral.integral_add_adjacent_intervals (intervalIntegrable_rfun le_rfl hT₀.le) hr1
  rw [hd, hW]
  unfold Rest
  rw [hN0, hg]
  linarith

lemma abs_sum_le_card_mul {s : Finset ℕ} {f : ℕ → ℝ} {M : ℝ} (h : ∀ n ∈ s, |f n| ≤ M) :
    |∑ n ∈ s, f n| ≤ s.card * M := by
  refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
  have := Finset.sum_le_card_nsmul s (fun n => |f n|) M h
  rwa [nsmul_eq_mul] at this

section
variable (φ : ℝ → ℝ) (B : ℕ) (N : ℝ)

lemma sum_mul_integral_cos_r {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) :
    ∑ n ∈ range (⌊(B : ℝ) * N⌋₊ + 1), φ (n / N) * (∫ t in a..b, Real.cos (n * liTheta t) * rfun t) =
      ∫ t in a..b, F0 φ B N t * rfun t := by
  have hi : ∀ n ∈ range (⌊(B : ℝ) * N⌋₊ + 1), IntervalIntegrable
      (fun t => φ (n / N) * (Real.cos (n * liTheta t) * rfun t)) volume a b := fun n _ =>
    ((intervalIntegrable_rfun ha hb).continuousOn_mul
      (Continuous.continuousOn (by have := continuous_liTheta; fun_prop))).const_mul _
  have e : ∀ n ∈ range (⌊(B : ℝ) * N⌋₊ + 1), φ (n / N) * (∫ t in a..b, Real.cos (n * liTheta t) * rfun t) =
      ∫ t in a..b, φ (n / N) * (Real.cos (n * liTheta t) * rfun t) := fun n _ =>
    (intervalIntegral.integral_const_mul _ _).symm
  rw [Finset.sum_congr rfl e, ← intervalIntegral.integral_finset_sum hi]
  refine intervalIntegral.integral_congr fun t _ => ?_
  simp only [F0, Finset.sum_mul]
  refine Finset.sum_congr rfl fun n _ => ?_
  ring

lemma sum_mul_integral_cos_log {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) :
    ∑ n ∈ range (⌊(B : ℝ) * N⌋₊ + 1), φ (n / N) *
        (1 / Real.pi * ∫ t in a..b, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))) =
      1 / Real.pi * ∫ t in a..b, F0 φ B N t * Real.log (t / (2 * Real.pi)) := by
  have hi : ∀ n ∈ range (⌊(B : ℝ) * N⌋₊ + 1), IntervalIntegrable
      (fun t => φ (n / N) * (Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi)))) volume a b := fun n _ =>
    ((intervalIntegrable_log_div_two_pi ha hb).continuousOn_mul
      (Continuous.continuousOn (by have := continuous_liTheta; fun_prop))).const_mul _
  have e : ∀ n ∈ range (⌊(B : ℝ) * N⌋₊ + 1), φ (n / N) *
      (1 / Real.pi * ∫ t in a..b, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))) =
      1 / Real.pi * ∫ t in a..b, φ (n / N) * (Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))) := by
    intro n _; rw [intervalIntegral.integral_const_mul]; ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, ← intervalIntegral.integral_finset_sum hi]
  congr 1
  refine intervalIntegral.integral_congr fun t _ => ?_
  simp only [F0, Finset.sum_mul]
  refine Finset.sum_congr rfl fun n _ => ?_
  ring

lemma sum_mul_aCoef {a b : ℝ} :
    ∑ n ∈ range (⌊(B : ℝ) * N⌋₊ + 1), φ (n / N) * (2 * n * aCoef a b n) = ∫ t in a..b, zetaS t * Gf φ B N t := by
  have hc : Continuous fun t : ℝ => 1 / (t ^ 2 + 1 / 4) := continuous_inv_sq_add
  have hi : ∀ n ∈ range (⌊(B : ℝ) * N⌋₊ + 1), IntervalIntegrable
      (fun t => φ (n / N) * (2 * n) * (zetaS t / (t ^ 2 + 1 / 4) * Real.sin (n * liTheta t))) volume a b := by
    intro n _
    apply IntervalIntegrable.const_mul
    have hcont : Continuous fun t : ℝ => 1 / (t ^ 2 + 1 / 4) * Real.sin (n * liTheta t) :=
      hc.mul (by have := continuous_liTheta; fun_prop)
    have h1 := (intervalIntegrable_zetaS a b).mul_continuousOn hcont.continuousOn
    refine h1.congr_uIoo fun t _ => ?_
    simp only [div_eq_mul_one_div (zetaS t)]; ring
  have e : ∀ n ∈ range (⌊(B : ℝ) * N⌋₊ + 1), φ (n / N) * (2 * n * aCoef a b n) =
      ∫ t in a..b, φ (n / N) * (2 * n) * (zetaS t / (t ^ 2 + 1 / 4) * Real.sin (n * liTheta t)) := by
    intro n _; unfold aCoef; rw [intervalIntegral.integral_const_mul]; ring
  rw [Finset.sum_congr rfl e, ← intervalIntegral.integral_finset_sum hi]
  refine intervalIntegral.integral_congr fun t _ => ?_
  simp only [Gf, F1, Finset.mul_sum]
  refine Finset.sum_congr rfl fun n _ => ?_
  have : 0 < t ^ 2 + 1 / 4 := by positivity
  field_simp

end

end LeanLi
