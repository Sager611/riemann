import Riemann.Deep.QRH.Bridge
import Riemann.Deep.WeilEF.Z23.FromPNTPlus.StrongPNTPrefix

/-!
# L(1, χ) ≥ e^{-101} / log q for every nonprincipal Dirichlet character (rung 3)

From OpenAI's zero-free half-plane Re s > 7/8 for all Dirichlet L-functions
(`OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re`, Riemann/Deep/QRH/) we get the
effective lower bound ‖L(1, χ)‖ ≥ e^{-101}/log q for every χ ≠ 1 mod q ≥ 3
(`norm_LFunction_one_ge`). For a real character χ ≠ 1, L(1, χ) is real and positive, so
L(1, χ) ≥ e^{-101}/log q (`LFunction_one_re_ge`, `LFunction_one_lower_real_primitive_re`): no
Siegel-type smallness of L(1, χ), with an explicit constant.

Steps.
1. Partial sums: ‖Σ_{k ≤ n} χ(k)‖ ≤ q (`norm_sum_Icc_chi_le`).
2. Continuation by partial summation: L(s, χ) = s ∫₁^∞ A(t) t^{-s-1} dt on Re s > 0, A(t) = Σ_{k≤t} χ(k)
   (Mathlib `LSeries_eq_mul_integral` on Re s > 1, then the identity theorem; the integral is a
   Mellin transform, hence holomorphic). So ‖L(s, χ)‖ ≤ q‖s‖/Re s (`norm_LFunction_le`).
3. ‖L(σ, χ)‖ ≥ (σ - 1)/σ for σ > 1, from L(s,χ)·L(s,χμ) = 1 and Σ n^{-σ} ≤ σ/(σ - 1).
4. Borel–Carathéodory for J = log L(c₀ + ·, χ) - log L(c₀, χ) on the disc |z| ≤ a + 1/10 around
   c₀ = 1 + a, a = 1/log q (zero-free: Re ≥ 9/10 > 7/8). Re J ≤ log(4q) - log(a/2) ≤ 2 log q + 3,
   so ‖J(-a)‖ ≤ 20 a (2 log q + 3) ≤ 100 and ‖L(1)‖ ≥ ‖L(1 + a)‖ e^{-100} ≥ (a/2) e^{-100}.
-/

open Complex Filter Topology Set MeasureTheory Asymptotics
open scoped Real

namespace Fayed.QRHCor

noncomputable section

private lemma one_lt_log_three' : 1 < Real.log 3 := by
  rw [Real.lt_log_iff_exp_lt (by norm_num)]
  have := Real.exp_one_lt_d9
  norm_num at this ⊢
  linarith

section Partial

variable {q : ℕ} [NeZero q]

/-- A full period of a nonprincipal character sums to 0, from any starting point. -/
lemma sum_range_shift_eq_zero {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) (m : ℕ) :
    ∑ k ∈ Finset.range q, χ ((m + k : ℕ) : ZMod q) = 0 := by
  have h : ∑ k ∈ Finset.range q, χ ((m + k : ℕ) : ZMod q) =
      ∑ x : ZMod q, χ ((m : ZMod q) + x) := by
    refine Finset.sum_nbij' (fun k => (k : ZMod q)) (fun x => x.val) ?_ ?_ ?_ ?_ ?_
    · intro a _; exact Finset.mem_univ _
    · intro a _; exact Finset.mem_range.mpr (ZMod.val_lt a)
    · intro a ha; exact ZMod.val_cast_of_lt (Finset.mem_range.mp ha)
    · intro a _; exact ZMod.natCast_zmod_val a
    · intro a _; simp [Nat.cast_add]
  rw [h]
  calc ∑ x : ZMod q, χ ((m : ZMod q) + x) = ∑ x : ZMod q, χ x :=
        Fintype.sum_equiv (Equiv.addLeft (m : ZMod q)) _ _ (fun x => by simp)
    _ = 0 := MulChar.sum_eq_zero_of_ne_one hχ

lemma norm_sum_range_chi_le {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) (n : ℕ) :
    ‖∑ k ∈ Finset.range n, χ (k : ZMod q)‖ ≤ (q : ℝ) - 1 := by
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    rcases lt_or_ge n q with h | h
    · have hn : (n : ℝ) + 1 ≤ q := by exact_mod_cast h
      calc ‖∑ k ∈ Finset.range n, χ (k : ZMod q)‖ ≤ ∑ k ∈ Finset.range n, ‖χ (k : ZMod q)‖ :=
            norm_sum_le _ _
        _ ≤ ∑ k ∈ Finset.range n, (1 : ℝ) :=
            Finset.sum_le_sum (fun k _ => DirichletCharacter.norm_le_one χ _)
        _ = n := by simp
        _ ≤ (q : ℝ) - 1 := by linarith
    · obtain ⟨m, rfl⟩ : ∃ m, n = m + q := ⟨n - q, by omega⟩
      have hsplit : ∑ k ∈ Finset.range (m + q), χ (k : ZMod q) =
          ∑ k ∈ Finset.range m, χ (k : ZMod q) + ∑ k ∈ Finset.range q, χ ((m + k : ℕ) : ZMod q) :=
        Finset.sum_range_add _ _ _
      rw [hsplit, sum_range_shift_eq_zero hχ, add_zero]
      exact ih m (by have := NeZero.pos q; omega)

lemma sum_Icc_eq_sum_range_sub (f : ℕ → ℂ) (n : ℕ) :
    ∑ k ∈ Finset.Icc 1 n, f k = ∑ k ∈ Finset.range (n + 1), f k - f 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_Icc_succ_top (by omega), ih, Finset.sum_range_succ _ (n + 1)]
    ring

/-- Partial sums of a nonprincipal character are bounded by the modulus. -/
lemma norm_sum_Icc_chi_le {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) (n : ℕ) :
    ‖∑ k ∈ Finset.Icc 1 n, χ (k : ZMod q)‖ ≤ q := by
  rw [sum_Icc_eq_sum_range_sub (fun k => χ (k : ZMod q))]
  calc _ ≤ ‖∑ k ∈ Finset.range (n + 1), χ (k : ZMod q)‖ + ‖χ ((0 : ℕ) : ZMod q)‖ :=
        norm_sub_le _ _
    _ ≤ ((q : ℝ) - 1) + 1 :=
        add_le_add (norm_sum_range_chi_le hχ _) (DirichletCharacter.norm_le_one χ _)
    _ = q := by ring

/-- A(t) = Σ_{1 ≤ k ≤ t} χ(k). -/
def Achi (χ : DirichletCharacter ℂ q) (t : ℝ) : ℂ := ∑ k ∈ Finset.Icc 1 ⌊t⌋₊, χ (k : ZMod q)

lemma norm_Achi_le {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) (t : ℝ) : ‖Achi χ t‖ ≤ q :=
  norm_sum_Icc_chi_le hχ _

omit [NeZero q] in
lemma measurable_Achi (χ : DirichletCharacter ℂ q) : Measurable (Achi χ) :=
  (measurable_from_nat (f := fun n : ℕ => ∑ k ∈ Finset.Icc 1 n, χ (k : ZMod q))).comp
    measurable_id.nat_floor

/-- A on (1, ∞), zero elsewhere. -/
def Fchi (χ : DirichletCharacter ℂ q) : ℝ → ℂ := (Ioi (1 : ℝ)).indicator (Achi χ)

lemma norm_Fchi_le {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) (t : ℝ) : ‖Fchi χ t‖ ≤ q := by
  unfold Fchi indicator
  split_ifs
  · exact norm_Achi_le hχ t
  · simp

omit [NeZero q] in
lemma measurable_Fchi (χ : DirichletCharacter ℂ q) : Measurable (Fchi χ) :=
  (measurable_Achi χ).indicator measurableSet_Ioi

/-- I(s) = ∫₁^∞ A(t) t^{-s-1} dt, as a Mellin transform. -/
def Ichi (χ : DirichletCharacter ℂ q) (s : ℂ) : ℂ := mellin (Fchi χ) (-s)

omit [NeZero q] in
lemma mellin_Fchi (χ : DirichletCharacter ℂ q) (z : ℂ) :
    mellin (Fchi χ) z = ∫ t in Ioi (1 : ℝ), (t : ℂ) ^ (z - 1) * Achi χ t := by
  rw [mellin]
  have : (fun t : ℝ => (t : ℂ) ^ (z - 1) • Fchi χ t) =
      (Ioi (1 : ℝ)).indicator (fun t => (t : ℂ) ^ (z - 1) * Achi χ t) := by
    ext t; simp only [Fchi, smul_eq_mul, indicator]; split_ifs <;> simp
  rw [this, setIntegral_indicator measurableSet_Ioi,
    inter_eq_right.mpr (Ioi_subset_Ioi zero_le_one)]

omit [NeZero q] in
lemma Ichi_eq (χ : DirichletCharacter ℂ q) (s : ℂ) :
    Ichi χ s = ∫ t in Ioi (1 : ℝ), (t : ℂ) ^ (-s - 1) * Achi χ t :=
  mellin_Fchi χ (-s)

lemma kernel_chi_norm_le {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) {s : ℂ} {t : ℝ} (ht : 1 < t) :
    ‖(t : ℂ) ^ (-s - 1) * Achi χ t‖ ≤ q * t ^ (-s.re - 1) := by
  have ht0 : 0 < t := lt_trans zero_lt_one ht
  rw [norm_mul, Complex.norm_cpow_eq_rpow_re_of_pos ht0]
  have h1 : (-s - 1).re = -s.re - 1 := by simp
  rw [h1, mul_comm]
  exact mul_le_mul_of_nonneg_right (norm_Achi_le hχ t) (Real.rpow_nonneg ht0.le _)

lemma norm_Ichi_le {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) {s : ℂ} (hs : 0 < s.re) :
    ‖Ichi χ s‖ ≤ q / s.re := by
  rw [Ichi_eq]
  refine (norm_integral_le_of_norm_le
    (Integrable.const_mul (integrableOn_Ioi_rpow_of_lt (a := -s.re - 1) (by linarith) one_pos)
      (q : ℝ))
    ((ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun t ht =>
      kernel_chi_norm_le hχ ht))).trans (le_of_eq ?_)
  rw [integral_const_mul, integral_Ioi_rpow_of_lt (by linarith) one_pos, Real.one_rpow]
  have : -s.re - 1 + 1 = -s.re := by ring
  rw [this, neg_div_neg_eq, mul_one_div]

lemma locallyIntegrable_Fchi {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) :
    LocallyIntegrable (Fchi χ) := by
  intro x
  refine ⟨Ioo (x - 1) (x + 1), Ioo_mem_nhds (by linarith) (by linarith), ?_⟩
  exact Measure.integrableOn_of_bounded (M := q) measure_Ioo_lt_top.ne
    (measurable_Fchi χ).aestronglyMeasurable (ae_of_all _ (norm_Fchi_le hχ))

lemma differentiableAt_Ichi {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) {s : ℂ} (hs : 0 < s.re) :
    DifferentiableAt ℂ (Ichi χ) s := by
  have hm : DifferentiableAt ℂ (mellin (Fchi χ)) (-s) := by
    refine mellin_differentiableAt_of_isBigO_rpow (a := 0) (b := (-s).re - 1)
      ((locallyIntegrable_Fchi hχ).locallyIntegrableOn _) ?_ (by simp; linarith) ?_ (by linarith)
    · refine IsBigO.of_bound (q : ℝ) ?_
      filter_upwards [eventually_gt_atTop (0 : ℝ)] with t _
      rw [neg_zero, Real.rpow_zero, norm_one, mul_one]
      exact norm_Fchi_le hχ t
    · have hz : Fchi χ =ᶠ[𝓝[>] 0] fun _ => (0 : ℂ) := by
        filter_upwards [nhdsWithin_le_nhds (Iio_mem_nhds (zero_lt_one : (0 : ℝ) < 1))] with t ht
        simp only [Fchi, indicator, mem_Ioi]
        rw [ite_eq_right_iff.mpr (fun h => absurd (lt_trans h ht) (lt_irrefl 1))]
      exact (isBigO_zero _ _).congr' hz.symm EventuallyEq.rfl
  exact hm.comp s differentiableAt_id.neg

/-- Partial summation for Re s > 1. -/
lemma LFunction_eq_mul_Ichi_of_one_lt {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) {s : ℂ}
    (hs : 1 < s.re) : χ.LFunction s = s * Ichi χ s := by
  have hO : (fun n : ℕ => ∑ k ∈ Finset.Icc 1 n, (fun m : ℕ => χ (m : ZMod q)) k) =O[atTop]
      fun n : ℕ => (n : ℝ) ^ (0 : ℝ) := by
    refine IsBigO.of_bound (q : ℝ) (Eventually.of_forall fun n => ?_)
    rw [Real.rpow_zero, norm_one, mul_one]
    exact norm_sum_Icc_chi_le hχ n
  rw [DirichletCharacter.LFunction_eq_LSeries χ hs,
    LSeries_eq_mul_integral _ (r := 0) le_rfl (by linarith)
      (DirichletCharacter.LSeriesSummable_of_one_lt_re χ hs) hO, Ichi_eq]
  congr 1
  refine setIntegral_congr_fun measurableSet_Ioi (fun t _ => ?_)
  simp only [Achi]
  rw [show -(s + 1) = -s - 1 by ring]
  ring

/-- **Continuation by partial summation:** L(s, χ) = s ∫₁^∞ A(t) t^{-s-1} dt on Re s > 0. -/
theorem LFunction_eq_mul_Ichi {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) {s : ℂ} (hs : 0 < s.re) :
    χ.LFunction s = s * Ichi χ s := by
  have hUo : IsOpen {z : ℂ | 0 < z.re} := isOpen_lt continuous_const Complex.continuous_re
  have hf : AnalyticOnNhd ℂ χ.LFunction {z : ℂ | 0 < z.re} := fun z _ =>
    (DirichletCharacter.differentiable_LFunction hχ).analyticAt z
  have hg : AnalyticOnNhd ℂ (fun z => z * Ichi χ z) {z : ℂ | 0 < z.re} := by
    refine DifferentiableOn.analyticOnNhd (fun z hz => ?_) hUo
    exact (DifferentiableAt.mul differentiableAt_id
      (differentiableAt_Ichi hχ hz)).differentiableWithinAt
  have hpc : IsPreconnected {z : ℂ | 0 < z.re} := (convex_halfSpace_re_gt 0).isPreconnected
  have h2 : (2 : ℂ) ∈ {z : ℂ | 0 < z.re} := by show (0 : ℝ) < (2 : ℂ).re; norm_num
  have hfg : χ.LFunction =ᶠ[𝓝 (2 : ℂ)] fun z => z * Ichi χ z := by
    have : {z : ℂ | 1 < z.re} ∈ 𝓝 (2 : ℂ) :=
      (isOpen_lt continuous_const Complex.continuous_re).mem_nhds (by norm_num)
    filter_upwards [this] with z hz using LFunction_eq_mul_Ichi_of_one_lt hχ hz
  exact hf.eqOn_of_preconnected_of_eventuallyEq hg hpc h2 hfg hs

/-- **Convexity-type bound:** ‖L(s, χ)‖ ≤ q‖s‖/Re s on Re s > 0, for χ ≠ 1. -/
theorem norm_LFunction_le {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1) {s : ℂ} (hs : 0 < s.re) :
    ‖χ.LFunction s‖ ≤ q * ‖s‖ / s.re := by
  rw [LFunction_eq_mul_Ichi hχ hs, norm_mul]
  calc ‖s‖ * ‖Ichi χ s‖ ≤ ‖s‖ * (q / s.re) :=
        mul_le_mul_of_nonneg_left (norm_Ichi_le hχ hs) (norm_nonneg _)
    _ = q * ‖s‖ / s.re := by ring

/-- An L-series with coefficients of norm ≤ 1 is ≤ σ/(σ - 1) in norm on Re s = σ > 1. -/
lemma norm_LSeries_le_re {f : ℕ → ℂ} (hf : ∀ n, ‖f n‖ ≤ 1) {s : ℂ} (hs : 1 < s.re) :
    ‖LSeries f s‖ ≤ s.re / (s.re - 1) := by
  set σ := s.re with hσdef
  have hσ : 1 < σ := hs
  set b : ℕ → ℝ := fun n => if n = 0 then 0 else 1 / (n : ℝ) ^ σ with hb
  have hle : ∀ n, ‖LSeries.term f s n‖ ≤ b n := by
    intro n
    rw [LSeries.norm_term_eq]
    simp only [hb]
    split_ifs
    · exact le_rfl
    · exact div_le_div_of_nonneg_right (hf n) (by positivity)
  have hbs : Summable b := by
    refine Summable.of_nonneg_of_le (fun n => ?_) (fun n => ?_)
      (Real.summable_one_div_nat_rpow.mpr hσ)
    · simp only [hb]
      split_ifs
      · exact le_rfl
      · positivity
    · simp only [hb]
      split_ifs
      · positivity
      · exact le_rfl
  have hsum : Summable (fun n => ‖LSeries.term f s n‖) :=
    Summable.of_nonneg_of_le (fun _ => norm_nonneg _) hle hbs
  have htb : ∑' n, b n = ∑' n : ℕ, 1 / ((n : ℝ) + 1) ^ σ := by
    rw [hbs.tsum_eq_zero_add]
    simp [hb]
  have hz : ∑' n : ℕ, 1 / ((n : ℝ) + 1) ^ σ ≤ σ / (σ - 1) := by
    have h := ZetaAsymptotics.zeta_limit_aux1 hσ
    have hT : 0 ≤ ZetaAsymptotics.termTSum σ :=
      tsum_nonneg (fun n => ZetaAsymptotics.term_nonneg (n + 1) σ)
    have hσ1 : (σ - 1) ≠ 0 := (by linarith : (0 : ℝ) < σ - 1).ne'
    have e : σ / (σ - 1) = 1 + 1 / (σ - 1) := by field_simp; ring
    have hσT := mul_nonneg (by linarith : (0 : ℝ) ≤ σ) hT
    rw [e]
    linarith
  calc ‖LSeries f s‖ = ‖∑' n, LSeries.term f s n‖ := rfl
    _ ≤ ∑' n, ‖LSeries.term f s n‖ := norm_tsum_le_tsum_norm hsum
    _ ≤ ∑' n, b n := hsum.tsum_le_tsum hle hbs
    _ = ∑' n : ℕ, 1 / ((n : ℝ) + 1) ^ σ := htb
    _ ≤ σ / (σ - 1) := hz

/-- ‖L(s, χ)‖ ≥ (σ - 1)/σ on Re s = σ > 1 (from L(s,χ) · L(s,χμ) = 1). -/
theorem norm_LFunction_ge_of_one_lt_re (χ : DirichletCharacter ℂ q) {s : ℂ} (hs : 1 < s.re) :
    (s.re - 1) / s.re ≤ ‖χ.LFunction s‖ := by
  set σ := s.re with hσdef
  have hσ : 1 < σ := hs
  have h1 := DirichletCharacter.LSeries.mul_mu_eq_one χ hs
  rw [← DirichletCharacter.LFunction_eq_LSeries χ hs] at h1
  have key : ∀ x : ℂ, ‖x‖ ≤ σ / (σ - 1) → χ.LFunction s * x = 1 →
      (σ - 1) / σ ≤ ‖χ.LFunction s‖ := by
    intro x hx hLx
    have hn := congrArg norm hLx
    rw [norm_mul, norm_one] at hn
    have hσ1 : 0 < σ - 1 := by linarith
    have hx' : (σ - 1) * ‖x‖ ≤ σ := by
      rw [le_div_iff₀ hσ1] at hx
      linarith
    rw [div_le_iff₀ (by linarith : (0 : ℝ) < σ)]
    nlinarith [mul_le_mul_of_nonneg_left hx' (norm_nonneg (χ.LFunction s)), norm_nonneg x]
  refine key _ (norm_LSeries_le_re (fun n => ?_) hs) h1
  simp only [Pi.mul_apply, norm_mul, Complex.norm_intCast]
  have ha := DirichletCharacter.norm_le_one χ (n : ZMod q)
  have hb : |((ArithmeticFunction.moebius n : ℤ) : ℝ)| ≤ 1 := by
    exact_mod_cast ArithmeticFunction.abs_moebius_le_one
  calc _ ≤ 1 * 1 := mul_le_mul ha hb (abs_nonneg _) zero_le_one
    _ = 1 := one_mul 1

/-- ‖L(σ, χ)‖ ≥ (σ - 1)/σ for real σ > 1. -/
theorem norm_LFunction_ge_of_one_lt (χ : DirichletCharacter ℂ q) {σ : ℝ} (hσ : 1 < σ) :
    (σ - 1) / σ ≤ ‖χ.LFunction (σ : ℂ)‖ := by
  have h := norm_LFunction_ge_of_one_lt_re χ (s := (σ : ℂ)) (by rw [Complex.ofReal_re]; exact hσ)
  rwa [Complex.ofReal_re] at h

end Partial

/-- **Rung 3: effective lower bound at s = 1.** For every nonprincipal Dirichlet character χ
mod q ≥ 3, ‖L(1, χ)‖ ≥ e^{-101}/log q. (Uses OpenAI's 7/8 zero-free half-plane.) -/
theorem norm_LFunction_one_ge {q : ℕ} [NeZero q] (hq : 3 ≤ q) {χ : DirichletCharacter ℂ q}
    (hχ : χ ≠ 1) : Real.exp (-101) / Real.log q ≤ ‖χ.LFunction 1‖ := by
  set L := Real.log q with hLdef
  have hq0 : (0 : ℝ) < q := by exact_mod_cast NeZero.pos q
  have hL1 : 1 < L :=
    lt_of_lt_of_le one_lt_log_three' (Real.log_le_log (by norm_num) (by exact_mod_cast hq))
  have hL0 : 0 < L := by linarith
  set a : ℝ := 1 / L with hadef
  have ha0 : 0 < a := by positivity
  have ha1 : a < 1 := by rw [hadef, div_lt_one hL0]; exact hL1
  set R : ℝ := a + 1 / 10 with hRdef
  have hR0 : 0 < R := by positivity
  have haR : a < R := by linarith
  set c₀ : ℂ := ((1 + a : ℝ) : ℂ) with hc₀
  have hc₀re : c₀.re = 1 + a := by simp [hc₀]
  have hc₀norm : ‖c₀‖ = 1 + a := by
    rw [hc₀, Complex.norm_real, Real.norm_of_nonneg (by linarith)]
  -- geometry of the disc |z| ≤ R around c₀
  have hgeom : ∀ z ∈ Metric.closedBall (0 : ℂ) R, 9 / 10 ≤ (c₀ + z).re ∧ ‖c₀ + z‖ ≤ 31 / 10 := by
    intro z hz
    rw [mem_closedBall_zero_iff] at hz
    have hre := (abs_le.mp ((Complex.abs_re_le_norm z).trans hz)).1
    constructor
    · rw [add_re, hc₀re]; linarith
    · calc ‖c₀ + z‖ ≤ ‖c₀‖ + ‖z‖ := norm_add_le _ _
        _ ≤ (1 + a) + R := by rw [hc₀norm]; linarith
        _ ≤ 31 / 10 := by linarith
  have hanal : AnalyticOnNhd ℂ (fun z => χ.LFunction (c₀ + z)) (Metric.closedBall 0 R) := by
    intro z _
    exact ((DirichletCharacter.differentiable_LFunction hχ).analyticAt (c₀ + z)).comp
      (analyticAt_const.add analyticAt_id)
  have hB0 : ∀ z ∈ Metric.closedBall (0 : ℂ) R, χ.LFunction (c₀ + z) ≠ 0 := fun z hz =>
    OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re χ
      (by linarith [(hgeom z hz).1]) (fun h => hχ h.1)
  obtain ⟨J, hJan, hJ0, -, hJre⟩ := LogOfAnalyticFunction (r := a) ha0 haR hanal hB0
  -- lower bound at the centre
  have hcen : a / 2 ≤ ‖χ.LFunction c₀‖ := by
    have h := norm_LFunction_ge_of_one_lt χ (σ := 1 + a) (by linarith)
    have e : (1 + a - 1) / (1 + a) = a / (1 + a) := by ring
    rw [e] at h
    refine le_trans ?_ h
    rw [div_le_div_iff₀ (by norm_num) (by linarith)]
    nlinarith
  -- upper bound on the disc
  have hup : ∀ z ∈ Metric.ball (0 : ℂ) R, ‖χ.LFunction (c₀ + z)‖ ≤ 4 * q := by
    intro z hz
    obtain ⟨g1, g2⟩ := hgeom z (Metric.ball_subset_closedBall hz)
    have h := norm_LFunction_le hχ (s := c₀ + z) (by linarith)
    refine h.trans ?_
    rw [div_le_iff₀ (by linarith)]
    have e1 : (q : ℝ) * ‖c₀ + z‖ ≤ q * (31 / 10) := mul_le_mul_of_nonneg_left g2 hq0.le
    have e2 : (q : ℝ) * (9 / 10) ≤ q * (c₀ + z).re := mul_le_mul_of_nonneg_left g1 hq0.le
    nlinarith
  -- Re J ≤ M on the disc
  set M : ℝ := Real.log (4 * q) - Real.log (a / 2) with hMdef
  have hJM : ∀ z ∈ Metric.ball (0 : ℂ) R, (J z).re ≤ M := by
    intro z hz
    rw [← hJre z hz]
    have hBz : 0 < ‖χ.LFunction (c₀ + z)‖ :=
      norm_pos_iff.mpr (hB0 z (Metric.ball_subset_closedBall hz))
    have l1 := Real.log_le_log hBz (hup z hz)
    have l2 := Real.log_le_log (by positivity) hcen
    show Real.log ‖χ.LFunction (c₀ + z)‖ - Real.log ‖χ.LFunction (c₀ + 0)‖ ≤ M
    rw [add_zero]
    linarith
  have hM0 : 0 < M := by
    rw [hMdef, sub_pos]
    apply Real.log_lt_log (by positivity)
    have : (3 : ℝ) ≤ q := by exact_mod_cast hq
    linarith
  have hMle : M ≤ 2 * L + 3 := by
    have e1 : Real.log (4 * q) = Real.log 4 + L := Real.log_mul (by norm_num) hq0.ne'
    have e2 : Real.log (a / 2) = Real.log a - Real.log 2 := Real.log_div ha0.ne' (by norm_num)
    have e3 : Real.log a = -Real.log L := by rw [hadef, one_div, Real.log_inv]
    have e4 : Real.log 4 = 2 * Real.log 2 := by
      rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; norm_num
    have h2 := Real.log_two_lt_d9
    have h3 := Real.log_le_sub_one_of_pos hL0
    rw [hMdef, e1, e2, e3, e4]
    linarith
  -- Borel–Carathéodory at z = -a
  have hz : (-(a : ℂ)) ∈ Metric.closedBall (0 : ℂ) a := by
    rw [mem_closedBall_zero_iff, norm_neg, Complex.norm_real, Real.norm_of_nonneg ha0.le]
  have hBC := borelCaratheodory' hM0 hR0 haR hJan.analyticOn hJ0 hJM hz
  have hRa : R - a = 1 / 10 := by rw [hRdef]; ring
  have haL : a * L = 1 := by rw [hadef]; exact one_div_mul_cancel hL0.ne'
  have hJ100 : ‖J (-(a : ℂ))‖ ≤ 100 := by
    refine hBC.trans ?_
    rw [hRa, div_le_iff₀ (by norm_num)]
    nlinarith [mul_le_mul_of_nonneg_right hMle ha0.le]
  -- conclude
  have hre := hJre (-(a : ℂ)) (Metric.closedBall_subset_ball haR hz)
  simp only [add_zero] at hre
  have hpt : c₀ + -(a : ℂ) = 1 := by rw [hc₀]; push_cast; ring
  rw [hpt] at hre
  have hreJ : -100 ≤ (J (-(a : ℂ))).re := by
    have h1 := Complex.abs_re_le_norm (J (-(a : ℂ)))
    have h2 := neg_abs_le (J (-(a : ℂ))).re
    linarith
  have h1pos : 0 < ‖χ.LFunction 1‖ :=
    norm_pos_iff.mpr (DirichletCharacter.LFunction_apply_one_ne_zero hχ)
  have hlog : Real.log (a / 2) - 100 ≤ Real.log ‖χ.LFunction 1‖ := by
    have := Real.log_le_log (by positivity) hcen
    linarith
  have hexp : (a / 2) * Real.exp (-100) ≤ ‖χ.LFunction 1‖ := by
    calc (a / 2) * Real.exp (-100) = Real.exp (Real.log (a / 2) - 100) := by
          rw [Real.exp_sub, Real.exp_log (by positivity), Real.exp_neg]
          ring
      _ ≤ Real.exp (Real.log ‖χ.LFunction 1‖) := Real.exp_le_exp.mpr hlog
      _ = ‖χ.LFunction 1‖ := Real.exp_log h1pos
  have he : Real.exp (-101) ≤ Real.exp (-100) / 2 := by
    have e : Real.exp (-101) = Real.exp (-100) * Real.exp (-1) := by
      rw [← Real.exp_add]; norm_num
    have h2 : Real.exp (-1) ≤ 1 / 2 := by
      have h3 : 2 ≤ Real.exp 1 := by have := Real.add_one_le_exp (1 : ℝ); linarith
      have h4 : Real.exp (-1) * Real.exp 1 = 1 := by rw [← Real.exp_add]; norm_num
      have h5 := Real.exp_pos (-1)
      nlinarith
    rw [e]
    have := Real.exp_pos (-100)
    nlinarith
  calc Real.exp (-101) / L ≤ (Real.exp (-100) / 2) / L := div_le_div_of_nonneg_right he hL0.le
    _ = (a / 2) * Real.exp (-100) := by rw [hadef]; ring
    _ ≤ ‖χ.LFunction 1‖ := hexp

/-- The ∃-form for real primitive characters (the brief's rung 3): there is an absolute,
effective c > 0 (c = e^{-101}) with ‖L(1, χ)‖ ≥ c/log q. -/
theorem LFunction_one_lower_real_primitive :
    ∃ c : ℝ, 0 < c ∧ ∀ (q : ℕ) [NeZero q], 3 ≤ q → ∀ χ : DirichletCharacter ℂ q,
      χ.IsPrimitive → χ ≠ 1 → (∀ a : ZMod q, (χ a).im = 0) →
        c / Real.log q ≤ ‖χ.LFunction 1‖ :=
  ⟨Real.exp (-101), Real.exp_pos _, fun _ _ hq _ _ hχ _ => norm_LFunction_one_ge hq hχ⟩

section RealChar

variable {q : ℕ} [NeZero q]

/-- For a real character, L(σ, χ) is real for real σ > 1. -/
lemma LFunction_im_eq_zero_of_one_lt {χ : DirichletCharacter ℂ q}
    (hreal : ∀ a : ZMod q, (χ a).im = 0) {σ : ℝ} (hσ : 1 < σ) :
    (χ.LFunction (σ : ℂ)).im = 0 := by
  rw [DirichletCharacter.LFunction_eq_LSeries χ (by rw [Complex.ofReal_re]; exact hσ)]
  have hterm : ∀ n : ℕ, (LSeries.term (fun n : ℕ => χ (n : ZMod q)) (σ : ℂ) n).im = 0 := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · simp
    · rw [LSeries.term_of_ne_zero hn, ← Complex.ofReal_natCast,
        ← Complex.ofReal_cpow (Nat.cast_nonneg (α := ℝ) n) σ, Complex.div_ofReal_im, hreal,
        zero_div]
  have hsum : ∑' n, LSeries.term (fun n : ℕ => χ (n : ZMod q)) (σ : ℂ) n =
      ∑' n, (((LSeries.term (fun n : ℕ => χ (n : ZMod q)) (σ : ℂ) n).re : ℝ) : ℂ) :=
    tsum_congr (fun n => Complex.ext (by simp) (by simp [hterm n]))
  show (∑' n, LSeries.term (fun n : ℕ => χ (n : ZMod q)) (σ : ℂ) n).im = 0
  rw [hsum, ← Complex.ofReal_tsum, Complex.ofReal_im]

/-- For a real character χ ≠ 1, L(1, χ) is real. -/
lemma LFunction_one_im_eq_zero {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1)
    (hreal : ∀ a : ZMod q, (χ a).im = 0) : (χ.LFunction 1).im = 0 := by
  have hc : Continuous (fun σ : ℝ => (χ.LFunction (σ : ℂ)).im) :=
    Complex.continuous_im.comp
      ((DirichletCharacter.differentiable_LFunction hχ).continuous.comp Complex.continuous_ofReal)
  have h1 : Tendsto (fun σ : ℝ => (χ.LFunction (σ : ℂ)).im) (𝓝[>] 1)
      (𝓝 ((χ.LFunction ((1 : ℝ) : ℂ)).im)) :=
    (hc.tendsto 1).mono_left nhdsWithin_le_nhds
  have h2 : Tendsto (fun σ : ℝ => (χ.LFunction (σ : ℂ)).im) (𝓝[>] 1) (𝓝 0) :=
    tendsto_const_nhds.congr' (eventually_nhdsWithin_of_forall fun σ hσ =>
      (LFunction_im_eq_zero_of_one_lt hreal hσ).symm)
  have := tendsto_nhds_unique h1 h2
  simpa using this

/-- For a real character χ ≠ 1, L(1, χ) > 0 (intermediate value theorem on [1, X], since L is
real and zero-free on [1, ∞) and L(X) → 1 as X → ∞). -/
lemma LFunction_one_re_pos {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1)
    (hreal : ∀ a : ZMod q, (χ a).im = 0) : 0 < (χ.LFunction 1).re := by
  have hc : Continuous (fun σ : ℝ => (χ.LFunction (σ : ℂ)).re) :=
    Complex.continuous_re.comp
      ((DirichletCharacter.differentiable_LFunction hχ).continuous.comp Complex.continuous_ofReal)
  have habs : LSeries.abscissaOfAbsConv (fun n : ℕ => χ (n : ZMod q)) < ⊤ :=
    lt_of_le_of_lt (LSeries.abscissaOfAbsConv_le_of_forall_lt_LSeriesSummable (x := 1)
      (fun y hy => DirichletCharacter.LSeriesSummable_of_one_lt_re χ (by simpa using hy)))
      (EReal.coe_lt_top 1)
  have hlim' : Tendsto (fun x : ℝ => (LSeries (fun n : ℕ => χ (n : ZMod q)) x).re) atTop
      (𝓝 1) := by
    have := (Complex.continuous_re.tendsto _).comp (LSeries.tendsto_atTop habs)
    simpa [Function.comp_def] using this
  obtain ⟨X, hXpos, hX1⟩ := ((hlim'.eventually (lt_mem_nhds (by norm_num : (0 : ℝ) < 1))).and
    (eventually_gt_atTop 1)).exists
  have hgX : 0 < (χ.LFunction (X : ℂ)).re := by
    rw [DirichletCharacter.LFunction_eq_LSeries χ (by rw [Complex.ofReal_re]; exact hX1)]
    exact hXpos
  by_contra hneg
  have hg1 : (χ.LFunction ((1 : ℝ) : ℂ)).re ≤ 0 := by simpa using not_lt.mp hneg
  obtain ⟨σ, hσ, hσ0⟩ := intermediate_value_Icc hX1.le hc.continuousOn ⟨hg1, hgX.le⟩
  have him : (χ.LFunction (σ : ℂ)).im = 0 := by
    rcases eq_or_lt_of_le hσ.1 with h | h
    · subst h
      simpa using LFunction_one_im_eq_zero hχ hreal
    · exact LFunction_im_eq_zero_of_one_lt hreal h
  have hzero : χ.LFunction (σ : ℂ) = 0 :=
    Complex.ext (by simpa using hσ0) (by simpa using him)
  exact DirichletCharacter.LFunction_ne_zero_of_one_le_re χ (Or.inl hχ)
    (by rw [Complex.ofReal_re]; exact hσ.1) hzero

end RealChar

/-- **Rung 3, real form:** for a real character χ ≠ 1 mod q ≥ 3, L(1, χ) is a real number
with L(1, χ) ≥ e^{-101}/log q. -/
theorem LFunction_one_re_ge {q : ℕ} [NeZero q] (hq : 3 ≤ q) {χ : DirichletCharacter ℂ q}
    (hχ : χ ≠ 1) (hreal : ∀ a : ZMod q, (χ a).im = 0) :
    (χ.LFunction 1).im = 0 ∧ Real.exp (-101) / Real.log q ≤ (χ.LFunction 1).re := by
  have him := LFunction_one_im_eq_zero hχ hreal
  have hre := LFunction_one_re_pos hχ hreal
  refine ⟨him, ?_⟩
  have hnorm : ‖χ.LFunction 1‖ = (χ.LFunction 1).re := by
    calc ‖χ.LFunction 1‖ = ‖(((χ.LFunction 1).re : ℝ) : ℂ)‖ := by
          congr 1
          exact Complex.ext (by simp) (by simp [him])
      _ = (χ.LFunction 1).re := Complex.norm_of_nonneg hre.le
  rw [← hnorm]
  exact norm_LFunction_one_ge hq hχ

/-- The brief's rung 3 verbatim (real primitive χ, real part): ∃ c > 0 (c = e^{-101}) with
c/log q ≤ L(1, χ). Compare OpenAI's siegel_LFunction_one_lower_bound: c_ε q^{-ε} ≤ L(1, χ). -/
theorem LFunction_one_lower_real_primitive_re :
    ∃ c : ℝ, 0 < c ∧ ∀ (q : ℕ) [NeZero q], 3 ≤ q → ∀ χ : DirichletCharacter ℂ q,
      χ.IsPrimitive → χ ≠ 1 → (∀ a : ZMod q, (χ a).im = 0) →
        c / Real.log q ≤ (χ.LFunction 1).re :=
  ⟨Real.exp (-101), Real.exp_pos _,
    fun _ _ hq _ _ hχ hreal => (LFunction_one_re_ge hq hχ hreal).2⟩

end

end Fayed.QRHCor
