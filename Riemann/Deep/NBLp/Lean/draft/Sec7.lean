/- UNCOMPILED DRAFT (nblplean, 2026-11-12). Not checked by Lean: every compile of it exceeds the 290 MB
   per-process cap (imports alone: Pow.Real 334 MB, LpSeminorm.Defs 343 MB, LSeries.RiemannZeta 342 MB).
   See ../../LEAN.md. Do not import from checked files. -/
import Riemann.Deep.NBLp.Lean.Toolkit
import Mathlib.Analysis.SumIntegralComparisons

/-!
# NBLp in Lean, part 4: §7 objects, A1 identities, [C8] facts, A5, A6 (L^p part), A7

PROOF.md §7. Definitions φ_ε, E_ε, Φ_ε, H_ε, H_{ε,N}, Φ_{ε,N}, R_ε, Ψ_ε and the constants
r, ε₀, C_Γ, K_p, A_p of §1. Proved here: A1's identities on (0,1] (A1b, A1c), the [C8] facts
used (1/ζ(1+ε) ≤ ε, ζ decreasing on (1,∞)), A5 and A6 from their L²_σ inputs (Lemma M),
and A7 (N → ∞ for fixed ε, then ε → 0).
-/

open Complex MeasureTheory Set Filter Topology
open scoped ENNReal

namespace NBLp

/-- μ(k) as a complex number. -/
noncomputable def mu (k : ℕ) : ℂ := (ArithmeticFunction.moebius k : ℂ)

/-- §7 A1: φ_ε(j) = Σ_{d|j} μ(d) d^{−ε}. -/
noncomputable def phiE (ε : ℝ) (j : ℕ) : ℂ := ∑ d ∈ j.divisors, mu d * (d : ℂ) ^ (-(ε : ℂ))

/-- §7 A1: E_ε(y) = Σ_{j≤y} φ_ε(j) − y/ζ(1+ε). -/
noncomputable def EE (ε y : ℝ) : ℂ :=
  ∑ j ∈ Finset.Icc 1 ⌊y⌋₊, phiE ε j - (y : ℂ) / riemannZeta (1 + ε)

/-- §7 A1: Φ_ε(x) = χ₊(x) − E_ε(1/x). -/
noncomputable def PhiE (ε x : ℝ) : ℂ := chiP x - EE ε (1 / x)

/-- §7 A1: H_ε(x) = −E_ε(⌊1/x⌋). -/
noncomputable def HE (ε x : ℝ) : ℂ := -EE ε (⌊1 / x⌋₊ : ℝ)

/-- μ(k) k^{−ε} (1 − k/N). -/
noncomputable def coef (ε : ℝ) (N k : ℕ) : ℂ := mu k * (k : ℂ) ^ (-(ε : ℂ)) * (1 - (k : ℂ) / N)

/-- Coefficients of Theorem A's approximant: χ − Σ cA_k f_k = χ + H_{ε,N}. -/
noncomputable def cA (ε : ℝ) (N k : ℕ) : ℂ := -coef ε N k

/-- §7 A1: H_{ε,N} = Σ_{k=2}^N μ(k) k^{−ε} (1 − k/N) f_k. -/
noncomputable def HEN (ε : ℝ) (N : ℕ) (x : ℝ) : ℂ := ∑ k ∈ Finset.Icc 2 N, coef ε N k * (fk k x : ℂ)

/-- §7 A1: Φ_{ε,N} = χ₊ + Σ_{k=1}^N μ(k) k^{−ε} (1 − k/N) ρ_k. -/
noncomputable def PhiEN (ε : ℝ) (N : ℕ) (x : ℝ) : ℂ :=
  chiP x + ∑ k ∈ Finset.Icc 1 N, coef ε N k * (rho k x : ℂ)

/-- §7 A2: R_ε(s) = ζ(s)/ζ(s+ε) (= 0 at s = 1 − ε). -/
noncomputable def RE (ε : ℝ) (s : ℂ) : ℂ := riemannZeta s * invZeta (s + ε)

/-- §7 A2: Ψ_ε(s) = (1 − R_ε(s))/s. -/
noncomputable def PsiE (ε : ℝ) (s : ℂ) : ℂ := (1 - RE ε s) / s

/-- §0: r = 2p/(2−p). -/
noncomputable def rr (p : ℝ) : ℝ := 2 * p / (2 - p)

/-- §1: ε₀ = min{1/2, (2/p − 1)/2, 1 − 1/p}. -/
noncomputable def eps0 (p : ℝ) : ℝ := min (1 / 2) (min ((2 / p - 1) / 2) (1 - 1 / p))

/-- §1: the integrand of C_Γ. -/
noncomputable def gCG (t : ℝ) : ℝ := (8 + Real.log (2 + |t|)) * (2 + |t|) ^ (1 / 2 : ℝ) * (1 / 4 + t ^ 2)⁻¹

/-- §1: C_Γ = e⁴ ∫ (8 + log(2+|t|)) (2+|t|)^{1/2} (1/4 + t²)^{−1} dt. -/
noncomputable def CG : ℝ := Real.exp 4 * ∫ t : ℝ, gCG t

/-- §1: K_p = 4p/(2−p) + C_Γ/(2π). -/
noncomputable def Kp (p : ℝ) : ℝ := 4 * p / (2 - p) + CG / (2 * Real.pi)

/-- §1: A_p = K_p^{1/2} (2/r)^{1/r}. -/
noncomputable def Ap (p : ℝ) : ℝ := Kp p ^ (1 / 2 : ℝ) * (2 / rr p) ^ (1 / rr p)

/-- The standing assumptions of §7: 1 < p < 2, 0 < ε ≤ ε₀. -/
structure EpsOK (p ε : ℝ) : Prop where
  p1 : 1 < p
  p2 : p < 2
  pos : 0 < ε
  le : ε ≤ eps0 p

namespace EpsOK
variable {p ε : ℝ} (h : EpsOK p ε)
include h

lemma half : ε ≤ 1 / 2 := le_trans h.le (min_le_left _ _)
lemma e2 : ε ≤ (2 / p - 1) / 2 := le_trans h.le (le_trans (min_le_right _ _) (min_le_left _ _))
lemma e3 : ε ≤ 1 - 1 / p := le_trans h.le (le_trans (min_le_right _ _) (min_le_right _ _))
lemma p0 : 0 < p := by linarith [h.p1]
lemma ip1 : 1 / p < 1 := by rw [div_lt_one h.p0]; exact h.p1
lemma ip2 : 1 / 2 < 1 / p := one_div_lt_one_div_of_lt h.p0 h.p2
lemma r0 : 0 < rr p := by unfold rr; apply div_pos <;> linarith [h.p1, h.p2]
lemma one_div_rr : 1 / p = 1 / 2 + 1 / rr p := by
  have : p ≠ 0 := h.p0.ne'
  have : 2 - p ≠ 0 := by linarith [h.p2]
  unfold rr; field_simp; ring
/-- σ' = 1/p − ε/2 ≥ 1/2. -/
lemma s1_ge : 1 / 2 ≤ 1 / p - ε / 2 := by
  have : 2 / p = 2 * (1 / p) := by ring
  linarith [h.e2, h.ip2]
/-- σ' + ε = 1/p + ε/2 < 1. -/
lemma s1e_lt : 1 / p + ε / 2 < 1 := by linarith [h.e3, h.pos]
/-- 2σ' − 1 ≥ (2 − p)/(2p) > 0. -/
lemma two_s1 : (2 - p) / (2 * p) ≤ 2 * (1 / p - ε / 2) - 1 := by
  have h2 := h.e2
  have e : (2 - p) / (2 * p) = (2 / p - 1) / 2 := by
    have : p ≠ 0 := h.p0.ne'
    field_simp
  have : 2 / p = 2 * (1 / p) := by ring
  rw [e]; linarith
lemma two_s1_pos : 0 < (2 - p) / (2 * p) := div_pos (by linarith [h.p2]) (by linarith [h.p0])

end EpsOK

/-! ## [C8] facts used in §7 -/

lemma zeta_real {a : ℝ} (ha : 1 < a) :
    riemannZeta a = ((∑' n : ℕ, 1 / (n : ℝ) ^ a : ℝ) : ℂ) := by
  rw [zeta_eq_tsum_one_div_nat_cpow (by simpa using ha), Complex.ofReal_tsum]
  congr 1; ext n
  rw [Complex.ofReal_div, Complex.ofReal_one, Complex.ofReal_cpow (Nat.cast_nonneg n),
    Complex.ofReal_natCast]

lemma tz_anti {a b : ℝ} (ha : 1 < a) (hab : a ≤ b) :
    (∑' n : ℕ, 1 / (n : ℝ) ^ b) ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ a := by
  refine Summable.tsum_le_tsum (fun n => ?_) (Real.summable_one_div_nat_rpow.mpr (by linarith))
    (Real.summable_one_div_nat_rpow.mpr ha)
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [Real.zero_rpow (by linarith : b ≠ 0), Real.zero_rpow (by linarith : a ≠ 0)]
  · have h1 : (1:ℝ) ≤ n := by exact_mod_cast hn
    exact one_div_le_one_div_of_le (by positivity) (Real.rpow_le_rpow_of_exponent_le h1 hab)

lemma tz_ge {a : ℝ} (ha : 1 < a) : 1 / (a - 1) ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ a := by
  have hs := Real.summable_one_div_nat_rpow.mpr ha
  have key : ∀ m : ℕ, (1 - (1 + m : ℝ) ^ (1 - a)) / (a - 1) ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ a := by
    intro m
    have hanti : AntitoneOn (fun x : ℝ => x ^ (-a)) (Icc 1 (1 + m)) := by
      intro x hx y hy hxy
      exact Real.rpow_le_rpow_of_nonpos (by linarith [hx.1]) hxy (by linarith)
    have h1 := hanti.integral_le_sum
    have h2 : ∫ x in (1:ℝ)..1 + m, x ^ (-a) = (1 - (1 + m : ℝ) ^ (1 - a)) / (a - 1) := by
      rw [integral_rpow (Or.inr ⟨by intro h; linarith, ?_⟩)]
      · have e : -a + 1 = 1 - a := by ring
        rw [e, Real.one_rpow]
        have : (1 - a) ≠ 0 := by linarith
        have : (a - 1) ≠ 0 := by linarith
        field_simp
        ring
      · intro h0
        rw [Set.mem_uIcc] at h0
        have := (Nat.cast_nonneg m : (0:ℝ) ≤ m)
        rcases h0 with h0 | h0 <;> linarith [h0.1, h0.2]
    have h3 : ∑ i ∈ Finset.range m, (1 + (i : ℝ)) ^ (-a) ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ a := by
      have e : ∀ i : ℕ, (1 + (i : ℝ)) ^ (-a) = 1 / ((i + 1 : ℕ) : ℝ) ^ a := by
        intro i
        rw [Real.rpow_neg (by positivity), one_div]; push_cast; ring_nf
      simp_rw [e]
      have hs1 : Summable (fun i : ℕ => 1 / ((i + 1 : ℕ) : ℝ) ^ a) := by
        exact (summable_nat_add_iff 1).mpr hs
      calc ∑ i ∈ Finset.range m, 1 / ((i + 1 : ℕ) : ℝ) ^ a ≤ ∑' i : ℕ, 1 / ((i + 1 : ℕ) : ℝ) ^ a :=
            hs1.sum_le_tsum _ (fun i _ => by positivity)
        _ ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ a := by
            rw [hs.tsum_eq_zero_add]
            simp [Real.zero_rpow (by linarith : a ≠ 0)]
    rw [h2] at h1
    exact h1.trans h3
  have ht : Tendsto (fun m : ℕ => (1 - (1 + m : ℝ) ^ (1 - a)) / (a - 1)) atTop
      (𝓝 ((1 - 0) / (a - 1))) := by
    refine ((tendsto_const_nhds.sub ?_).div_const _)
    have h1 : Tendsto (fun x : ℝ => x ^ (-(a - 1))) atTop (𝓝 0) := tendsto_rpow_neg_atTop (by linarith)
    have h2 : Tendsto (fun m : ℕ => (1 + m : ℝ)) atTop atTop :=
      tendsto_atTop_add_const_left _ (1:ℝ) tendsto_natCast_atTop_atTop
    have := h1.comp h2
    refine this.congr (fun m => ?_)
    simp only [Function.comp_apply, neg_sub]
  rw [sub_zero] at ht
  exact le_of_tendsto' ht key

/-- [C8]: 1/ζ(1+ε) ≤ ε. -/
lemma norm_inv_zeta_le {ε : ℝ} (hε : 0 < ε) : ‖(riemannZeta (1 + ε))⁻¹‖ ≤ ε := by
  have h := tz_ge (a := 1 + ε) (by linarith)
  have e : ((1 + ε : ℝ) : ℂ) = 1 + (ε : ℂ) := by push_cast; ring
  rw [← e, zeta_real (by linarith), norm_inv, Complex.norm_real, Real.norm_eq_abs]
  have hpos : 0 < ∑' n : ℕ, 1 / (n : ℝ) ^ (1 + ε) := lt_of_lt_of_le (by
    rw [show 1 + ε - 1 = ε by ring]; exact one_div_pos.mpr hε) h
  rw [abs_of_pos hpos]
  rw [show 1 + ε - 1 = ε by ring] at h
  rw [inv_le_comm₀ hpos hε]
  simpa using h

/-! ## A1: identities on (0,1] -/

lemma chiP_of_mem {x : ℝ} (hx : x ∈ Ioc (0:ℝ) 1) : chiP x = 1 := by simp [chiP, hx.2]

/-- A1: H_ε = Φ_ε − 1 − {1/x}/ζ(1+ε) on (0,1]. -/
lemma A1b (ε : ℝ) {x : ℝ} (hx : x ∈ Ioc (0:ℝ) 1) :
    1 + HE ε x = PhiE ε x - ((Int.fract (1 / x) : ℝ) : ℂ) * (riemannZeta (1 + ε))⁻¹ := by
  have hx0 : 0 ≤ 1 / x := by have := hx.1; positivity
  have hfl : ⌊((⌊1 / x⌋₊ : ℕ) : ℝ)⌋₊ = ⌊1 / x⌋₊ := Nat.floor_natCast _
  have hfr : ((Int.fract (1 / x) : ℝ) : ℂ) = ((1 / x : ℝ) : ℂ) - (((⌊1 / x⌋₊ : ℕ) : ℝ) : ℂ) := by
    rw [Int.fract, ← Int.natCast_floor_eq_floor hx0, Int.cast_natCast, Complex.ofReal_sub]
  simp only [HE, PhiE, EE, chiP_of_mem hx, hfl, hfr]
  ring

lemma Icc_one_eq (N : ℕ) (hN : 1 ≤ N) : Finset.Icc 1 N = insert 1 (Finset.Icc 2 N) := by
  ext k; simp only [Finset.mem_Icc, Finset.mem_insert]; omega

/-- A1: H_{ε,N} = Φ_{ε,N} − 1 − P_N(1+ε){1/x} on (0,1]. -/
lemma A1c (ε : ℝ) {N : ℕ} (hN : 1 ≤ N) {x : ℝ} (hx : x ∈ Ioc (0:ℝ) 1) :
    HEN ε N x = PhiEN ε N x - 1 - PN N (1 + ε) * ((Int.fract (1 / x) : ℝ) : ℂ) := by
  have h1 : HEN ε N x = ∑ k ∈ Finset.Icc 1 N, coef ε N k * (fk k x : ℂ) := by
    rw [HEN, Icc_one_eq N hN, Finset.sum_insert (by simp)]
    simp [fk]
  have h2 : ∀ k ∈ Finset.Icc 1 N, coef ε N k * (fk k x : ℂ) = coef ε N k * (rho k x : ℂ) -
      (ArithmeticFunction.moebius k : ℂ) * (1 - (k : ℂ) / N) * (k : ℂ) ^ (-(1 + (ε : ℂ))) *
        ((Int.fract (1 / x) : ℝ) : ℂ) := by
    intro k hk
    have hk0 : (k : ℂ) ≠ 0 := by
      have := (Finset.mem_Icc.mp hk).1
      exact_mod_cast (by omega : k ≠ 0)
    have hcp : (k : ℂ) ^ (-(1 + (ε : ℂ))) = (k : ℂ) ^ (-(ε : ℂ)) * (k : ℂ)⁻¹ := by
      rw [neg_add, cpow_add _ _ hk0, cpow_neg_one, mul_comm]
    rw [hcp]
    simp only [fk, rho, coef, mu]
    push_cast
    ring
  rw [h1, PhiEN, chiP_of_mem hx, PN, Finset.sum_congr rfl h2, Finset.sum_sub_distrib,
    ← Finset.sum_mul]
  ring

lemma approx_cA (ε : ℝ) (N : ℕ) (x : ℝ) : approx N (cA ε N) x = 1 + HEN ε N x := by
  simp only [approx, cA, HEN, neg_mul, Finset.sum_neg_distrib, sub_neg_eq_add]

/-! ## measurability -/

lemma measurable_fract_inv : Measurable (fun x : ℝ => ((Int.fract (1 / x) : ℝ) : ℂ)) :=
  Complex.measurable_ofReal.comp (show Measurable (fun x : ℝ => Int.fract (1 / x)) from
    measurable_fract.comp (measurable_const.div measurable_id))

lemma measurable_HE (ε : ℝ) : Measurable (HE ε) := by
  show Measurable (fun x => (fun n : ℕ => -EE ε (n : ℝ)) ⌊1 / x⌋₊)
  exact measurable_from_nat.comp (Nat.measurable_floor.comp (measurable_const.div measurable_id))

lemma measurable_HEN (ε : ℝ) (N : ℕ) : Measurable (HEN ε N) :=
  Finset.measurable_sum _ fun k _ =>
    measurable_const.mul (Complex.measurable_ofReal.comp (measurable_fk k))

/-! ## L² on (0,1] from L²_a -/

lemma vol01_univ : vol01 Set.univ = 1 := by simp

lemma eLpNorm_weight_le {a : ℝ} {F : ℝ → ℂ} (hF : MemL2w a F) :
    eLpNorm (fun x => ((x ^ (a - 1 / 2) : ℝ) : ℂ) * F x) 2 vol01 ≤
      ENNReal.ofReal (Real.sqrt (L2w a F)) := by
  have hFm : AEStronglyMeasurable F vol01 :=
    hF.1.mono_measure (Measure.restrict_mono Ioc_subset_Ioi_self le_rfl)
  have hm : AEStronglyMeasurable (fun x => ((x ^ (a - 1 / 2) : ℝ) : ℂ) * F x) vol01 :=
    (Complex.measurable_ofReal.comp (measurable_id.pow_const _)).aestronglyMeasurable.mul hFm
  rw [eLpNorm_eq_lintegral_rpow_enorm_toReal two_ne_zero ENNReal.ofNat_ne_top hm]
  have heq : ∀ x ∈ Ioi (0:ℝ), ‖((x ^ (a - 1 / 2) : ℝ) : ℂ) * F x‖ₑ ^ (ENNReal.toReal 2) =
      ENNReal.ofReal (‖F x‖ ^ 2 * x ^ (2 * a - 1)) := by
    intro x hx
    rw [← ofReal_norm, norm_mul, Complex.norm_real,
      Real.norm_of_nonneg (Real.rpow_nonneg (le_of_lt hx) _), ENNReal.toReal_ofNat,
      ENNReal.ofReal_rpow_of_nonneg (mul_nonneg (Real.rpow_nonneg (le_of_lt hx) _) (norm_nonneg _))
        (by norm_num)]
    congr 1
    rw [Real.rpow_two, mul_pow, ← Real.rpow_natCast (x ^ (a - 1 / 2)) 2,
      ← Real.rpow_mul (le_of_lt hx)]
    have e : (a - 1 / 2) * ((2 : ℕ) : ℝ) = 2 * a - 1 := by push_cast; ring
    rw [e, mul_comm]
  have hle : ∫⁻ x in Ioc (0:ℝ) 1, ‖((x ^ (a - 1 / 2) : ℝ) : ℂ) * F x‖ₑ ^ (ENNReal.toReal 2) ≤
      ENNReal.ofReal (L2w a F) := by
    calc ∫⁻ x in Ioc (0:ℝ) 1, ‖((x ^ (a - 1 / 2) : ℝ) : ℂ) * F x‖ₑ ^ (ENNReal.toReal 2)
        ≤ ∫⁻ x in Ioi (0:ℝ), ‖((x ^ (a - 1 / 2) : ℝ) : ℂ) * F x‖ₑ ^ (ENNReal.toReal 2) :=
          lintegral_mono_set Ioc_subset_Ioi_self
      _ = ∫⁻ x in Ioi (0:ℝ), ENNReal.ofReal (‖F x‖ ^ 2 * x ^ (2 * a - 1)) :=
          setLIntegral_congr_fun measurableSet_Ioi heq
      _ = ENNReal.ofReal (L2w a F) := by
          rw [L2w, ofReal_integral_eq_lintegral_ofReal hF.2
            (ae_restrict_of_forall_mem measurableSet_Ioi fun x hx =>
              mul_nonneg (by positivity) (Real.rpow_nonneg (le_of_lt hx) _))]
  have h0 : 0 ≤ L2w a F := setIntegral_nonneg measurableSet_Ioi fun x hx =>
    mul_nonneg (by positivity) (Real.rpow_nonneg (le_of_lt hx) _)
  calc (∫⁻ x in Ioc (0:ℝ) 1, ‖((x ^ (a - 1 / 2) : ℝ) : ℂ) * F x‖ₑ ^ (ENNReal.toReal 2)) ^
        (1 / ENNReal.toReal 2) ≤ (ENNReal.ofReal (L2w a F)) ^ (1 / ENNReal.toReal 2) :=
        ENNReal.rpow_le_rpow hle (by norm_num)
    _ = ENNReal.ofReal (Real.sqrt (L2w a F)) := by
        rw [ENNReal.toReal_ofNat, Real.sqrt_eq_rpow, ENNReal.ofReal_rpow_of_nonneg h0 (by norm_num)]

/-! ## constants -/

lemma gCG_nonneg (t : ℝ) : 0 ≤ gCG t := by
  unfold gCG
  have h1 : 0 ≤ Real.log (2 + |t|) := Real.log_nonneg (by linarith [abs_nonneg t])
  have h2 : 0 ≤ (2 + |t|) ^ (1 / 2 : ℝ) := Real.rpow_nonneg (by linarith [abs_nonneg t]) _
  have h3 : 0 ≤ (1 / 4 + t ^ 2)⁻¹ := inv_nonneg.mpr (by positivity)
  exact mul_nonneg (mul_nonneg (by linarith) h2) h3

lemma CG_nonneg : 0 ≤ CG := mul_nonneg (Real.exp_pos 4).le (integral_nonneg gCG_nonneg)

lemma EpsOK.Kp_nonneg {p ε : ℝ} (h : EpsOK p ε) : 0 ≤ Kp p :=
  add_nonneg (div_nonneg (by linarith [h.p1]) (by linarith [h.p2]))
    (div_nonneg CG_nonneg (by positivity))

lemma EpsOK.Ap_nonneg {p ε : ℝ} (h : EpsOK p ε) : 0 ≤ Ap p * ε ^ (1 - 1 / p) :=
  mul_nonneg (mul_nonneg (Real.rpow_nonneg h.Kp_nonneg _)
    (Real.rpow_nonneg (div_pos two_pos h.r0).le _)) (Real.rpow_nonneg h.pos.le _)

lemma eps0_pos {p : ℝ} (hp1 : 1 < p) (hp2 : p < 2) : 0 < eps0 p := by
  have h1 : 1 < 2 / p := by rw [lt_div_iff₀ (by linarith)]; linarith
  have h2 : 1 / p < 1 := by rw [div_lt_one (by linarith)]; exact hp1
  unfold eps0
  exact lt_min (by norm_num) (lt_min (by linarith) (by linarith))

/-- (1 − r(σ' − ½))^{−1/r} √(K_p ε) = A_p ε^{1−1/p}, σ' = 1/p − ε/2 (A5's constant). -/
lemma lemmaM_const {p ε : ℝ} (h : EpsOK p ε) :
    (1 - rr p * ((1 / p - ε / 2) - 1 / 2)) ^ (-1 / rr p) * Real.sqrt (Kp p * ε) =
      Ap p * ε ^ (1 - 1 / p) := by
  have hr := h.r0
  have hε := h.pos
  have hK : 0 ≤ Kp p := h.Kp_nonneg
  have hr2 : 0 ≤ rr p / 2 := (div_pos hr two_pos).le
  have e2 : rr p * (1 / rr p) = 1 := by field_simp
  have e1 : 1 - rr p * ((1 / p - ε / 2) - 1 / 2) = rr p / 2 * ε := by
    have : 1 / p - ε / 2 - 1 / 2 = 1 / rr p - ε / 2 := by rw [h.one_div_rr]; ring
    rw [this, mul_sub, e2]; ring
  have e3 : (rr p / 2) ^ (-1 / rr p) = (2 / rr p) ^ (1 / rr p) := by
    rw [neg_div, Real.rpow_neg hr2, ← Real.inv_rpow hr2, inv_div]
  have e4 : 1 - 1 / p = -1 / rr p + 1 / 2 := by rw [h.one_div_rr]; ring
  rw [e1, Real.mul_rpow hr2 hε.le, e3, Real.sqrt_eq_rpow, Real.mul_rpow hK hε.le, Ap, e4,
    Real.rpow_add hε]
  ring

/-! ## A5 -/

lemma one_le_ofReal {p : ℝ} (hp : 1 < p) : (1 : ℝ≥0∞) ≤ ENNReal.ofReal p := by
  rw [← ENNReal.ofReal_one]; exact ENNReal.ofReal_le_ofReal hp.le

lemma eLpNorm_fract_mul_le {p C : ℝ} (z : ℂ) (hz : ‖z‖ ≤ C) :
    eLpNorm (fun x => ((Int.fract (1 / x) : ℝ) : ℂ) * z) (ENNReal.ofReal p) vol01 ≤
      ENNReal.ofReal C := by
  have := eLpNorm_le_of_ae_bound (μ := vol01) (p := ENNReal.ofReal p) (C := C)
    (measurable_fract_inv.mul_const z).aestronglyMeasurable (Eventually.of_forall fun x => by
      rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg (Int.fract_nonneg _)]
      calc Int.fract (1 / x) * ‖z‖ ≤ 1 * C :=
            mul_le_mul (Int.fract_lt_one _).le hz (norm_nonneg _) zero_le_one
        _ = C := one_mul C)
  rwa [vol01_univ, ENNReal.one_rpow, one_mul] at this

/-- §7 A5: ‖χ + H_ε‖_p ≤ A_p ε^{1−1/p} + ε, from B(ε) = ‖Φ_ε‖²_{L²_σ'} ≤ K_p ε (A4). -/
theorem A5 {p ε : ℝ} (h : EpsOK p ε) (hB : MemL2w (1 / p - ε / 2) (PhiE ε))
    (hB2 : L2w (1 / p - ε / 2) (PhiE ε) ≤ Kp p * ε) :
    eLpNorm (fun x => 1 + HE ε x) (ENNReal.ofReal p) vol01 ≤
      ENNReal.ofReal (Ap p * ε ^ (1 - 1 / p) + ε) := by
  have hPm : AEStronglyMeasurable (PhiE ε) vol01 :=
    hB.1.mono_measure (Measure.restrict_mono Ioc_subset_Ioi_self le_rfl)
  set G : ℝ → ℂ := fun x => ((Int.fract (1 / x) : ℝ) : ℂ) * (riemannZeta (1 + ε))⁻¹ with hG
  have hG_le : eLpNorm G (ENNReal.ofReal p) vol01 ≤ ENNReal.ofReal ε :=
    eLpNorm_fract_mul_le _ (norm_inv_zeta_le h.pos)
  have hM := lemmaM (p := p) (r := rr p) (a := 1 / p - ε / 2) h.p1.le h.r0 h.one_div_rr
    (by linarith [h.pos]) hPm
  have hW := eLpNorm_weight_le hB
  have hsq : Real.sqrt (L2w (1 / p - ε / 2) (PhiE ε)) ≤ Real.sqrt (Kp p * ε) :=
    Real.sqrt_le_sqrt hB2
  calc eLpNorm (fun x => 1 + HE ε x) (ENNReal.ofReal p) vol01
      = eLpNorm (PhiE ε - G) (ENNReal.ofReal p) vol01 :=
        eLpNorm_congr_ae ((ae_restrict_iff' measurableSet_Ioc).mpr
          (Eventually.of_forall fun x hx => by simp only [Pi.sub_apply, hG]; exact A1b ε hx))
    _ ≤ eLpNorm (PhiE ε) (ENNReal.ofReal p) vol01 + eLpNorm G (ENNReal.ofReal p) vol01 :=
        eLpNorm_sub_le (one_le_ofReal h.p1)
    _ ≤ ENNReal.ofReal ((1 - rr p * ((1 / p - ε / 2) - 1 / 2)) ^ (-1 / rr p)) *
          ENNReal.ofReal (Real.sqrt (Kp p * ε)) + ENNReal.ofReal ε :=
        add_le_add (hM.trans (mul_le_mul' le_rfl (hW.trans (ENNReal.ofReal_le_ofReal hsq)))) hG_le
    _ = ENNReal.ofReal (Ap p * ε ^ (1 - 1 / p) + ε) := by
        rw [← ENNReal.ofReal_mul' (Real.sqrt_nonneg _), lemmaM_const h,
          ← ENNReal.ofReal_add h.Ap_nonneg h.pos.le]

/-! ## A6 (L^p part) -/

/-- §7 A6: ‖H_ε − H_{ε,N}‖_p ≤ C' N^{−ε/2}, from ‖Φ_ε − Φ_{ε,N}‖²_{L²_σ''} ≤ C N^{−ε} and
Lemma P at s = 1 + ε. -/
theorem A6 {p ε : ℝ} (h : EpsOK p ε) {C : ℝ}
    (hD : ∀ N : ℕ, 1 ≤ N → MemL2w (1 / p - ε / 4) (PhiE ε - PhiEN ε N) ∧
      L2w (1 / p - ε / 4) (PhiE ε - PhiEN ε N) ≤ C * (N : ℝ) ^ (-ε))
    (hP : ∀ N : ℕ, 1 ≤ N → ‖PN N (1 + ε) - invZeta (1 + ε)‖ ≤
      (1 / (1 - ε / 2) + 1 / (ε / 2)) * (N : ℝ) ^ (-(ε / 2))) :
    ∃ C' : ℝ, ∀ N : ℕ, 1 ≤ N →
      eLpNorm (fun x => HEN ε N x - HE ε x) (ENNReal.ofReal p) vol01 ≤
        ENNReal.ofReal (C' * (N : ℝ) ^ (-(ε / 2))) := by
  set K2 : ℝ := (1 - rr p * ((1 / p - ε / 4) - 1 / 2)) ^ (-1 / rr p) with hK2
  set cP : ℝ := 1 / (1 - ε / 2) + 1 / (ε / 2) with hcP
  have hK2n : 0 ≤ K2 := by
    have e2 : rr p * (1 / rr p) = 1 := by field_simp [h.r0.ne']
    have : 1 - rr p * ((1 / p - ε / 4) - 1 / 2) = rr p / 4 * ε := by
      have : 1 / p - ε / 4 - 1 / 2 = 1 / rr p - ε / 4 := by rw [h.one_div_rr]; ring
      rw [this, mul_sub, e2]; ring
    rw [hK2, this]
    exact Real.rpow_nonneg (mul_nonneg (div_pos h.r0 (by norm_num)).le h.pos.le) _
  have hcPn : 0 ≤ cP := by
    have : 0 < 1 - ε / 2 := by linarith [h.half]
    have : 0 < ε / 2 := by linarith [h.pos]
    positivity
  refine ⟨K2 * Real.sqrt (max C 0) + cP, fun N hN => ?_⟩
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
  have hNr : 0 ≤ (N : ℝ) ^ (-(ε / 2)) := Real.rpow_nonneg hNpos.le _
  obtain ⟨hDm, hDb⟩ := hD N hN
  set D := PhiE ε - PhiEN ε N with hD'
  have hDm' : AEStronglyMeasurable D vol01 :=
    hDm.1.mono_measure (Measure.restrict_mono Ioc_subset_Ioi_self le_rfl)
  set z : ℂ := PN N (1 + ε) - (riemannZeta (1 + ε))⁻¹ with hz
  have hinv : invZeta (1 + ε) = (riemannZeta (1 + ε))⁻¹ := by
    have : (1 : ℂ) + ε ≠ 1 := by
      intro h0
      have := congrArg Complex.re h0
      simp at this
      linarith [h.pos]
    simp [invZeta, this]
  have hzb : ‖z‖ ≤ cP * (N : ℝ) ^ (-(ε / 2)) := by rw [hz, ← hinv]; exact hP N hN
  have hident : (fun x => HEN ε N x - HE ε x) =ᵐ[vol01]
      (-D) - (fun x => ((Int.fract (1 / x) : ℝ) : ℂ) * z) := by
    refine (ae_restrict_iff' measurableSet_Ioc).mpr (Eventually.of_forall fun x hx => ?_)
    have h1 := A1b ε hx
    have h2 := A1c ε hN hx
    simp only [Pi.sub_apply, Pi.neg_apply, hD', hz]
    rw [h2]
    linear_combination (-1 : ℂ) * h1
  have hsq : Real.sqrt (L2w (1 / p - ε / 4) D) ≤ Real.sqrt (max C 0) * (N : ℝ) ^ (-(ε / 2)) := by
    calc Real.sqrt (L2w (1 / p - ε / 4) D) ≤ Real.sqrt (max C 0 * (N : ℝ) ^ (-ε)) :=
          Real.sqrt_le_sqrt (hDb.trans (mul_le_mul_of_nonneg_right (le_max_left _ _)
            (Real.rpow_nonneg hNpos.le _)))
      _ = Real.sqrt (max C 0) * (N : ℝ) ^ (-(ε / 2)) := by
          rw [Real.sqrt_mul (le_max_right _ _), Real.sqrt_eq_rpow ((N : ℝ) ^ (-ε)),
            ← Real.rpow_mul hNpos.le]
          congr 2; ring
  have hM := lemmaM (p := p) (r := rr p) (a := 1 / p - ε / 4) h.p1.le h.r0 h.one_div_rr
    (by linarith [h.pos]) hDm'
  have hW := eLpNorm_weight_le hDm
  calc eLpNorm (fun x => HEN ε N x - HE ε x) (ENNReal.ofReal p) vol01
      = eLpNorm ((-D) - (fun x => ((Int.fract (1 / x) : ℝ) : ℂ) * z)) (ENNReal.ofReal p) vol01 :=
        eLpNorm_congr_ae hident
    _ ≤ eLpNorm (-D) (ENNReal.ofReal p) vol01 +
          eLpNorm (fun x => ((Int.fract (1 / x) : ℝ) : ℂ) * z) (ENNReal.ofReal p) vol01 :=
        eLpNorm_sub_le (one_le_ofReal h.p1)
    _ ≤ ENNReal.ofReal K2 * ENNReal.ofReal (Real.sqrt (max C 0) * (N : ℝ) ^ (-(ε / 2))) +
          ENNReal.ofReal (cP * (N : ℝ) ^ (-(ε / 2))) := by
        rw [eLpNorm_neg]
        exact add_le_add (hM.trans (mul_le_mul' le_rfl (hW.trans (ENNReal.ofReal_le_ofReal hsq))))
          (eLpNorm_fract_mul_le _ hzb)
    _ = ENNReal.ofReal ((K2 * Real.sqrt (max C 0) + cP) * (N : ℝ) ^ (-(ε / 2))) := by
        rw [← ENNReal.ofReal_mul hK2n, ← ENNReal.ofReal_add
          (mul_nonneg hK2n (mul_nonneg (Real.sqrt_nonneg _) hNr)) (mul_nonneg hcPn hNr)]
        congr 1; ring

/-! ## A7: N → ∞ for fixed ε, then ε → 0 -/

/-- §7 A7, first half: A5 + A6 give ‖χ + H_{ε,N}‖_p ≤ A_p ε^{1−1/p} + 2ε for N ≥ N₀(p,ε). -/
theorem bound_of_A5_A6 {p ε : ℝ} (h : EpsOK p ε)
    (h5 : eLpNorm (fun x => 1 + HE ε x) (ENNReal.ofReal p) vol01 ≤
      ENNReal.ofReal (Ap p * ε ^ (1 - 1 / p) + ε))
    (h6 : ∃ C' : ℝ, ∀ N : ℕ, 1 ≤ N →
      eLpNorm (fun x => HEN ε N x - HE ε x) (ENNReal.ofReal p) vol01 ≤
        ENNReal.ofReal (C' * (N : ℝ) ^ (-(ε / 2)))) :
    ∃ N₀ : ℕ, ∀ N ≥ N₀, eLpNorm (approx N (cA ε N)) (ENNReal.ofReal p) vol01 ≤
      ENNReal.ofReal (Ap p * ε ^ (1 - 1 / p) + 2 * ε) := by
  obtain ⟨C', hC'⟩ := h6
  have ht : Tendsto (fun N : ℕ => C' * (N : ℝ) ^ (-(ε / 2))) atTop (𝓝 (C' * 0)) :=
    ((tendsto_rpow_neg_atTop (by linarith [h.pos])).comp tendsto_natCast_atTop_atTop).const_mul C'
  rw [mul_zero] at ht
  obtain ⟨N₁, hN₁⟩ := eventually_atTop.mp (ht.eventually (Iic_mem_nhds h.pos))
  refine ⟨max N₁ 1, fun N hN => ?_⟩
  have hN1 : 1 ≤ N := le_trans (le_max_right _ _) hN
  have hb : C' * (N : ℝ) ^ (-(ε / 2)) ≤ ε := hN₁ N (le_trans (le_max_left _ _) hN)
  have e : approx N (cA ε N) = (fun x => 1 + HE ε x) + (fun x => HEN ε N x - HE ε x) := by
    funext x; simp only [Pi.add_apply, approx_cA]; ring
  rw [e]
  calc eLpNorm ((fun x => 1 + HE ε x) + (fun x => HEN ε N x - HE ε x)) (ENNReal.ofReal p) vol01
      ≤ eLpNorm (fun x => 1 + HE ε x) (ENNReal.ofReal p) vol01 +
          eLpNorm (fun x => HEN ε N x - HE ε x) (ENNReal.ofReal p) vol01 :=
        eLpNorm_add_le (one_le_ofReal h.p1)
    _ ≤ ENNReal.ofReal (Ap p * ε ^ (1 - 1 / p) + ε) + ENNReal.ofReal ε :=
        add_le_add h5 ((hC' N hN1).trans (ENNReal.ofReal_le_ofReal hb))
    _ = ENNReal.ofReal (Ap p * ε ^ (1 - 1 / p) + 2 * ε) := by
        rw [← ENNReal.ofReal_add (add_nonneg h.Ap_nonneg h.pos.le) h.pos.le]
        congr 1; ring

/-- §7 A7, second half: if the bound holds for every ε ∈ (0, ε₀], then d_{N,p} → 0. -/
theorem chiInW_of_eps {p A ε₀ : ℝ} (hp : 1 < p) (hε₀ : 0 < ε₀)
    (h : ∀ ε, 0 < ε → ε ≤ ε₀ → ∃ N₀ : ℕ, ∀ N ≥ N₀,
      eLpNorm (approx N (cA ε N)) (ENNReal.ofReal p) vol01 ≤
        ENNReal.ofReal (A * ε ^ (1 - 1 / p) + 2 * ε)) :
    ChiInW p := by
  intro δ hδ
  have hq : 0 < 1 - 1 / p := by
    have : 1 / p < 1 := by rw [div_lt_one (by linarith)]; exact hp
    linarith
  have ht : Tendsto (fun ε : ℝ => A * ε ^ (1 - 1 / p) + 2 * ε) (𝓝[>] 0) (𝓝 0) := by
    have h1 : Tendsto (fun ε : ℝ => ε ^ (1 - 1 / p)) (𝓝 0) (𝓝 0) := by
      have := (Real.continuousAt_rpow_const 0 (1 - 1 / p) (Or.inr hq.le)).tendsto
      rwa [Real.zero_rpow hq.ne'] at this
    have h2 : Tendsto (fun ε : ℝ => A * ε ^ (1 - 1 / p) + 2 * ε) (𝓝 0) (𝓝 (A * 0 + 2 * 0)) :=
      (h1.const_mul A).add (tendsto_id.const_mul 2)
    simpa using h2.mono_left nhdsWithin_le_nhds
  obtain ⟨ε, hε1, hε2⟩ := ((ht.eventually (Iio_mem_nhds hδ)).and (Ioc_mem_nhdsGT hε₀)).exists
  obtain ⟨N₀, hN₀⟩ := h ε hε2.1 hε2.2
  exact ⟨N₀, cA ε N₀, lt_of_le_of_lt (hN₀ N₀ le_rfl) ((ENNReal.ofReal_lt_ofReal_iff hδ).mpr hε1)⟩

/-- §7 assembled: from A4's output (B(ε) ≤ K_p ε) and A6's L²_σ'' input, for every ε ∈ (0, ε₀]. -/
theorem chiInW_of_A4_A6 {p : ℝ} (hp1 : 1 < p) (hp2 : p < 2)
    (hA4 : ∀ ε, EpsOK p ε → MemL2w (1 / p - ε / 2) (PhiE ε) ∧
      L2w (1 / p - ε / 2) (PhiE ε) ≤ Kp p * ε)
    (hA6 : ∀ ε, EpsOK p ε → ∃ C : ℝ, ∀ N : ℕ, 1 ≤ N →
      MemL2w (1 / p - ε / 4) (PhiE ε - PhiEN ε N) ∧
        L2w (1 / p - ε / 4) (PhiE ε - PhiEN ε N) ≤ C * (N : ℝ) ^ (-ε))
    (hP : ∀ ε, EpsOK p ε → ∀ N : ℕ, 1 ≤ N → ‖PN N (1 + ε) - invZeta (1 + ε)‖ ≤
      (1 / (1 - ε / 2) + 1 / (ε / 2)) * (N : ℝ) ^ (-(ε / 2))) :
    ChiInW p := by
  refine chiInW_of_eps (A := Ap p) hp1 (eps0_pos hp1 hp2) fun ε hε hεle => ?_
  have h : EpsOK p ε := ⟨hp1, hp2, hε, hεle⟩
  obtain ⟨C, hC⟩ := hA6 ε h
  exact bound_of_A5_A6 h (A5 h (hA4 ε h).1 (hA4 ε h).2) (A6 h hC (hP ε h))

end NBLp


