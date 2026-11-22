module

public import Riemann.Deep.LeanLi.Thm4c

/-!
# LeanLi: Montgomery's pair correlation ⇒ `Σ_{γ,γ'≤T} (1 + log T·|γ - γ'|)^{-2} ≪ T log T`

Fejér kernel `k(y) = 2∫_0^{1/2} (1/2 - α) cos(αy) dα = 4 sin²(y/4)/y² ≥ 0` with Fourier support `[-1/2, 1/2]`:
`Σ_{γ,γ'} w(γ - γ') k(L(γ - γ') - c) = 2∫_0^{1/2} (1/2 - α) Re(e^{-iαc} (TL/2π) F(α)) dα ≪ TL`, and
`(1 + |x|)^{-2} ≤ π²(k(x) + k(x - 2π))`.
-/

@[expose] public section

open Real Filter Set Finset MeasureTheory Topology Complex

namespace LeanLi

/-- The Fejér kernel. -/
noncomputable def kern (y : ℝ) : ℝ := 2 * ∫ α in (0 : ℝ)..(1 / 2), (1 / 2 - α) * Real.cos (α * y)

lemma hasDerivAt_kprim {y : ℝ} (hy : y ≠ 0) (α : ℝ) :
    HasDerivAt (fun α : ℝ => (1 / 2 - α) * Real.sin (α * y) / y - Real.cos (α * y) / y ^ 2)
      ((1 / 2 - α) * Real.cos (α * y)) α := by
  have hs : HasDerivAt (fun α : ℝ => Real.sin (α * y)) (Real.cos (α * y) * y) α := (hasDerivAt_mul_const y).sin
  have hc : HasDerivAt (fun α : ℝ => Real.cos (α * y)) (-Real.sin (α * y) * y) α := (hasDerivAt_mul_const y).cos
  have hl : HasDerivAt (fun α : ℝ => 1 / 2 - α) (-1) α := (hasDerivAt_id' α).const_sub (1 / 2)
  have h := ((hl.mul hs).div_const y).sub (hc.div_const (y ^ 2))
  convert h using 1
  field_simp
  ring

lemma kern_eq {y : ℝ} (hy : y ≠ 0) : kern y = 4 * Real.sin (y / 4) ^ 2 / y ^ 2 := by
  unfold kern
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun α _ => hasDerivAt_kprim hy α)
    ((by fun_prop : Continuous fun α : ℝ => (1 / 2 - α) * Real.cos (α * y)).intervalIntegrable _ _)]
  have h2 : Real.cos (1 / 2 * y) = 1 - 2 * Real.sin (y / 4) ^ 2 := by
    rw [show 1 / 2 * y = 2 * (y / 4) by ring, Real.cos_two_mul, Real.cos_sq']; ring
  simp only [sub_self, zero_mul, Real.sin_zero, Real.cos_zero, mul_zero, zero_div, zero_sub]
  rw [h2]
  field_simp
  ring

lemma kern_zero : kern 0 = 1 / 4 := by
  unfold kern
  simp only [mul_zero, Real.cos_zero, mul_one]
  rw [intervalIntegral.integral_sub intervalIntegrable_const intervalIntegral.intervalIntegrable_id,
    intervalIntegral.integral_const, integral_id]
  norm_num

lemma kern_nonneg (y : ℝ) : 0 ≤ kern y := by
  rcases eq_or_ne y 0 with h | h
  · rw [h, kern_zero]; norm_num
  · rw [kern_eq h]; positivity

lemma kern_ge {y D : ℝ} (hD : 0 < D) (hyD : y ^ 2 ≤ D) : 4 * Real.sin (y / 4) ^ 2 / D ≤ kern y := by
  rcases eq_or_ne y 0 with h | h
  · rw [h]; simp only [zero_div, Real.sin_zero]; norm_num; exact kern_nonneg 0
  · rw [kern_eq h]
    exact div_le_div_of_nonneg_left (by positivity) (by positivity) hyD

/-- `(1 + |x|)^{-2} ≤ π² (k(x) + k(x - 2π))`. -/
lemma inv_sq_le_kern (x : ℝ) : 1 / (1 + |x|) ^ 2 ≤ Real.pi ^ 2 * (kern x + kern (x - 2 * Real.pi)) := by
  have hπ := Real.pi_gt_three
  set D := (|x| + 2 * Real.pi) ^ 2 with hD
  have hD0 : 0 < D := by positivity
  have h1 : x ^ 2 ≤ D := by
    rw [hD, ← sq_abs x]; exact pow_le_pow_left₀ (abs_nonneg _) (by linarith) 2
  have h2 : (x - 2 * Real.pi) ^ 2 ≤ D := by
    rw [hD, ← sq_abs (x - 2 * Real.pi)]
    apply pow_le_pow_left₀ (abs_nonneg _)
    have := abs_sub x (2 * Real.pi)
    rw [abs_of_pos (by positivity : (0 : ℝ) < 2 * Real.pi)] at this
    exact this
  have k1 := kern_ge hD0 h1
  have k2 := kern_ge hD0 h2
  have hsc : Real.sin ((x - 2 * Real.pi) / 4) ^ 2 = Real.cos (x / 4) ^ 2 := by
    rw [show (x - 2 * Real.pi) / 4 = x / 4 - Real.pi / 2 by ring, Real.sin_sub_pi_div_two, neg_sq]
  rw [hsc] at k2
  have hsum : 4 / D ≤ kern x + kern (x - 2 * Real.pi) := by
    have e : 4 / D = 4 * Real.sin (x / 4) ^ 2 / D + 4 * Real.cos (x / 4) ^ 2 / D := by
      rw [← add_div, ← mul_add, Real.sin_sq_add_cos_sq, mul_one]
    linarith
  have h3 : 1 / (1 + |x|) ^ 2 ≤ Real.pi ^ 2 * (4 / D) := by
    rw [mul_div_assoc', div_le_div_iff₀ (by positivity) hD0, hD]
    have : |x| + 2 * Real.pi ≤ 2 * Real.pi * (1 + |x|) := by nlinarith [abs_nonneg x]
    have := pow_le_pow_left₀ (by positivity) this 2
    nlinarith
  exact h3.trans (mul_le_mul_of_nonneg_left hsum (by positivity))

/-- Pair weight `w(u) = 4/(4 + u²)`. -/
noncomputable def wPC (u : ℝ) : ℝ := 4 / (4 + u ^ 2)

lemma wPC_pos (u : ℝ) : 0 < wPC u := by unfold wPC; positivity

/-- Double sums over indexed zeros as weighted double sums over zeros. -/
lemma sum_ZK_ZK (T : ℝ) (f : ℝ → ℝ → ℝ) :
    ∑ j ∈ ZK T, ∑ k ∈ ZK T, f (ordK j) (ordK k) =
      ∑ ρ ∈ zset T, ∑ ρ' ∈ zset T, (mult ρ : ℝ) * mult ρ' * f ρ.im ρ'.im := by
  have h1 : ∀ j : ZeroIdx, ∑ k ∈ ZK T, f (ordK j) (ordK k) = ∑ ρ' ∈ zset T, (mult ρ' : ℝ) * f (ordK j) ρ'.im :=
    fun j => sum_ZK T (f (ordK j))
  simp only [h1]
  have h2 := sum_ZK T (fun x => ∑ ρ' ∈ zset T, (mult ρ' : ℝ) * f x ρ'.im)
  rw [h2]
  refine Finset.sum_congr rfl fun ρ _ => ?_
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun ρ' _ => by ring

/-- `P(c) = Σ_{j,k} w(Δ) k(LΔ - c)`. -/
noncomputable def Ppair (T c : ℝ) : ℝ :=
  ∑ j ∈ ZK T, ∑ k ∈ ZK T, wPC (ordK j - ordK k) * kern (Real.log T * (ordK j - ordK k) - c)

/-- `G(α) = Σ_{j,k} w(Δ) cos(α(LΔ - c))`. -/
noncomputable def Gpair (T c α : ℝ) : ℝ :=
  ∑ j ∈ ZK T, ∑ k ∈ ZK T, wPC (ordK j - ordK k) * Real.cos (α * (Real.log T * (ordK j - ordK k) - c))

lemma continuous_Gpair (T c : ℝ) : Continuous (Gpair T c) := by
  unfold Gpair
  exact continuous_finset_sum _ fun j _ => continuous_finset_sum _ fun k _ => by fun_prop

lemma Ppair_eq (T c : ℝ) : Ppair T c = 2 * ∫ α in (0 : ℝ)..(1 / 2), (1 / 2 - α) * Gpair T c α := by
  have hterm : ∀ a y : ℝ, a * kern y = 2 * ∫ α in (0 : ℝ)..(1 / 2), (1 / 2 - α) * (a * Real.cos (α * y)) := by
    intro a y
    unfold kern
    rw [show a * (2 * ∫ α in (0 : ℝ)..(1 / 2), (1 / 2 - α) * Real.cos (α * y)) =
      2 * (a * ∫ α in (0 : ℝ)..(1 / 2), (1 / 2 - α) * Real.cos (α * y)) by ring,
      ← intervalIntegral.integral_const_mul]
    congr 1
    exact intervalIntegral.integral_congr fun α _ => by ring
  unfold Ppair Gpair
  simp only [hterm, ← Finset.mul_sum]
  congr 1
  have e1 : ∀ α : ℝ, (1 / 2 - α) * ∑ j ∈ ZK T, ∑ k ∈ ZK T,
      wPC (ordK j - ordK k) * Real.cos (α * (Real.log T * (ordK j - ordK k) - c)) =
      ∑ j ∈ ZK T, ∑ k ∈ ZK T, (1 / 2 - α) *
        (wPC (ordK j - ordK k) * Real.cos (α * (Real.log T * (ordK j - ordK k) - c))) := fun α => by
    simp only [Finset.mul_sum]
  simp only [e1]
  rw [intervalIntegral.integral_finset_sum (fun j _ => (continuous_finset_sum _ fun k _ =>
    (by fun_prop : Continuous fun α : ℝ => (1 / 2 - α) *
      (wPC (ordK j - ordK k) * Real.cos (α * (Real.log T * (ordK j - ordK k) - c))))).intervalIntegrable _ _)]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [intervalIntegral.integral_finset_sum (fun k _ => (by fun_prop : Continuous fun α : ℝ => (1 / 2 - α) *
      (wPC (ordK j - ordK k) * Real.cos (α * (Real.log T * (ordK j - ordK k) - c)))).intervalIntegrable _ _)]

lemma finsum_zerosUpTo {M : Type*} [AddCommMonoid M] (T : ℝ) (f : ℂ → M) :
    ∑ᶠ ρ ∈ zerosUpTo T, f ρ = ∑ ρ ∈ zset T, f ρ :=
  finsum_mem_eq_finite_toFinset_sum _ (finite_zerosUpTo T)

lemma Gpair_le (T c α : ℝ) (hT : 1 < T) :
    Gpair T c α ≤ T * Real.log T / (2 * Real.pi) * ‖montgomeryF α T‖ := by
  have hπ := Real.pi_pos
  have hL : 0 < Real.log T := Real.log_pos hT
  have hT0 : 0 < T := by linarith
  set Z : ℂ := ∑ ρ ∈ zset T, ∑ ρ' ∈ zset T, ((mult ρ * mult ρ' : ℕ) : ℂ) *
    Complex.exp (I * (α * (ρ.im - ρ'.im) * Real.log T)) * ((4 / (4 + (ρ.im - ρ'.im) ^ 2) : ℝ) : ℂ) with hZ
  have hF : montgomeryF α T = ((2 * Real.pi / (T * Real.log T) : ℝ) : ℂ) * Z := by
    unfold montgomeryF
    simp only [finsum_zerosUpTo]
    rw [hZ]
  have hG : Gpair T c α = (Complex.exp (-(I * (α * c))) * Z).re := by
    unfold Gpair
    rw [sum_ZK_ZK T (fun x y => wPC (x - y) * Real.cos (α * (Real.log T * (x - y) - c))), hZ, Finset.mul_sum,
      Complex.re_sum]
    refine Finset.sum_congr rfl fun ρ _ => ?_
    rw [Finset.mul_sum, Complex.re_sum]
    refine Finset.sum_congr rfl fun ρ' _ => ?_
    have h1 : Complex.exp (-(I * (α * c))) * Complex.exp (I * (α * (ρ.im - ρ'.im) * Real.log T)) =
        Complex.exp (((α * (Real.log T * (ρ.im - ρ'.im) - c)) : ℝ) * I) := by
      rw [← Complex.exp_add]; congr 1; push_cast; ring
    have e : Complex.exp (-(I * (α * c))) * (((mult ρ * mult ρ' : ℕ) : ℂ) *
        Complex.exp (I * (α * (ρ.im - ρ'.im) * Real.log T)) * ((4 / (4 + (ρ.im - ρ'.im) ^ 2) : ℝ) : ℂ)) =
        (((mult ρ : ℝ) * mult ρ' * wPC (ρ.im - ρ'.im) : ℝ) : ℂ) *
          (Complex.exp (-(I * (α * c))) * Complex.exp (I * (α * (ρ.im - ρ'.im) * Real.log T))) := by
      unfold wPC; push_cast; ring
    rw [e, h1, Complex.re_ofReal_mul, Complex.exp_ofReal_mul_I_re]
    ring
  rw [hG]
  have hn : ‖Complex.exp (-(I * (α * c)))‖ = 1 := by
    rw [show -(I * (α * c : ℂ)) = ((-(α * c) : ℝ) : ℂ) * I by push_cast; ring, Complex.norm_exp_ofReal_mul_I]
  calc (Complex.exp (-(I * (α * c))) * Z).re ≤ ‖Complex.exp (-(I * (α * c))) * Z‖ := Complex.re_le_norm _
    _ = ‖Z‖ := by rw [norm_mul, hn, one_mul]
    _ = T * Real.log T / (2 * Real.pi) * ‖montgomeryF α T‖ := by
        rw [hF, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
        field_simp

/-- **Montgomery ⇒ `P(c) ≤ T log T/2`** uniformly in `c`. -/
lemma Ppair_le (L : LiteratureInputs) (hRH : RiemannHypothesis) :
    ∃ T₀ : ℝ, 3 ≤ T₀ ∧ ∀ T, T₀ ≤ T → ∀ c : ℝ, Ppair T c ≤ T * Real.log T / 2 := by
  obtain ⟨T₁, hT₁⟩ := L.montgomery hRH 1 one_pos
  refine ⟨max T₁ 3, le_max_right _ _, fun T hT c => ?_⟩
  have hT3 : 3 ≤ T := le_trans (le_max_right _ _) hT
  have hT0 : 0 < T := by linarith
  have hπ := Real.pi_gt_three
  have hL : 0 < Real.log T := Real.log_pos (by linarith)
  set ℓ := Real.log T with hℓ
  set A := T * ℓ / (2 * Real.pi) with hA
  have hA0 : 0 ≤ A := by positivity
  set bnd : ℝ → ℝ := fun α => 1 / 2 * (A * (α + 2 * (ℓ * Real.exp (-(2 * ℓ) * α)) + 1)) with hbnd
  have hpt : ∀ α ∈ Set.Icc (0 : ℝ) (1 / 2), (1 / 2 - α) * Gpair T c α ≤ bnd α := by
    intro α hα
    have hG := Gpair_le T c α (by linarith)
    have hM := hT₁ T (le_trans (le_max_left _ _) hT) α hα.1 hα.2
    have hrp : T ^ (-2 * α) = Real.exp (-(2 * ℓ) * α) := by
      rw [Real.rpow_def_of_pos hT0, hℓ]; ring_nf
    rw [hrp] at hM
    have hE0 : 0 ≤ Real.exp (-(2 * ℓ) * α) * ℓ := by positivity
    have h2 : ‖(α : ℂ) + ((Real.exp (-(2 * ℓ) * α) * ℓ : ℝ) : ℂ)‖ = α + Real.exp (-(2 * ℓ) * α) * ℓ := by
      rw [← Complex.ofReal_add, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (by linarith [hα.1])]
    have h3 := norm_sub_norm_le (montgomeryF α T) ((α : ℂ) + ((Real.exp (-(2 * ℓ) * α) * ℓ : ℝ) : ℂ))
    have hF : ‖montgomeryF α T‖ ≤ α + 2 * (ℓ * Real.exp (-(2 * ℓ) * α)) + 1 := by
      rw [h2] at h3; nlinarith
    have hGb : Gpair T c α ≤ A * (α + 2 * (ℓ * Real.exp (-(2 * ℓ) * α)) + 1) :=
      hG.trans (mul_le_mul_of_nonneg_left hF hA0)
    have hB0 : 0 ≤ A * (α + 2 * (ℓ * Real.exp (-(2 * ℓ) * α)) + 1) := by
      have : 0 ≤ ℓ * Real.exp (-(2 * ℓ) * α) := by positivity
      have : 0 ≤ α := hα.1
      positivity
    calc (1 / 2 - α) * Gpair T c α ≤ (1 / 2 - α) * (A * (α + 2 * (ℓ * Real.exp (-(2 * ℓ) * α)) + 1)) :=
          mul_le_mul_of_nonneg_left hGb (by linarith [hα.2])
      _ ≤ bnd α := by
          simp only [hbnd]
          exact mul_le_mul_of_nonneg_right (by linarith [hα.1]) hB0
  have hcb : Continuous bnd := by simp only [hbnd]; fun_prop
  have hint1 : IntervalIntegrable (fun α => (1 / 2 - α) * Gpair T c α) volume 0 (1 / 2) :=
    ((continuous_const.sub continuous_id).mul (continuous_Gpair T c)).intervalIntegrable _ _
  have hmono := intervalIntegral.integral_mono_on (by norm_num) hint1 (hcb.intervalIntegrable _ _) hpt
  have hH : ∀ α ∈ Set.uIcc (0 : ℝ) (1 / 2), HasDerivAt
      (fun α => 1 / 2 * (A * (α ^ 2 / 2 - Real.exp (-(2 * ℓ) * α) + α))) (bnd α) α := by
    intro α _
    have h1 : HasDerivAt (fun α : ℝ => Real.exp (-(2 * ℓ) * α)) (Real.exp (-(2 * ℓ) * α) * (-(2 * ℓ) * 1)) α :=
      ((hasDerivAt_id' α).const_mul (-(2 * ℓ))).exp
    have h2 : HasDerivAt (fun α : ℝ => α ^ 2 / 2) α α := by
      have := (hasDerivAt_pow 2 α).div_const 2
      convert this using 1; push_cast; ring
    have h := (((h2.sub h1).add (hasDerivAt_id' α)).const_mul A).const_mul (1 / 2 : ℝ)
    simp only [hbnd]
    convert h using 1
    ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hH (hcb.intervalIntegrable _ _)] at hmono
  rw [Ppair_eq]
  have hexp : 0 ≤ Real.exp (-(2 * ℓ) * (1 / 2)) := Real.exp_pos _ |>.le
  simp only [mul_zero, Real.exp_zero] at hmono
  generalize Real.exp (-(2 * ℓ) * (1 / 2)) = e at hmono hexp
  norm_num at hmono
  have h2A : 2 * A ≤ T * ℓ / 2 := by
    have e : 2 * A = T * ℓ / Real.pi := by rw [hA]; field_simp
    rw [e]
    exact div_le_div_of_nonneg_left (by positivity) (by norm_num) (by linarith)
  have hAe := mul_nonneg hA0 hexp
  linarith

section
variable (L : LiteratureInputs)
include L

/-- **Montgomery ⇒ few close pairs**: `Σ_{γ,γ'≤T} (1 + log T·|γ - γ'|)^{-2} ≤ K T log T`. -/
theorem pair_sum_bound (hRH : RiemannHypothesis) : ∃ K T₀ : ℝ, 0 ≤ K ∧ 3 ≤ T₀ ∧ ∀ T, T₀ ≤ T →
    ∑ j ∈ ZK T, ∑ k ∈ ZK T, 1 / (1 + Real.log T * |ordK j - ordK k|) ^ 2 ≤ K * (T * Real.log T) := by
  obtain ⟨T₁, hT₁3, hP⟩ := Ppair_le L hRH
  obtain ⟨C, T₂, hT₂, hC⟩ := local_count L
  obtain ⟨T₃, hT₃, hc⟩ := count_le L
  refine ⟨5 / 4 * Real.pi ^ 2 + 24 * |C|, max T₁ (max T₂ T₃), by positivity,
    le_trans hT₁3 (le_max_left _ _), fun T hT => ?_⟩
  have hT1 : T₁ ≤ T := le_trans (le_max_left _ _) hT
  have hT2 : T₂ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hT
  have hT3' : T₃ ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) hT
  have hT3 : 3 ≤ T := le_trans hT₁3 hT1
  have hπ := Real.pi_gt_three
  have hℓ1 : 1 ≤ Real.log T := by
    rw [Real.le_log_iff_exp_le (by linarith)]; have := Real.exp_one_lt_d9; linarith
  set ℓ := Real.log T with hℓ
  have hℓ0 : 0 < ℓ := by linarith
  have hpt : ∀ j k : ZeroIdx, 1 / (1 + ℓ * |ordK j - ordK k|) ^ 2 ≤
      5 / 4 * Real.pi ^ 2 * (wPC (ordK j - ordK k) * kern (ℓ * (ordK j - ordK k) - 0) +
        wPC (ordK j - ordK k) * kern (ℓ * (ordK j - ordK k) - 2 * Real.pi)) +
      4 / ℓ ^ 2 * (1 / (1 + |ordK k - ordK j|)) := by
    intro j k
    set Δ := ordK j - ordK k with hΔ
    have hk1 := kern_nonneg (ℓ * Δ)
    have hk2 := kern_nonneg (ℓ * Δ - 2 * Real.pi)
    have hw := wPC_pos Δ
    have hsym : |ordK k - ordK j| = |Δ| := abs_sub_comm _ _
    rw [hsym, sub_zero]
    have hfar0 : 0 ≤ 4 / ℓ ^ 2 * (1 / (1 + |Δ|)) := by positivity
    have hnear0 : 0 ≤ 5 / 4 * Real.pi ^ 2 * (wPC Δ * kern (ℓ * Δ) + wPC Δ * kern (ℓ * Δ - 2 * Real.pi)) := by
      positivity
    rcases lt_or_ge |Δ| 1 with h | h
    · have h1 := inv_sq_le_kern (ℓ * Δ)
      rw [abs_mul, abs_of_pos hℓ0] at h1
      have hw45 : 4 / 5 ≤ wPC Δ := by
        unfold wPC; rw [div_le_div_iff₀ (by norm_num) (by positivity)]
        have : Δ ^ 2 ≤ 1 := by rw [← sq_abs]; nlinarith [abs_nonneg Δ]
        nlinarith
      have h3 : kern (ℓ * Δ) + kern (ℓ * Δ - 2 * Real.pi) ≤
          5 / 4 * (wPC Δ * kern (ℓ * Δ) + wPC Δ * kern (ℓ * Δ - 2 * Real.pi)) := by
        have e1 : kern (ℓ * Δ) ≤ 5 / 4 * (wPC Δ * kern (ℓ * Δ)) := by nlinarith
        have e2 : kern (ℓ * Δ - 2 * Real.pi) ≤ 5 / 4 * (wPC Δ * kern (ℓ * Δ - 2 * Real.pi)) := by nlinarith
        linarith
      have h2 := mul_le_mul_of_nonneg_left h3 (sq_nonneg Real.pi)
      have e : Real.pi ^ 2 * (5 / 4 * (wPC Δ * kern (ℓ * Δ) + wPC Δ * kern (ℓ * Δ - 2 * Real.pi))) =
          5 / 4 * Real.pi ^ 2 * (wPC Δ * kern (ℓ * Δ) + wPC Δ * kern (ℓ * Δ - 2 * Real.pi)) := by ring
      linarith
    · have hD0 : 0 < |Δ| := by linarith
      have h1 : 1 / (1 + ℓ * |Δ|) ^ 2 ≤ 1 / (ℓ * |Δ|) ^ 2 :=
        one_div_le_one_div_of_le (by positivity) (pow_le_pow_left₀ (by positivity) (by linarith) 2)
      have h2 : 1 / (ℓ * |Δ|) ^ 2 ≤ 4 / ℓ ^ 2 * (1 / (1 + |Δ|)) := by
        rw [div_mul_div_comm, mul_one, div_le_div_iff₀ (by positivity) (by positivity)]
        have h4 : 1 + |Δ| ≤ 4 * |Δ| ^ 2 := by nlinarith
        have := mul_le_mul_of_nonneg_left h4 (sq_nonneg ℓ)
        nlinarith
      linarith
  have hB : ∀ u : ℝ, (((ZK T).filter (fun k => u ≤ ordK k ∧ ordK k ≤ u + 1)).card : ℝ) ≤ |C| * ℓ := by
    intro u
    rw [card_ZK_window]
    exact (hC T hT2 u).trans (mul_le_mul_of_nonneg_right (le_abs_self C) hℓ0.le)
  set M := ⌊T⌋₊ with hM
  have hMT : (M : ℝ) ≤ T := Nat.floor_le (by linarith)
  have hTM : T < M + 1 := Nat.lt_floor_add_one T
  have hrow : ∀ j ∈ ZK T, ∑ k ∈ ZK T, 1 / (1 + |ordK k - ordK j|) ≤ 2 * (|C| * ℓ) * (1 + Real.log (T + 1)) := by
    intro j hj
    have hfar := sum_inv_one_add_le (ZK T) ordK (ordK j) M (fun k hk => by
      have h1 := mem_ZK.1 hk; have h2 := mem_ZK.1 hj
      have h3 : 0 < ordK k := k.1.2.2; have h4 : 0 < ordK j := j.1.2.2
      rw [abs_lt]; constructor <;> linarith) (|C| * ℓ) hB
    exact hfar.trans (mul_le_mul_of_nonneg_left (harmonic_sum_le M hMT)
      (mul_nonneg (by norm_num) (mul_nonneg (abs_nonneg _) hℓ0.le)))
  have hfarsum : ∑ j ∈ ZK T, ∑ k ∈ ZK T, 1 / (1 + |ordK k - ordK j|) ≤
      (zeroCount T : ℝ) * (2 * (|C| * ℓ) * (1 + Real.log (T + 1))) := by
    refine (Finset.sum_le_sum hrow).trans (le_of_eq ?_)
    rw [Finset.sum_const, nsmul_eq_mul, card_ZK]
  have hNT : (zeroCount T : ℝ) ≤ T * ℓ := hc T hT3'
  have hH : 1 + Real.log (T + 1) ≤ 3 * ℓ := by
    have h1 : Real.log (T + 1) ≤ Real.log (2 * T) := Real.log_le_log (by linarith) (by linarith)
    rw [Real.log_mul (by norm_num) (by linarith)] at h1
    have h2 : Real.log 2 ≤ 1 := by
      have := Real.log_le_sub_one_of_pos (show (0 : ℝ) < 2 by norm_num); linarith
    linarith
  have hfar2 : 4 / ℓ ^ 2 * ∑ j ∈ ZK T, ∑ k ∈ ZK T, 1 / (1 + |ordK k - ordK j|) ≤ 24 * |C| * (T * ℓ) := by
    have hCl : 0 ≤ 2 * (|C| * ℓ) := mul_nonneg (by norm_num) (mul_nonneg (abs_nonneg _) hℓ0.le)
    have hH0 : 0 ≤ 1 + Real.log (T + 1) := by have := Real.log_nonneg (show (1 : ℝ) ≤ T + 1 by linarith); linarith
    have hTl : 0 ≤ T * ℓ := mul_nonneg (by linarith) hℓ0.le
    have h1 : (zeroCount T : ℝ) * (2 * (|C| * ℓ) * (1 + Real.log (T + 1))) ≤ (T * ℓ) * (2 * (|C| * ℓ) * (3 * ℓ)) :=
      mul_le_mul hNT (mul_le_mul_of_nonneg_left hH hCl) (mul_nonneg hCl hH0) hTl
    have h2 := mul_le_mul_of_nonneg_left (hfarsum.trans h1) (show (0 : ℝ) ≤ 4 / ℓ ^ 2 by positivity)
    have e : 4 / ℓ ^ 2 * ((T * ℓ) * (2 * (|C| * ℓ) * (3 * ℓ))) = 24 * |C| * (T * ℓ) := by
      field_simp; ring
    linarith
  have hP0 := hP T hT1 0
  have hP2 := hP T hT1 (2 * Real.pi)
  calc ∑ j ∈ ZK T, ∑ k ∈ ZK T, 1 / (1 + ℓ * |ordK j - ordK k|) ^ 2 ≤
        ∑ j ∈ ZK T, ∑ k ∈ ZK T, (5 / 4 * Real.pi ^ 2 * (wPC (ordK j - ordK k) * kern (ℓ * (ordK j - ordK k) - 0) +
          wPC (ordK j - ordK k) * kern (ℓ * (ordK j - ordK k) - 2 * Real.pi)) +
          4 / ℓ ^ 2 * (1 / (1 + |ordK k - ordK j|))) :=
        Finset.sum_le_sum fun j _ => Finset.sum_le_sum fun k _ => hpt j k
    _ = 5 / 4 * Real.pi ^ 2 * (Ppair T 0 + Ppair T (2 * Real.pi)) +
        4 / ℓ ^ 2 * ∑ j ∈ ZK T, ∑ k ∈ ZK T, 1 / (1 + |ordK k - ordK j|) := by
        unfold Ppair
        simp only [← hℓ]
        simp only [Finset.sum_add_distrib, ← Finset.mul_sum, mul_add]
    _ ≤ 5 / 4 * Real.pi ^ 2 * (T * ℓ / 2 + T * ℓ / 2) + 24 * |C| * (T * ℓ) := by
        have := mul_le_mul_of_nonneg_left (add_le_add hP0 hP2) (show (0 : ℝ) ≤ 5 / 4 * Real.pi ^ 2 by positivity)
        linarith
    _ = (5 / 4 * Real.pi ^ 2 + 24 * |C|) * (T * ℓ) := by ring

end

end LeanLi
