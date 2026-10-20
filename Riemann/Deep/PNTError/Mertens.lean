import Riemann.Deep.PNTError.General
import Riemann.Deep.PNTError.InvZeta

/-!
# Power-saving PNT, part 5: M(x) = Σ_{n ≤ x} μ(n) = O(x^{σ₀+ε}) from a zero-free half-plane

`mertens_error_of_zeroFree`: if ζ has no zero in Re s > σ₀ (1/2 ≤ σ₀ < 1), then for every ε > 0
there is C with |M(x)| ≤ C x^{σ₀+ε} for x ≥ 2.

Proof: truncated Perron for 1/ζ(s) = Σ μ(n) n^{-s} on Re s = 4 (`perron_gen`, the ψ argument with
general coefficients |aₙ| ≤ n), contour shift to Re s = σ₀ + ε/2 with no residue (1/ζ is regular at
s = 1, `shift_mu`), and 1/ζ(s) ≪ |t|^{ε/20} there (`inv_zeta_bound`, three circles), T = x⁵.
-/

open Complex Set MeasureTheory intervalIntegral Filter Topology Asymptotics
open scoped Real Interval LSeries.notation ArithmeticFunction.Moebius

namespace Fayed.PNTError

open Fayed.InputsProof.PerronAux

noncomputable section

/-- The Mertens function M(x) = Σ_{n ≤ x} μ(n). -/
def mertens (x : ℝ) : ℝ := ∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, (μ n : ℝ)

/-- A Perron integrand G(s) x^s / s. -/
def FG (G : ℂ → ℂ) (x : ℝ) (s : ℂ) : ℂ := G s * (x : ℂ) ^ s / s

lemma term_mul_eq_gen {a : ℕ → ℂ} (ha0 : a 0 = 0) {x : ℝ} (hx : 0 < x) (s : ℂ) (n : ℕ) :
    LSeries.term a s n * (x : ℂ) ^ s / s = a n * kern (Lx x n) s := by
  rcases eq_or_ne n 0 with rfl | hn
  · simp [ha0]
  · rw [LSeries.term_of_ne_zero hn, kern, Lx]
    have hn' : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn
    have hx' : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr hx.ne'
    rw [cpow_def_of_ne_zero hn', cpow_def_of_ne_zero hx', ← Complex.natCast_log,
      ← Complex.ofReal_log hx.le]
    push_cast
    rw [sub_mul, Complex.exp_sub]
    ring

lemma FG_eq_tsum {a : ℕ → ℂ} (ha0 : a 0 = 0) {G : ℂ → ℂ}
    (hG : ∀ s : ℂ, 1 < s.re → G s = LSeries a s) {x : ℝ} (hx : 0 < x) {s : ℂ} (hs : 1 < s.re) :
    FG G x s = ∑' n : ℕ, a n * kern (Lx x n) s := by
  rw [FG, hG s hs, LSeries, ← tsum_mul_right, ← tsum_div_const]
  exact tsum_congr (fun n => term_mul_eq_gen ha0 hx s n)

lemma hasSum_PG {a : ℕ → ℂ} (ha0 : a 0 = 0) (ha : ∀ n : ℕ, ‖a n‖ ≤ n) {G : ℂ → ℂ}
    (hG : ∀ s : ℂ, 1 < s.re → G s = LSeries a s) {x T : ℝ} (hx : 0 < x) (hT : 0 < T) :
    HasSum (fun n : ℕ => a n *
        (((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, kern (Lx x n) ((4 : ℝ) + y * I)))
      (((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, FG G x ((4 : ℝ) + y * I)) := by
  set c : ℂ := ((1 / (2 * π) : ℝ) : ℂ)
  set F : ℕ → ℝ → ℂ := fun n y => a n * kern (Lx x n) ((4 : ℝ) + y * I) with hF
  have hint : ∀ n, Integrable (F n) (volume.restrict (Ioc (-T) T)) := fun n =>
    ((continuous_const.mul (kern_cont' _ 4 (by norm_num))).integrableOn_Icc).mono_set
      Ioc_subset_Icc_self
  have hbd : ∀ n, ∫ y in Ioc (-T) T, ‖F n y‖ ≤ 2 * T * x ^ 4 * (1 / (n : ℝ) ^ 3) := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · simp [F, ha0]
    · have hn0 : (0 : ℝ) < n := Nat.cast_pos.mpr (Nat.pos_of_ne_zero hn)
      have hn0' : (n : ℝ) ≠ 0 := hn0.ne'
      have h1 : ∀ y, ‖F n y‖ ≤ x ^ 4 * (1 / (n : ℝ) ^ 3) := by
        intro y
        simp only [F]
        rw [norm_mul]
        have hk := norm_kern_four (Lx x n) y
        rw [exp_four_Lx hx hn] at hk
        calc ‖a n‖ * ‖kern (Lx x n) ((4 : ℝ) + y * I)‖ ≤ n * (x ^ 4 / (n : ℝ) ^ 4 / 4) :=
              mul_le_mul (ha n) hk (norm_nonneg _) (Nat.cast_nonneg n)
          _ = x ^ 4 * (1 / (n : ℝ) ^ 3) / 4 := by field_simp
          _ ≤ x ^ 4 * (1 / (n : ℝ) ^ 3) := by
              have : 0 ≤ x ^ 4 * (1 / (n : ℝ) ^ 3) := by positivity
              linarith
      have hBint : IntegrableOn (fun _ : ℝ => x ^ 4 * (1 / (n : ℝ) ^ 3)) (Ioc (-T) T) :=
        (continuous_const.integrableOn_Icc).mono_set Ioc_subset_Icc_self
      calc ∫ y in Ioc (-T) T, ‖F n y‖ ≤ ∫ y in Ioc (-T) T, x ^ 4 * (1 / (n : ℝ) ^ 3) :=
            setIntegral_mono_on (hint n).norm hBint measurableSet_Ioc (fun y _ => h1 y)
        _ = (T - -T) * (x ^ 4 * (1 / (n : ℝ) ^ 3)) := by
            rw [setIntegral_const, Real.volume_real_Ioc_of_le (by linarith), smul_eq_mul]
        _ = 2 * T * x ^ 4 * (1 / (n : ℝ) ^ 3) := by ring
  have hsum : Summable (fun n => ∫ y in Ioc (-T) T, ‖F n y‖) :=
    Summable.of_nonneg_of_le (fun n => integral_nonneg (fun y => norm_nonneg _)) hbd
      ((Real.summable_one_div_nat_pow.mpr (by norm_num : 1 < 3)).mul_left (2 * T * x ^ 4))
  have hswap := integral_tsum_of_summable_integral_norm hint hsum
  have hsumF : Summable (fun n => ∫ y in Ioc (-T) T, F n y) :=
    Summable.of_norm_bounded hsum (fun n => norm_integral_le_integral_norm _)
  have hf : ∀ n, a n * (c * ∫ y in (-T)..T, kern (Lx x n) ((4 : ℝ) + y * I)) =
      c * ∫ y in Ioc (-T) T, F n y := by
    intro n
    rw [intervalIntegral.integral_of_le (by linarith)]
    simp only [F]
    rw [MeasureTheory.integral_const_mul]
    ring
  have hFz : ∫ y in (-T)..T, FG G x ((4 : ℝ) + y * I) = ∫ y in Ioc (-T) T, ∑' n, F n y := by
    rw [intervalIntegral.integral_of_le (by linarith)]
    congr 1
    funext y
    exact FG_eq_tsum ha0 hG hx (by norm_num)
  rw [hFz, ← hswap]
  simp_rw [hf]
  exact hsumF.hasSum.mul_left c

/-- **Truncated Perron formula** for a Dirichlet series with |aₙ| ≤ n (line Re s = 4). -/
theorem perron_gen {a : ℕ → ℂ} (ha0 : a 0 = 0) (ha : ∀ n : ℕ, ‖a n‖ ≤ n) {G : ℂ → ℂ}
    (hG : ∀ s : ℂ, 1 < s.re → G s = LSeries a s) {x T : ℝ} (hx : 1 ≤ x) (hT : 0 < T)
    (hfar : ∀ n : ℕ, 1 / 2 ≤ |x - n|) :
    ‖(∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, a n) -
        ((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, FG G x ((4 : ℝ) + y * I)‖ ≤
      4 * x ^ 4 * (x + 1) / (π * T) * (π ^ 2 / 6) := by
  have hx0 : 0 < x := by linarith
  have hS : HasSum (fun n : ℕ => a n * (if 0 < Lx x n then 1 else 0))
      (∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, a n) := by
    have hfin : ∀ n ∉ Finset.Ioc 0 ⌊x⌋₊, a n * (if 0 < Lx x n then 1 else 0) = 0 := by
      intro n hn
      rcases eq_or_ne n 0 with rfl | hn0
      · simp [ha0]
      · have h1 : ⌊x⌋₊ < n := by
          simp only [Finset.mem_Ioc, not_and, not_le] at hn
          exact hn (Nat.pos_of_ne_zero hn0)
        have hxn : x < n := (Nat.floor_lt hx0.le).mp h1
        have h2 : Lx x n < 0 := by
          simp only [Lx]
          linarith [Real.log_lt_log hx0 hxn]
        simp [not_lt.mpr h2.le]
    have hh : HasSum (fun n : ℕ => a n * (if 0 < Lx x n then 1 else 0))
        (∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, a n * (if 0 < Lx x n then 1 else 0)) :=
      hasSum_sum_of_ne_finset_zero hfin
    have e : (∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, a n * (if 0 < Lx x n then 1 else 0)) =
        ∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, a n := by
      refine Finset.sum_congr rfl (fun n hn => ?_)
      have hn' := Finset.mem_Ioc.mp hn
      have hnx : (n : ℝ) < x := by
        have h1 : (n : ℝ) ≤ x := (Nat.le_floor_iff hx0.le).mp hn'.2
        rcases h1.lt_or_eq with h | h
        · exact h
        · have := hfar n
          rw [h, sub_self, abs_zero] at this
          norm_num at this
      have h3 : 0 < Lx x n := by
        simp only [Lx]
        have hn0 : (0 : ℝ) < n := Nat.cast_pos.mpr hn'.1
        linarith [Real.log_lt_log hn0 hnx]
      simp [h3]
    rw [e] at hh
    exact hh
  have hd := hS.sub (hasSum_PG ha0 ha hG hx0 hT)
  have hg : HasSum (fun n : ℕ => 4 * x ^ 4 * (x + 1) / (π * T) * (1 / (n : ℝ) ^ 2))
      (4 * x ^ 4 * (x + 1) / (π * T) * (π ^ 2 / 6)) := hasSum_zeta_two.mul_left _
  refine hd.norm_le_of_bounded hg (fun n => ?_)
  rcases eq_or_ne n 0 with rfl | hn
  · simp [ha0]
  · have hn1 : (1 : ℝ) ≤ n := Nat.one_le_cast.mpr (Nat.one_le_iff_ne_zero.mpr hn)
    have hn0 : (0 : ℝ) < n := by linarith
    have hLabs : 1 / 2 / (x + n) ≤ |Lx x n| := log_sub_ge' hx0 hn0 (by norm_num) (hfar n)
    have hLpos : 0 < |Lx x n| := lt_of_lt_of_le (by positivity) hLabs
    have hk := perron_kernel (abs_pos.mp hLpos) (by norm_num : (0 : ℝ) < 4) hT
    rw [← mul_sub, norm_mul, norm_sub_rev]
    rw [exp_four_Lx hx0 hn] at hk
    exact (mul_le_mul_of_nonneg_left hk (norm_nonneg _)).trans
      (key_ineq (ha n) hn1 hx hT hLabs)

/-! ### The Möbius case -/

/-- μ as a complex sequence. -/
abbrev aμ : ℕ → ℂ := fun n => ((μ n : ℤ) : ℂ)

lemma aμ_zero : aμ 0 = 0 := by simp [aμ]

lemma norm_aμ_le (n : ℕ) : ‖aμ n‖ ≤ n := by
  rcases eq_or_ne n 0 with rfl | hn
  · simp [aμ]
  · have hn1 : (1 : ℝ) ≤ n := Nat.one_le_cast.mpr (Nat.one_le_iff_ne_zero.mpr hn)
    have h : ‖aμ n‖ ≤ 1 := by
      simp only [aμ]
      rw [Complex.norm_intCast]
      exact_mod_cast ArithmeticFunction.abs_moebius_le_one
    linarith

lemma inv_zeta_eq_LSeries {s : ℂ} (hs : 1 < s.re) : (riemannZeta s)⁻¹ = LSeries aμ s := by
  have h := LSeries_one_mul_Lseries_moebius hs
  rw [LSeries_one_eq_riemannZeta hs] at h
  exact (eq_inv_of_mul_eq_one_right h).symm

/-- The Perron integrand for M: x^s / (s ζ(s)). -/
abbrev Fμ (x : ℝ) : ℂ → ℂ := FG (fun s => (riemannZeta s)⁻¹) x

lemma mertens_cast (x : ℝ) : (mertens x : ℂ) = ∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, aμ n := by
  simp [mertens, aμ]

theorem perron_mu {x T : ℝ} (hx : 1 ≤ x) (hT : 0 < T) (hfar : ∀ n : ℕ, 1 / 2 ≤ |x - n|) :
    ‖(mertens x : ℂ) - ((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, Fμ x ((4 : ℝ) + y * I)‖ ≤
      4 * x ^ 4 * (x + 1) / (π * T) * (π ^ 2 / 6) := by
  rw [mertens_cast]
  exact perron_gen aμ_zero norm_aμ_le (fun s hs => inv_zeta_eq_LSeries hs) hx hT hfar

lemma Fμ_differentiableAt {σ₀ : ℝ} (hZ : ZeroFree σ₀) {x : ℝ} (hx : 0 < x) {s : ℂ}
    (hs : σ₀ < s.re) (hs0 : s ≠ 0) (hs1 : s ≠ 1) : DifferentiableAt ℂ (Fμ x) s := by
  have hx' : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr hx.ne'
  have h1 : DifferentiableAt ℂ (fun s => (riemannZeta s)⁻¹) s :=
    (differentiableAt_riemannZeta hs1).inv (hZ s hs)
  have h2 : DifferentiableAt ℂ (fun s : ℂ => (x : ℂ) ^ s) s :=
    DifferentiableAt.const_cpow differentiableAt_id (Or.inl hx')
  exact (h1.mul h2).div differentiableAt_id hs0

/-- 1/ζ vanishes at s = 1, so x^s/(sζ(s)) is bounded near 1. -/
lemma Fμ_near_one {x : ℝ} (hx : 0 < x) :
    (Fμ x - fun s => (0 : ℂ) / (s - 1)) =O[𝓝[≠] 1] (1 : ℂ → ℂ) := by
  have hx' : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr hx.ne'
  have hc : ContinuousAt (fun s : ℂ => (s - 1) * (riemannZeta₁ s)⁻¹ * ((x : ℂ) ^ s / s)) 1 := by
    have h1 : ContinuousAt riemannZeta₁ 1 := differentiable_riemannZeta₁.continuous.continuousAt
    have h2 : ContinuousAt (fun s => (riemannZeta₁ s)⁻¹) 1 := h1.inv₀ (by simp)
    have h3 : ContinuousAt (fun s : ℂ => (x : ℂ) ^ s / s) 1 :=
      ((DifferentiableAt.const_cpow differentiableAt_id (Or.inl hx')).div differentiableAt_id
        one_ne_zero).continuousAt
    exact ((continuousAt_id.sub continuousAt_const).mul h2).mul h3
  have hO : (fun s : ℂ => (s - 1) * (riemannZeta₁ s)⁻¹ * ((x : ℂ) ^ s / s))
      =O[𝓝[≠] 1] (fun _ => (1 : ℂ)) :=
    (hc.tendsto.mono_left nhdsWithin_le_nhds).isBigO_one ℂ
  refine hO.congr' ?_ (Eventually.of_forall fun s => rfl)
  filter_upwards [inv_riemannZeta_eq_sub_mul] with s hs
  simp only [Pi.sub_apply, Fμ, FG, zero_div, sub_zero, hs]
  ring

/-- **Contour shift for x^s/(sζ(s))**: no residue in [σ₁, 4] × [-T, T]. -/
theorem shift_mu {σ₀ σ₁ x T : ℝ} (hZ : ZeroFree σ₀) (h01 : σ₀ < σ₁) (hσ₁0 : 0 < σ₁)
    (hσ₁ : σ₁ < 1) (hx : 0 < x) (hT : 0 < T) :
    ((1 / (2 * π) : ℝ) : ℂ) * (∫ y in (-T)..T, Fμ x ((4 : ℝ) + y * I)) -
      ((1 / (2 * π) : ℝ) : ℂ) * (∫ y in (-T)..T, Fμ x ((σ₁ : ℂ) + y * I)) =
    ((1 / (2 * π) : ℝ) : ℂ) * (I * (∫ t in σ₁..4, Fμ x (t + (-T : ℝ) * I)) -
      I * (∫ t in σ₁..4, Fμ x (t + (T : ℝ) * I))) := by
  have hres : RectangleIntegral' (Fμ x) (⟨σ₁, -T⟩ : ℂ) ⟨4, T⟩ = 0 := by
    refine ResidueTheoremOnRectangleWithSimplePole' (p := 1) (A := 0)
      (by simp only; linarith) (by simp only; linarith) ?_ ?_ (Fμ_near_one hx)
    · rw [rectangle_mem_nhds_iff]
      simp only [mem_reProdIm, one_re, one_im]
      rw [uIoo_of_le (by linarith), uIoo_of_le (by linarith)]
      exact ⟨⟨hσ₁, by norm_num⟩, ⟨by linarith, hT⟩⟩
    · intro z hz
      have hzre : z.re ∈ [[σ₁, 4]] := hz.1.1
      rw [uIcc_of_le (by linarith)] at hzre
      have hz0 : z ≠ 0 := by
        intro h
        have h0 : σ₁ ≤ z.re := hzre.1
        rw [h, zero_re] at h0
        linarith
      exact (Fμ_differentiableAt hZ hx (by linarith [hzre.1]) hz0 hz.2).differentiableWithinAt
  simp only [RectangleIntegral', RectangleIntegral, HIntegral, VIntegral, smul_eq_mul] at hres
  set Hb := ∫ t in σ₁..4, Fμ x (t + (-T : ℝ) * I)
  set Ht := ∫ t in σ₁..4, Fμ x (t + (T : ℝ) * I)
  set V4 := ∫ y in (-T)..T, Fμ x ((4 : ℝ) + y * I)
  set Vs := ∫ y in (-T)..T, Fμ x ((σ₁ : ℂ) + y * I)
  have h' : Hb - Ht + I * V4 - I * Vs = 0 := by
    rcases mul_eq_zero.mp hres with h | h
    · exfalso
      have hπ : (π : ℂ) ≠ 0 := by exact_mod_cast Real.pi_pos.ne'
      simp [hπ] at h
    · linear_combination h
  rw [← mul_sub]
  congr 1
  linear_combination (-I) * h' + (V4 - Vs) * I_sq

/-! ### Bounds -/

lemma norm_Fμ {x : ℝ} (hx : 0 < x) (s : ℂ) :
    ‖Fμ x s‖ = ‖(riemannZeta s)⁻¹‖ * x ^ s.re / ‖s‖ := by
  simp only [Fμ, FG]
  rw [norm_div, norm_mul, Complex.norm_cpow_eq_rpow_re_of_pos hx]

lemma horiz_mu {σ₁ x T C η : ℝ} (hC : ∀ s : ℂ, σ₁ ≤ s.re → 4 ≤ |s.im| →
      ‖(riemannZeta s)⁻¹‖ ≤ C * |s.im| ^ η) (hC0 : 0 ≤ C)
    (hσ₁0 : 0 ≤ σ₁) (hσ₁ : σ₁ < 4) (hx : 1 ≤ x) (hT : 4 ≤ T) {y : ℝ} (hy : |y| = T) :
    ‖∫ t in σ₁..4, Fμ x (t + y * I)‖ ≤ 4 * (C * T ^ η * x ^ 4 / T) := by
  have hx0 : 0 < x := by linarith
  have hT0 : 0 < T := by linarith
  have hTη : 0 ≤ T ^ η := Real.rpow_nonneg hT0.le η
  have hK0 : 0 ≤ C * T ^ η * x ^ 4 / T :=
    div_nonneg (mul_nonneg (mul_nonneg hC0 hTη) (by positivity)) hT0.le
  have hpt : ∀ t ∈ Ι σ₁ 4, ‖Fμ x (t + y * I)‖ ≤ C * T ^ η * x ^ 4 / T := by
    intro t ht
    rw [uIoc_of_le hσ₁.le] at ht
    obtain ⟨ht1, ht2⟩ := ht
    have hre : ((t : ℂ) + y * I).re = t := by simp
    have him : ((t : ℂ) + y * I).im = y := by simp
    rw [norm_Fμ hx0, hre]
    have h1 := hC ((t : ℂ) + y * I) (by rw [hre]; linarith) (by rw [him, hy]; exact hT)
    rw [him, hy] at h1
    have h2 : x ^ t ≤ x ^ (4 : ℕ) := by
      have := Real.rpow_le_rpow_of_exponent_le hx ht2
      rwa [show (4 : ℝ) = ((4 : ℕ) : ℝ) by norm_num, Real.rpow_natCast] at this
    have h3 : T ≤ ‖(t : ℂ) + y * I‖ := by
      rw [← hy]
      have := Complex.abs_im_le_norm ((t : ℂ) + y * I)
      rwa [him] at this
    apply div_le_div₀ (mul_nonneg (mul_nonneg hC0 hTη) (by positivity)) _ hT0 h3
    exact mul_le_mul h1 h2 (Real.rpow_nonneg hx0.le t) (mul_nonneg hC0 hTη)
  refine (intervalIntegral.norm_integral_le_of_norm_le_const hpt).trans ?_
  have hσ₁' : |(4 : ℝ) - σ₁| ≤ 4 := by
    rw [abs_le]
    constructor <;> linarith
  calc C * T ^ η * x ^ 4 / T * |4 - σ₁| ≤ C * T ^ η * x ^ 4 / T * 4 :=
        mul_le_mul_of_nonneg_left hσ₁' hK0
    _ = 4 * (C * T ^ η * x ^ 4 / T) := by ring

lemma vert_mu {σ₁ x T C B η : ℝ} (hC : ∀ s : ℂ, σ₁ ≤ s.re → 4 ≤ |s.im| →
      ‖(riemannZeta s)⁻¹‖ ≤ C * |s.im| ^ η) (hC0 : 0 ≤ C) (hη : 0 ≤ η)
    (hB : ∀ y : ℝ, |y| ≤ 4 → ‖(riemannZeta (σ₁ + y * I))⁻¹‖ ≤ B)
    (hB0 : 0 ≤ B) (hσ₁0 : 0 < σ₁) (hx : 1 ≤ x) (hT : 4 ≤ T) :
    ‖∫ y in (-T)..T, Fμ x ((σ₁ : ℂ) + y * I)‖ ≤
      x ^ σ₁ * (5 * B / σ₁ + 2 * (C * T ^ η)) * (2 * Real.log (1 + T)) := by
  have hx0 : 0 < x := by linarith
  have hT0 : 0 < T := by linarith
  have hTη : 0 ≤ T ^ η := Real.rpow_nonneg hT0.le η
  have hxs : 0 < x ^ σ₁ := Real.rpow_pos_of_pos hx0 σ₁
  have hB' : 0 ≤ 5 * B / σ₁ := div_nonneg (mul_nonneg (by norm_num) hB0) hσ₁0.le
  have hCT : 0 ≤ 2 * (C * T ^ η) := mul_nonneg (by norm_num) (mul_nonneg hC0 hTη)
  have hK1 : 5 * B / σ₁ ≤ 5 * B / σ₁ + 2 * (C * T ^ η) := by linarith
  have hK2 : 2 * (C * T ^ η) ≤ 5 * B / σ₁ + 2 * (C * T ^ η) := by linarith
  have hpt : ∀ y ∈ Ioc (-T) T, ‖Fμ x ((σ₁ : ℂ) + y * I)‖ ≤
      x ^ σ₁ * (5 * B / σ₁ + 2 * (C * T ^ η)) * (1 + |y|)⁻¹ := by
    intro y hy
    have hre : ((σ₁ : ℂ) + y * I).re = σ₁ := by simp
    have him : ((σ₁ : ℂ) + y * I).im = y := by simp
    rw [norm_Fμ hx0, hre]
    have hn1 : σ₁ ≤ ‖(σ₁ : ℂ) + y * I‖ := by
      have := Complex.abs_re_le_norm ((σ₁ : ℂ) + y * I)
      rwa [hre, abs_of_pos hσ₁0] at this
    have hn2 : |y| ≤ ‖(σ₁ : ℂ) + y * I‖ := by
      have := Complex.abs_im_le_norm ((σ₁ : ℂ) + y * I)
      rwa [him] at this
    have hy1 : 0 < 1 + |y| := by positivity
    by_cases h4 : 4 ≤ |y|
    · have h1 := hC ((σ₁ : ℂ) + y * I) (by rw [hre]) (by rw [him]; exact h4)
      rw [him] at h1
      have hyT : |y| ≤ T := abs_le.mpr ⟨by linarith [hy.1], hy.2⟩
      have hpow : |y| ^ η ≤ T ^ η := Real.rpow_le_rpow (abs_nonneg y) hyT hη
      have h1' : ‖(riemannZeta ((σ₁ : ℂ) + y * I))⁻¹‖ ≤ C * T ^ η :=
        h1.trans (mul_le_mul_of_nonneg_left hpow hC0)
      have hs2 : (1 + |y|) / 2 ≤ ‖(σ₁ : ℂ) + y * I‖ := by linarith
      calc ‖(riemannZeta ((σ₁ : ℂ) + y * I))⁻¹‖ * x ^ σ₁ / ‖(σ₁ : ℂ) + y * I‖
          ≤ C * T ^ η * x ^ σ₁ / ((1 + |y|) / 2) :=
            div_le_div₀ (mul_nonneg (mul_nonneg hC0 hTη) hxs.le)
              (mul_le_mul_of_nonneg_right h1' hxs.le) (by positivity) hs2
        _ = x ^ σ₁ * (2 * (C * T ^ η)) * (1 + |y|)⁻¹ := by
            field_simp
        _ ≤ x ^ σ₁ * (5 * B / σ₁ + 2 * (C * T ^ η)) * (1 + |y|)⁻¹ :=
            mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hK2 hxs.le)
              (inv_nonneg.mpr hy1.le)
    · rw [not_le] at h4
      have h1 := hB y h4.le
      calc ‖(riemannZeta ((σ₁ : ℂ) + y * I))⁻¹‖ * x ^ σ₁ / ‖(σ₁ : ℂ) + y * I‖
          ≤ B * x ^ σ₁ / σ₁ :=
            div_le_div₀ (mul_nonneg hB0 hxs.le) (mul_le_mul_of_nonneg_right h1 hxs.le) hσ₁0 hn1
        _ = x ^ σ₁ * (5 * B / σ₁) * (1 / 5) := by ring
        _ ≤ x ^ σ₁ * (5 * B / σ₁ + 2 * (C * T ^ η)) * (1 + |y|)⁻¹ := by
            apply mul_le_mul (mul_le_mul_of_nonneg_left hK1 hxs.le) _ (by norm_num)
              (mul_nonneg hxs.le (le_trans hB' hK1))
            rw [one_div]
            exact inv_anti₀ hy1 (by linarith)
  have hint : IntervalIntegrable
      (fun y : ℝ => x ^ σ₁ * (5 * B / σ₁ + 2 * (C * T ^ η)) * (1 + |y|)⁻¹) volume (-T) T :=
    (continuous_const.mul (Continuous.inv₀ (by fun_prop) (fun y => by positivity))).intervalIntegrable
      _ _
  refine (intervalIntegral.norm_integral_le_of_norm_le (by linarith) (ae_of_all _ hpt) hint).trans
    (le_of_eq ?_)
  rw [intervalIntegral.integral_const_mul, integral_inv_one_add_abs hT0.le]

/-- 1/ζ is bounded on the segment Re s = σ₁, |Im s| ≤ 4. -/
lemma vert_compact_mu {σ₀ σ₁ : ℝ} (hZ : ZeroFree σ₀) (h01 : σ₀ < σ₁) (hσ₁ : σ₁ < 1) :
    ∃ B : ℝ, 0 ≤ B ∧ ∀ y : ℝ, |y| ≤ 4 → ‖(riemannZeta (σ₁ + y * I))⁻¹‖ ≤ B := by
  have hK : IsCompact ((fun y : ℝ => (σ₁ : ℂ) + y * I) '' Icc (-4) 4) :=
    isCompact_Icc.image (by fun_prop)
  have hc : ContinuousOn (fun s => (riemannZeta s)⁻¹) ((fun y : ℝ => (σ₁ : ℂ) + y * I) '' Icc (-4) 4) := by
    rintro _ ⟨y, -, rfl⟩
    have h1 : σ₀ < ((σ₁ : ℂ) + y * I).re := by simpa using h01
    have h2 : (σ₁ : ℂ) + y * I ≠ 1 := by
      intro h
      have := congrArg Complex.re h
      simp at this
      linarith
    exact ((differentiableAt_riemannZeta h2).continuousAt.inv₀ (hZ _ h1)).continuousWithinAt
  obtain ⟨B, hB⟩ := hK.exists_bound_of_continuousOn hc
  refine ⟨max B 0, le_max_right _ _, fun y hy => ?_⟩
  exact (hB _ ⟨y, abs_le.mp hy, rfl⟩).trans (le_max_left _ _)

end

end Fayed.PNTError
