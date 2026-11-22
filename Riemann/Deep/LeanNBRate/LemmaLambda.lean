import Riemann.Deep.LeanNBRate.Bernstein
import Riemann.Deep.LeanNBRate.FourierDeriv

/-!
# [N] Lemma Λ: frequency decay of a windowed group component

For `r` with nodes in `[−w, w]` and the time window `ψ = twin ℓ`,
`F(ζ) = ∫ r(u) ψ(u) e^{−iζu} du = Λ_r(Θ_ζ)`, `Θ_ζ(θ) = ∫ ψ(u) e^{−i(ζ−θ)u} du`; the multiplier lemma
with `G = Θ_ζ · cutoff(−w, w, η)` gives decay `(K/(ℓ d))^N` in `d = dist(ζ, [−w, w]) ≥ 2η`.
-/

open MeasureTheory Filter Topology Set Complex FourierTransform

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The time window as a complex function. -/
def psiC (ℓ : ℝ) (u : ℝ) : ℂ := (twin ℓ u : ℂ)

theorem twin_hasCompactSupport {ℓ : ℝ} (hℓ : 0 < ℓ) : HasCompactSupport (twin ℓ) :=
  cutoff_hasCompactSupport (by positivity)

theorem psiC_isTest {ℓ : ℝ} (hℓ : 0 < ℓ) : IsTest (psiC ℓ) :=
  ⟨ofRealCLM.contDiff.comp (twin_gevB hℓ).smooth,
    (twin_hasCompactSupport hℓ).comp_left ofReal_zero⟩

/-- The modulated window `u ↦ ψ(u) e^{−iζu}`. -/
def modPsi (ℓ ζ : ℝ) (u : ℝ) : ℂ := psiC ℓ u * cexp (-(I * ζ * u))

theorem modPsi_isTest {ℓ : ℝ} (hℓ : 0 < ℓ) (ζ : ℝ) : IsTest (modPsi ℓ ζ) := by
  refine ⟨(psiC_isTest hℓ).smooth.mul ?_, (psiC_isTest hℓ).supp.mul_right⟩
  exact (contDiff_const.mul ofRealCLM.contDiff).neg.cexp

/-- `Λ_r(G)` at `x = 0` is `Σ a_{jk} G^{(k)}(δ_j)`. -/
theorem ExpSum.Lam_zero (r : ExpSum) (G : ℝ → ℂ) :
    r.Lam 0 G = ∑ j, ∑ k ∈ Finset.range (r.mult j), r.coef j k * iteratedDeriv k G (r.node j) := by
  unfold ExpSum.Lam
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
  rw [Finset.sum_range_succ]
  have hz : ∀ l ∈ Finset.range k, (Nat.choose k l : ℂ) * (I * ((0 : ℝ) : ℂ)) ^ (k - l) *
      iteratedDeriv l G (r.node j) = 0 := by
    intro l hl
    rw [Finset.mem_range] at hl
    rw [ofReal_zero, mul_zero, zero_pow (by omega), mul_zero, zero_mul]
  rw [Finset.sum_eq_zero hz]
  simp

/-- `F(ζ) = ∫ r ψ e^{−iζu} = Λ_r(Θ_ζ)`, `Θ_ζ = fint (modPsi ℓ ζ)`. -/
theorem ExpSum.integral_mul_modPsi (r : ExpSum) {ℓ : ℝ} (hℓ : 0 < ℓ) (ζ : ℝ) :
    ∫ u : ℝ, r.eval u * modPsi ℓ ζ u = r.Lam 0 (fint (modPsi ℓ ζ)) := by
  have hT := modPsi_isTest hℓ ζ
  rw [r.Lam_zero]
  simp_rw [hT.iteratedDeriv_fint]
  unfold fint pmul ExpSum.eval
  have hint : ∀ (j : Fin r.J) (k : ℕ), Integrable (fun u : ℝ =>
      r.coef j k * ((I * u) ^ k * modPsi ℓ ζ u * cexp (I * r.node j * u))) := by
    intro j k
    refine Integrable.const_mul ?_ _
    have hb := (hT.pmul_test k).integrable
    refine (hb.mul_bdd (g := fun u : ℝ => cexp (I * r.node j * u)) (c := 1) ?_ ?_).congr
      (Eventually.of_forall fun u => ?_)
    · exact (by fun_prop : Continuous fun u : ℝ => cexp (I * r.node j * u)).aestronglyMeasurable
    · refine Eventually.of_forall fun u => ?_
      rw [show (I * r.node j * u : ℂ) = ((r.node j * u : ℝ) : ℂ) * I by push_cast; ring,
        Complex.norm_exp_ofReal_mul_I]
    · simp only [pmul]
  have hexp : ∀ u : ℝ, (∑ j, ∑ k ∈ Finset.range (r.mult j),
      r.coef j k * (I * u) ^ k * cexp (I * r.node j * u)) * modPsi ℓ ζ u =
      ∑ j, ∑ k ∈ Finset.range (r.mult j),
        r.coef j k * ((I * u) ^ k * modPsi ℓ ζ u * cexp (I * r.node j * u)) := by
    intro u
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun k _ => ?_
    ring
  simp_rw [hexp]
  rw [integral_finset_sum _ fun j _ => integrable_finset_sum _ fun k _ => hint j k]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [integral_finset_sum _ fun k _ => hint j k]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [integral_const_mul]

/-- `Θ_ζ^{(a)}(θ) = 2π · kern ((iu)^a ψ) (ζ − θ)`. -/
theorem iteratedDeriv_fint_modPsi {ℓ : ℝ} (hℓ : 0 < ℓ) (ζ : ℝ) (a : ℕ) (θ : ℝ) :
    iteratedDeriv a (fint (modPsi ℓ ζ)) θ =
      (2 * Real.pi : ℂ) * kern (pmul a (psiC ℓ)) (ζ - θ) := by
  rw [(modPsi_isTest hℓ ζ).iteratedDeriv_fint, kern_eq]
  have hπ : (2 * Real.pi : ℂ) ≠ 0 := by exact_mod_cast (show (0:ℝ) < 2 * Real.pi by positivity).ne'
  rw [← mul_assoc, mul_inv_cancel₀ hπ, one_mul]
  unfold fint pmul modPsi
  refine integral_congr_ae (Eventually.of_forall fun u => ?_)
  simp only
  have : cexp (-(I * ζ * u)) * cexp (I * θ * u) = cexp (-(I * ((ζ - θ : ℝ) : ℂ) * u)) := by
    rw [← Complex.exp_add]; push_cast; ring_nf
  calc (I * u) ^ a * (psiC ℓ u * cexp (-(I * ζ * u))) * cexp (I * θ * u)
      = (I * u) ^ a * psiC ℓ u * (cexp (-(I * ζ * u)) * cexp (I * θ * u)) := by ring
    _ = cexp (-(I * ((ζ - θ : ℝ) : ℂ) * u)) * ((I * u) ^ a * psiC ℓ u) := by rw [this]; ring
    _ = cexp (-(I * u * ((ζ - θ : ℝ) : ℂ))) * ((I * u) ^ a * psiC ℓ u) := by ring_nf


theorem iteratedDeriv_twin_eq_zero {ℓ u : ℝ} (hℓ : 0 < ℓ) (hu : u ∉ Icc 0 ℓ) (n : ℕ) :
    iteratedDeriv n (twin ℓ) u = 0 := by
  unfold twin
  refine iteratedDeriv_cutoff_eq_zero (by positivity) ?_ n
  rwa [show ℓ / 2 - ℓ / 2 = (0 : ℝ) by ring, show ℓ / 2 + ℓ / 2 = ℓ by ring]

/-- Jets of `u ↦ (iu)^a`. -/
theorem iteratedDeriv_Ipow (a i : ℕ) (u : ℝ) :
    iteratedDeriv i (fun u : ℝ => (I * (u : ℂ)) ^ a) u =
      I ^ a * ((a.descFactorial i : ℝ) * u ^ (a - i) : ℝ) := by
  have e : (fun u : ℝ => (I * (u : ℂ)) ^ a) = fun u : ℝ => I ^ a * ((u ^ a : ℝ) : ℂ) := by
    ext u; push_cast; ring
  rw [e, iteratedDeriv_const_mul_field]
  have h := iteratedDeriv_ofReal_comp (f := fun u : ℝ => u ^ a) (contDiff_id.pow a) i
  rw [show (fun u : ℝ => ((u ^ a : ℝ) : ℂ)) = fun u : ℝ => (((fun u : ℝ => u ^ a) u : ℝ) : ℂ) from rfl,
    h]
  simp

theorem norm_iteratedDeriv_Ipow_le {a i : ℕ} {u ℓ : ℝ} (hℓ : 0 < ℓ) (hu : u ∈ Icc 0 ℓ) :
    ‖iteratedDeriv i (fun u : ℝ => (I * (u : ℂ)) ^ a) u‖ ≤ (a : ℝ) ^ i * ℓ ^ a / ℓ ^ i := by
  rw [iteratedDeriv_Ipow, norm_mul, norm_pow, Complex.norm_I, one_pow, one_mul, Complex.norm_real,
    Real.norm_eq_abs, abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ a.descFactorial i),
    abs_of_nonneg (pow_nonneg hu.1 _)]
  rcases le_or_gt i a with hia | hia
  · have h1 : (a.descFactorial i : ℝ) ≤ (a : ℝ) ^ i := by exact_mod_cast Nat.descFactorial_le_pow a i
    have h2 : u ^ (a - i) ≤ ℓ ^ (a - i) := pow_le_pow_left₀ hu.1 hu.2 _
    have h3 : ℓ ^ (a - i) = ℓ ^ a / ℓ ^ i := by
      rw [eq_div_iff (pow_ne_zero _ hℓ.ne'), ← pow_add, Nat.sub_add_cancel hia]
    calc (a.descFactorial i : ℝ) * u ^ (a - i) ≤ (a : ℝ) ^ i * ℓ ^ (a - i) :=
          mul_le_mul h1 h2 (pow_nonneg hu.1 _) (by positivity)
      _ = (a : ℝ) ^ i * ℓ ^ a / ℓ ^ i := by rw [h3]; ring
  · rw [Nat.descFactorial_eq_zero_iff_lt.mpr hia]; simp only [Nat.cast_zero, zero_mul]; positivity

theorem natpow_self_le {m n : ℕ} (h : m ≤ n) : (m : ℝ) ^ (2 * m) ≤ (n : ℝ) ^ (2 * n) := by
  rcases Nat.eq_zero_or_pos m with rfl | hm
  · simp only [Nat.cast_zero, mul_zero, pow_zero]
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp
    · exact one_le_pow₀ (by exact_mod_cast hn)
  · calc (m : ℝ) ^ (2 * m) ≤ (n : ℝ) ^ (2 * m) := by gcongr
      _ ≤ (n : ℝ) ^ (2 * n) := by
          apply pow_le_pow_right₀ (by exact_mod_cast (hm.trans_le h)) (by omega)

/-- `C_ψ = B₀²`. Pointwise Leibniz bound for `∂^N ((iu)^a ψ)`. -/
theorem norm_iteratedDeriv_pmul_psiC_le {ℓ : ℝ} (hℓ : 0 < ℓ) (a N : ℕ) (u : ℝ) :
    ‖iteratedDeriv N (pmul a (psiC ℓ)) u‖ ≤
      (B₀ * B₀) * ℓ ^ a / ℓ ^ N * ((a : ℝ) + 8 * A₀) ^ N * (N : ℝ) ^ (2 * N) := by
  have hψ := twin_gevB hℓ
  have hsm : ContDiff ℝ (⊤ : ℕ∞) (fun u : ℝ => (I * (u : ℂ)) ^ a) :=
    (contDiff_const.mul ofRealCLM.contDiff).pow a
  have hpsi : ContDiff ℝ (⊤ : ℕ∞) (psiC ℓ) := (psiC_isTest hℓ).smooth
  have e : pmul a (psiC ℓ) = (fun u : ℝ => (I * (u : ℂ)) ^ a) * psiC ℓ := rfl
  rw [e, iteratedDeriv_mul (hsm.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt
    (hpsi.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt]
  have hpsiD : ∀ k, iteratedDeriv k (psiC ℓ) = fun u => ((iteratedDeriv k (twin ℓ) u : ℝ) : ℂ) :=
    fun k => iteratedDeriv_ofReal_comp hψ.smooth k
  simp_rw [hpsiD]
  have hA0 := A₀_pos
  by_cases hu : u ∈ Icc 0 ℓ
  · calc ‖∑ i ∈ Finset.range (N + 1), (N.choose i : ℂ) *
          iteratedDeriv i (fun u : ℝ => (I * (u : ℂ)) ^ a) u *
            ((iteratedDeriv (N - i) (twin ℓ) u : ℝ) : ℂ)‖
        ≤ ∑ i ∈ Finset.range (N + 1), (N.choose i : ℝ) * ((a : ℝ) ^ i * ℓ ^ a / ℓ ^ i) *
            ((B₀ * B₀) * (2 * (2 / (ℓ / 2) * A₀)) ^ (N - i) * (N : ℝ) ^ (2 * N)) := by
          refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun i hi => ?_)
          rw [Finset.mem_range] at hi
          rw [norm_mul, norm_mul, Complex.norm_natCast, Complex.norm_real, Real.norm_eq_abs]
          gcongr
          · exact norm_iteratedDeriv_Ipow_le hℓ hu
          · refine (hψ.bound (N - i) u).trans ?_
            have hBB : 0 ≤ B₀ * B₀ * (2 * (2 / (ℓ / 2) * A₀)) ^ (N - i) := by
              have := B₀_pos; positivity
            exact mul_le_mul_of_nonneg_left (natpow_self_le (Nat.sub_le N i)) hBB
      _ = (B₀ * B₀) * ℓ ^ a / ℓ ^ N * (N : ℝ) ^ (2 * N) *
            ∑ i ∈ Finset.range (N + 1), (N.choose i : ℝ) * ((a : ℝ) ^ i * (8 * A₀) ^ (N - i)) := by
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun i hi => ?_
          rw [Finset.mem_range] at hi
          have h1 : (2 * (2 / (ℓ / 2) * A₀)) ^ (N - i) = (8 * A₀) ^ (N - i) / ℓ ^ (N - i) := by
            rw [← div_pow]; congr 1; field_simp; ring
          have h2 : ℓ ^ i * ℓ ^ (N - i) = ℓ ^ N := by rw [← pow_add]; congr 1; omega
          rw [h1]
          field_simp
          rw [← h2]; ring
      _ = (B₀ * B₀) * ℓ ^ a / ℓ ^ N * ((a : ℝ) + 8 * A₀) ^ N * (N : ℝ) ^ (2 * N) := by
          have hsum : ∑ i ∈ Finset.range (N + 1), (N.choose i : ℝ) * ((a : ℝ) ^ i * (8 * A₀) ^ (N - i)) =
              ((a : ℝ) + 8 * A₀) ^ N := by
            rw [add_pow]
            refine Finset.sum_congr rfl fun i _ => ?_
            ring
          rw [hsum]; ring
  · have hz : ∀ i ∈ Finset.range (N + 1), (N.choose i : ℂ) *
        iteratedDeriv i (fun u : ℝ => (I * (u : ℂ)) ^ a) u *
          ((iteratedDeriv (N - i) (twin ℓ) u : ℝ) : ℂ) = 0 := by
      intro i _
      rw [iteratedDeriv_twin_eq_zero hℓ hu]; simp
    rw [Finset.sum_eq_zero hz, norm_zero]
    have := B₀_pos
    positivity


theorem integral_norm_iteratedDeriv_pmul_psiC_le {ℓ : ℝ} (hℓ : 0 < ℓ) (a N : ℕ) :
    ∫ u : ℝ, ‖iteratedDeriv N (pmul a (psiC ℓ)) u‖ ≤
      ℓ * ((B₀ * B₀) * ℓ ^ a / ℓ ^ N * ((a : ℝ) + 8 * A₀) ^ N * (N : ℝ) ^ (2 * N)) := by
  have hT := (psiC_isTest hℓ).pmul_test a
  have hzero : ∀ u ∉ Icc 0 ℓ, iteratedDeriv N (pmul a (psiC ℓ)) u = 0 := by
    intro u hu
    have hψ := twin_gevB hℓ
    have hsm : ContDiff ℝ (⊤ : ℕ∞) (fun u : ℝ => (I * (u : ℂ)) ^ a) :=
      (contDiff_const.mul ofRealCLM.contDiff).pow a
    have e : pmul a (psiC ℓ) = (fun u : ℝ => (I * (u : ℂ)) ^ a) * psiC ℓ := rfl
    rw [e, iteratedDeriv_mul (hsm.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt
      ((psiC_isTest hℓ).smooth.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt]
    refine Finset.sum_eq_zero fun i _ => ?_
    have hpd : iteratedDeriv (N - i) (psiC ℓ) u = ((iteratedDeriv (N - i) (twin ℓ) u : ℝ) : ℂ) :=
      congrFun (iteratedDeriv_ofReal_comp hψ.smooth (N - i)) u
    rw [hpd, iteratedDeriv_twin_eq_zero hℓ hu]; simp
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero (s := Icc 0 ℓ)
    (fun u hu => by rw [hzero u hu, norm_zero])]
  calc ∫ u in Icc 0 ℓ, ‖iteratedDeriv N (pmul a (psiC ℓ)) u‖
      ≤ ∫ u in Icc 0 ℓ, (B₀ * B₀) * ℓ ^ a / ℓ ^ N * ((a : ℝ) + 8 * A₀) ^ N * (N : ℝ) ^ (2 * N) :=
        setIntegral_mono_on ((hT.iterDeriv N).continuous.norm.continuousOn.integrableOn_compact
          isCompact_Icc) (integrableOn_const (by rw [Real.volume_Icc]; exact ENNReal.ofReal_ne_top))
          measurableSet_Icc fun u _ => norm_iteratedDeriv_pmul_psiC_le hℓ a N u
    _ = _ := by
        rw [setIntegral_const, Real.volume_real_Icc_of_le hℓ.le, smul_eq_mul, sub_zero]

/-- Gevrey decay of the jets of `Θ_ζ`: `|Θ_ζ^{(a)}(θ)| ≤ B₀² ℓ^{a+1} ((a+8A₀)N²/(ℓ|ζ−θ|))^N`. -/
theorem norm_iteratedDeriv_fint_modPsi_le {ℓ : ℝ} (hℓ : 0 < ℓ) (ζ θ : ℝ) (hζθ : ζ ≠ θ) (a N : ℕ) :
    ‖iteratedDeriv a (fint (modPsi ℓ ζ)) θ‖ ≤
      (B₀ * B₀) * ℓ ^ (a + 1) * (((a : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * |ζ - θ|)) ^ N := by
  rw [iteratedDeriv_fint_modPsi hℓ, norm_mul]
  have hT := (psiC_isTest hℓ).pmul_test a
  have hk := hT.norm_pow_mul_kern_le N (ζ - θ)
  have hint := integral_norm_iteratedDeriv_pmul_psiC_le hℓ a N
  have hd : 0 < |ζ - θ| := abs_pos.mpr (sub_ne_zero.mpr hζθ)
  have h2π : ‖(2 * Real.pi : ℂ)‖ = 2 * Real.pi := by
    rw [show (2 * Real.pi : ℂ) = ((2 * Real.pi : ℝ) : ℂ) by push_cast; ring, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos (by positivity)]
  rw [h2π]
  have hkern : ‖kern (pmul a (psiC ℓ)) (ζ - θ)‖ ≤
      (2 * Real.pi)⁻¹ * (∫ u : ℝ, ‖iteratedDeriv N (pmul a (psiC ℓ)) u‖) / |ζ - θ| ^ N := by
    rw [le_div_iff₀ (pow_pos hd N), mul_comm]; exact hk
  calc 2 * Real.pi * ‖kern (pmul a (psiC ℓ)) (ζ - θ)‖
      ≤ 2 * Real.pi * ((2 * Real.pi)⁻¹ * (∫ u : ℝ, ‖iteratedDeriv N (pmul a (psiC ℓ)) u‖) /
          |ζ - θ| ^ N) := by gcongr
    _ = (∫ u : ℝ, ‖iteratedDeriv N (pmul a (psiC ℓ)) u‖) / |ζ - θ| ^ N := by
        field_simp
    _ ≤ ℓ * ((B₀ * B₀) * ℓ ^ a / ℓ ^ N * ((a : ℝ) + 8 * A₀) ^ N * (N : ℝ) ^ (2 * N)) /
          |ζ - θ| ^ N := by gcongr
    _ = (B₀ * B₀) * ℓ ^ (a + 1) * (((a : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * |ζ - θ|)) ^ N := by
        rw [div_pow, mul_pow, mul_pow, ← pow_mul, pow_succ]
        field_simp


theorem fint_eq_fourier (g : ℝ → ℂ) : fint g = fun θ => 𝓕 g ((-(2 * Real.pi)⁻¹) * θ) := by
  have hπc : (2 * Real.pi : ℂ) ≠ 0 := by
    exact_mod_cast (show (0:ℝ) < 2 * Real.pi by positivity).ne'
  ext θ
  rw [fint_eq_kern]; unfold kern
  rw [← mul_assoc, mul_inv_cancel₀ hπc, one_mul]
  congr 1; field_simp

theorem IsTest.contDiff_fint {f : ℝ → ℂ} (hf : IsTest f) : ContDiff ℝ (⊤ : ℕ∞) (fint f) := by
  rw [fint_eq_fourier]
  have : 𝓕 f = ((𝓕 hf.toS : SchwartzMap ℝ ℂ) : ℝ → ℂ) := by rw [SchwartzMap.fourier_coe]; rfl
  rw [this]
  exact ((𝓕 hf.toS).smooth ⊤).comp (contDiff_const.mul contDiff_id)

/-- The Lemma Λ multiplier `G = Θ_ζ · cutoff(−w, w, η)`. -/
def lamG (ℓ ζ w η : ℝ) (θ : ℝ) : ℂ := fint (modPsi ℓ ζ) θ * (cutoff (-w) w η θ : ℂ)

theorem lamG_isTest {ℓ ζ w η : ℝ} (hℓ : 0 < ℓ) (hη : 0 < η) : IsTest (lamG ℓ ζ w η) := by
  refine ⟨(modPsi_isTest hℓ ζ).contDiff_fint.mul
    (ofRealCLM.contDiff.comp (cutoff_gevB (-w) w η hη).smooth), ?_⟩
  have : HasCompactSupport (fun θ : ℝ => (cutoff (-w) w η θ : ℂ)) :=
    (cutoff_hasCompactSupport hη).comp_left ofReal_zero
  exact this.mul_left

theorem ExpSum.Lam_zero_lamG (r : ExpSum) {ℓ ζ w η : ℝ} (hη : 0 < η)
    (hr : ∀ j, r.node j ∈ Icc (-w) w) :
    r.Lam 0 (lamG ℓ ζ w η) = r.Lam 0 (fint (modPsi ℓ ζ)) := by
  rw [r.Lam_zero, r.Lam_zero]
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
  congr 1
  have hev : lamG ℓ ζ w η =ᶠ[𝓝 (r.node j)] fint (modPsi ℓ ζ) := by
    filter_upwards [cutoff_eventuallyEq_one hη (hr j)] with θ hθ
    simp only [lamG, hθ, ofReal_one, mul_one]
  exact hev.iteratedDeriv_eq k

/-- `∫ |G^{(n)}| ≤ 2(w+η) B₀⁴ ℓ E_n (ℓ + A_χ n²)^n`, `E_n = (2(n+8A₀)N²/(ℓd))^N`, if `|ζ| ≥ w + d`, `d ≥ 2η`. -/
theorem integral_norm_iteratedDeriv_lamG_le {ℓ ζ w η d : ℝ} (hℓ : 0 < ℓ) (hw : 0 ≤ w) (hη : 0 < η)
    (hd : 2 * η ≤ d) (hζ : w + d ≤ |ζ|) (n N : ℕ) :
    ∫ θ : ℝ, ‖iteratedDeriv n (lamG ℓ ζ w η) θ‖ ≤
      2 * (w + η) * ((B₀ * B₀) * (B₀ * B₀) * ℓ *
        (2 * ((n : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * d)) ^ N *
        (ℓ + 2 * (2 / η * A₀) * (n : ℝ) ^ 2) ^ n) := by
  set s := w + η
  set Aχ := 2 * (2 / η * A₀)
  have hs : 0 < s := by positivity
  have hA0 := A₀_pos
  have hB0 := B₀_pos
  have hAχ : 0 < Aχ := by positivity
  have hd0 : 0 < d := by linarith
  have hgev := cutoff_gevB (-w) w η hη
  have hΘ := (modPsi_isTest hℓ ζ).contDiff_fint
  have hχc : ContDiff ℝ (⊤ : ℕ∞) (fun θ : ℝ => (cutoff (-w) w η θ : ℂ)) :=
    ofRealCLM.contDiff.comp hgev.smooth
  have hLeib : ∀ θ, iteratedDeriv n (lamG ℓ ζ w η) θ = ∑ b ∈ Finset.range (n + 1),
      (n.choose b : ℂ) * iteratedDeriv b (fint (modPsi ℓ ζ)) θ *
        ((iteratedDeriv (n - b) (cutoff (-w) w η) θ : ℝ) : ℂ) := by
    intro θ
    have e : lamG ℓ ζ w η = fint (modPsi ℓ ζ) * fun θ : ℝ => (cutoff (-w) w η θ : ℂ) := rfl
    rw [e, iteratedDeriv_mul (hΘ.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt
      (hχc.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt]
    refine Finset.sum_congr rfl fun b _ => ?_
    rw [iteratedDeriv_ofReal_comp hgev.smooth]
  set E := (2 * ((n : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * d)) ^ N
  have hpt : ∀ θ ∈ Icc (-s) s, ‖iteratedDeriv n (lamG ℓ ζ w η) θ‖ ≤
      (B₀ * B₀) * (B₀ * B₀) * ℓ * E * (ℓ + Aχ * (n : ℝ) ^ 2) ^ n := by
    intro θ hθ
    have hθabs : |θ| ≤ s := abs_le.mpr hθ
    have hζθ : d / 2 ≤ |ζ - θ| := by
      have : |ζ| - |θ| ≤ |ζ - θ| := abs_sub_abs_le_abs_sub ζ θ
      linarith
    have hζθ0 : ζ ≠ θ := by
      intro h; rw [h, sub_self, abs_zero] at hζθ; linarith
    rw [hLeib θ]
    calc ‖∑ b ∈ Finset.range (n + 1), (n.choose b : ℂ) * iteratedDeriv b (fint (modPsi ℓ ζ)) θ *
          ((iteratedDeriv (n - b) (cutoff (-w) w η) θ : ℝ) : ℂ)‖
        ≤ ∑ b ∈ Finset.range (n + 1), (n.choose b : ℝ) * ((B₀ * B₀) * ℓ ^ (b + 1) * E) *
            ((B₀ * B₀) * (Aχ * (n : ℝ) ^ 2) ^ (n - b)) := by
          refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun b hb => ?_)
          rw [Finset.mem_range] at hb
          rw [norm_mul, norm_mul, Complex.norm_natCast, Complex.norm_real, Real.norm_eq_abs]
          gcongr
          · refine (norm_iteratedDeriv_fint_modPsi_le hℓ ζ θ hζθ0 b N).trans ?_
            have hX : ((b : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * |ζ - θ|) ≤
                2 * ((n : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * d) := by
              calc ((b : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * |ζ - θ|)
                  ≤ ((n : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * (d / 2)) := by
                    gcongr
                    exact_mod_cast Nat.le_of_lt_succ hb
                _ = 2 * ((n : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * d) := by field_simp
            have hXN := pow_le_pow_left₀ (by positivity) hX N
            exact mul_le_mul_of_nonneg_left hXN (by positivity)
          · refine (hgev.bound (n - b) θ).trans ?_
            have hp : ((n - b : ℕ) : ℝ) ^ (2 * (n - b)) ≤ ((n : ℝ) ^ 2) ^ (n - b) := by
              calc ((n - b : ℕ) : ℝ) ^ (2 * (n - b)) = ((((n - b : ℕ) : ℝ)) ^ 2) ^ (n - b) := by
                    rw [← pow_mul, mul_comm]
                _ ≤ ((n : ℝ) ^ 2) ^ (n - b) := by
                    gcongr; exact_mod_cast Nat.sub_le n b
            calc B₀ * B₀ * (2 * (2 / η * A₀)) ^ (n - b) * ((n - b : ℕ) : ℝ) ^ (2 * (n - b))
                ≤ B₀ * B₀ * (2 * (2 / η * A₀)) ^ (n - b) * ((n : ℝ) ^ 2) ^ (n - b) := by gcongr
              _ = (B₀ * B₀) * (Aχ * (n : ℝ) ^ 2) ^ (n - b) := by rw [mul_assoc, ← mul_pow]
      _ = (B₀ * B₀) * (B₀ * B₀) * ℓ * E *
            ∑ b ∈ Finset.range (n + 1), (n.choose b : ℝ) * (ℓ ^ b * (Aχ * (n : ℝ) ^ 2) ^ (n - b)) := by
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun b _ => ?_
          rw [pow_succ]; ring
      _ = (B₀ * B₀) * (B₀ * B₀) * ℓ * E * (ℓ + Aχ * (n : ℝ) ^ 2) ^ n := by
          congr 1
          rw [add_pow]
          refine Finset.sum_congr rfl fun b _ => ?_
          ring
  have hzero : ∀ θ ∉ Icc (-s) s, iteratedDeriv n (lamG ℓ ζ w η) θ = 0 := by
    intro θ hθ
    have hθ' : θ ∉ Icc (-w - η) (w + η) := by
      simpa [s, show -w - η = -(w + η) by ring] using hθ
    rw [hLeib θ]
    refine Finset.sum_eq_zero fun b _ => ?_
    rw [iteratedDeriv_cutoff_eq_zero hη hθ']; simp
  have hcont : Continuous fun θ => ‖iteratedDeriv n (lamG ℓ ζ w η) θ‖ :=
    ((lamG_isTest hℓ hη).iterDeriv n).continuous.norm
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero (s := Icc (-s) s)
    (fun θ hθ => by rw [hzero θ hθ, norm_zero])]
  calc ∫ θ in Icc (-s) s, ‖iteratedDeriv n (lamG ℓ ζ w η) θ‖
      ≤ ∫ θ in Icc (-s) s, (B₀ * B₀) * (B₀ * B₀) * ℓ * E * (ℓ + Aχ * (n : ℝ) ^ 2) ^ n :=
        setIntegral_mono_on (hcont.continuousOn.integrableOn_compact isCompact_Icc)
          (integrableOn_const (by rw [Real.volume_Icc]; exact ENNReal.ofReal_ne_top))
          measurableSet_Icc hpt
    _ = 2 * s * ((B₀ * B₀) * (B₀ * B₀) * ℓ * E * (ℓ + Aχ * (n : ℝ) ^ 2) ^ n) := by
        rw [setIntegral_const, Real.volume_real_Icc_of_le (by linarith), smul_eq_mul]; ring


theorem one_add_pow_le_exp {x : ℝ} (hx : 0 ≤ x) (M : ℕ) : (1 + x) ^ M ≤ Real.exp (M * x) := by
  calc (1 + x) ^ M ≤ (Real.exp x) ^ M := by
        gcongr; linarith [Real.add_one_le_exp x]
    _ = Real.exp (M * x) := by rw [← Real.exp_nat_mul]

/-- **[N] Lemma Λ**: for `|ζ| ≥ w + d`, `d ≥ 2η`,
`|∫ r ψ e^{−iζu}| ≤ (2κ)^{M−1} · 16 (w+η) B₀⁴ ℓ (ℓ + 4A₀/η) (2(M+1+8A₀)N²/(ℓd))^N · sup_{[0,ℓ]}|r|`. -/
theorem lemmaLambda (hT : Nazarov1993_Turan) : ∃ κ : ℝ, 0 < κ ∧
    ∀ (r : ExpSum) (w η ℓ d ζ K : ℝ) (N : ℕ), 0 ≤ w → 0 < η → 0 < ℓ → 1 ≤ r.M →
      (∀ j, r.node j ∈ Icc (-w) w) → 2 * η ≤ d → w + d ≤ |ζ| →
      (r.M : ℝ) * (4 * A₀ / η) * ((r.M : ℝ) + 1) ^ 2 ≤ ℓ →
      (∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ K) →
      ‖∫ u : ℝ, r.eval u * modPsi ℓ ζ u‖ ≤
        (2 * κ) ^ (r.M - 1) * (16 * (w + η) * ((B₀ * B₀) * (B₀ * B₀)) * ℓ * (ℓ + 4 * A₀ / η) *
          (2 * ((r.M : ℝ) + 1 + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * d)) ^ N) * K := by
  obtain ⟨κ, hκ, hml⟩ := ml_bound hT
  refine ⟨κ, hκ, fun r w η ℓ d ζ K N hw hη hℓ hM hnodes hd hζ hℓbig hK => ?_⟩
  set M := r.M
  set Aχ := 2 * (2 / η * A₀)
  have hAχeq : Aχ = 4 * A₀ / η := by simp only [Aχ]; ring
  have hA0 := A₀_pos
  have hB0 := B₀_pos
  have hAχ : 0 < Aχ := by positivity
  have hd0 : 0 < d := by linarith
  have hK0 : 0 ≤ K := (norm_nonneg _).trans (hK 0 ⟨le_rfl, hℓ.le⟩)
  rw [r.integral_mul_modPsi hℓ ζ, ← r.Lam_zero_lamG hη hnodes]
  set R₀ := ℓ + 4 * Aχ
  set Y := ℓ + Aχ * ((M : ℝ) + 1) ^ 2
  have hR₀ : 0 < R₀ := by positivity
  have hY : 0 < Y := by positivity
  have hmain := hml r (lamG_isTest (ζ := ζ) (w := w) hℓ hη) ℓ K 0 R₀ Y hℓ hM ⟨le_rfl, hℓ.le⟩ hR₀ hY hK
  set E := (2 * ((M : ℝ) + 1 + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * d)) ^ N
  set Λ := 2 * (w + η) * ((B₀ * B₀) * (B₀ * B₀) * ℓ * E)
  have hΛ : 0 ≤ Λ := by positivity
  have hEn : ∀ n : ℕ, n ≤ M + 1 →
      (2 * ((n : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * d)) ^ N ≤ E := by
    intro n hn
    apply pow_le_pow_left₀ (by positivity)
    gcongr
    have : (n : ℝ) ≤ (M : ℝ) + 1 := by exact_mod_cast hn
    linarith
  have hA : ∀ n : ℕ, n ≤ M + 1 → ∫ θ : ℝ, ‖iteratedDeriv n (lamG ℓ ζ w η) θ‖ ≤
      Λ * (ℓ + Aχ * (n : ℝ) ^ 2) ^ n := by
    intro n hn
    refine (integral_norm_iteratedDeriv_lamG_le hℓ hw hη hd hζ n N).trans ?_
    simp only [Λ]
    have := hEn n hn
    have hp : 0 ≤ (ℓ + Aχ * (n : ℝ) ^ 2) ^ n := by positivity
    calc 2 * (w + η) * (B₀ * B₀ * (B₀ * B₀) * ℓ *
          (2 * ((n : ℝ) + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * d)) ^ N * (ℓ + Aχ * (n : ℝ) ^ 2) ^ n)
        ≤ 2 * (w + η) * (B₀ * B₀ * (B₀ * B₀) * ℓ * E * (ℓ + Aχ * (n : ℝ) ^ 2) ^ n) := by gcongr
      _ = _ := by ring
  have hM1 : 1 ≤ (M : ℝ) := by exact_mod_cast hM
  have hA0' : (∫ θ : ℝ, ‖lamG ℓ ζ w η θ‖) ≤ Λ := by
    have := hA 0 (by omega); simpa using this
  have hA2 : (∫ θ : ℝ, ‖iteratedDeriv 2 (lamG ℓ ζ w η) θ‖) ≤ Λ * R₀ ^ 2 := by
    refine (hA 2 (by omega)).trans ?_
    have hbase : ℓ + Aχ * ((2 : ℕ) : ℝ) ^ 2 ≤ R₀ := le_of_eq (by simp only [R₀]; push_cast; ring)
    exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by positivity) hbase _) hΛ
  have hAm1 : (∫ θ : ℝ, ‖iteratedDeriv (M - 1) (lamG ℓ ζ w η) θ‖) ≤ Λ * Y ^ (M - 1) := by
    refine (hA (M - 1) (by omega)).trans ?_
    have hle : ((M - 1 : ℕ) : ℝ) ≤ (M : ℝ) + 1 := by
      have : ((M - 1 : ℕ) : ℝ) ≤ (M : ℝ) := by exact_mod_cast Nat.sub_le M 1
      linarith
    have hsq : ((M - 1 : ℕ) : ℝ) ^ 2 ≤ ((M : ℝ) + 1) ^ 2 := pow_le_pow_left₀ (by positivity) hle 2
    have hbase : ℓ + Aχ * ((M - 1 : ℕ) : ℝ) ^ 2 ≤ Y := by
      simp only [Y]; have := mul_le_mul_of_nonneg_left hsq hAχ.le; linarith
    exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by positivity) hbase _) hΛ
  have hAp1 : (∫ θ : ℝ, ‖iteratedDeriv (M - 1 + 2) (lamG ℓ ζ w η) θ‖) ≤ Λ * Y ^ (M + 1) := by
    have e : M - 1 + 2 = M + 1 := by omega
    rw [e]
    refine (hA (M + 1) le_rfl).trans ?_
    simp only [Y]; push_cast; exact le_rfl
  have hYM : (ℓ ^ (M - 1))⁻¹ * (2 * Λ * Y ^ M) ≤ 2 * Real.exp 1 * Λ * ℓ := by
    have hYeq : Y = ℓ * (1 + Aχ * ((M : ℝ) + 1) ^ 2 / ℓ) := by simp only [Y]; field_simp
    have hx : 0 ≤ Aχ * ((M : ℝ) + 1) ^ 2 / ℓ := by positivity
    have hMx : (M : ℝ) * (Aχ * ((M : ℝ) + 1) ^ 2 / ℓ) ≤ 1 := by
      rw [← mul_div_assoc, div_le_one hℓ, ← mul_assoc, hAχeq]; exact hℓbig
    have hexp : (1 + Aχ * ((M : ℝ) + 1) ^ 2 / ℓ) ^ M ≤ Real.exp 1 :=
      (one_add_pow_le_exp hx M).trans (Real.exp_le_exp.mpr hMx)
    have hℓM : (ℓ ^ (M - 1))⁻¹ * ℓ ^ M = ℓ := by
      obtain ⟨m, hm⟩ : ∃ m, M = m + 1 := ⟨M - 1, by omega⟩
      rw [hm, Nat.add_sub_cancel, pow_succ]; field_simp
    rw [hYeq, mul_pow]
    calc (ℓ ^ (M - 1))⁻¹ * (2 * Λ * (ℓ ^ M * (1 + Aχ * ((M : ℝ) + 1) ^ 2 / ℓ) ^ M))
        = 2 * Λ * ((ℓ ^ (M - 1))⁻¹ * ℓ ^ M) * (1 + Aχ * ((M : ℝ) + 1) ^ 2 / ℓ) ^ M := by ring
      _ = 2 * Λ * ℓ * (1 + Aχ * ((M : ℝ) + 1) ^ 2 / ℓ) ^ M := by rw [hℓM]
      _ ≤ 2 * Λ * ℓ * Real.exp 1 := by gcongr
      _ = 2 * Real.exp 1 * Λ * ℓ := by ring
  have he : Real.exp 1 ≤ 3 := by have := Real.exp_one_lt_d9; linarith
  have hbr : R₀ * (∫ θ : ℝ, ‖lamG ℓ ζ w η θ‖) + (∫ θ : ℝ, ‖iteratedDeriv 2 (lamG ℓ ζ w η) θ‖) / R₀ +
      (ℓ ^ (M - 1))⁻¹ * (Y * (∫ θ : ℝ, ‖iteratedDeriv (M - 1) (lamG ℓ ζ w η) θ‖) +
        (∫ θ : ℝ, ‖iteratedDeriv (M - 1 + 2) (lamG ℓ ζ w η) θ‖) / Y) ≤ 8 * Λ * (ℓ + Aχ) := by
    have h1 : R₀ * (∫ θ : ℝ, ‖lamG ℓ ζ w η θ‖) ≤ R₀ * Λ := by gcongr
    have h2 : (∫ θ : ℝ, ‖iteratedDeriv 2 (lamG ℓ ζ w η) θ‖) / R₀ ≤ Λ * R₀ := by
      rw [div_le_iff₀ hR₀]
      calc (∫ θ : ℝ, ‖iteratedDeriv 2 (lamG ℓ ζ w η) θ‖) ≤ Λ * R₀ ^ 2 := hA2
        _ = Λ * R₀ * R₀ := by ring
    have h3 : Y * (∫ θ : ℝ, ‖iteratedDeriv (M - 1) (lamG ℓ ζ w η) θ‖) +
        (∫ θ : ℝ, ‖iteratedDeriv (M - 1 + 2) (lamG ℓ ζ w η) θ‖) / Y ≤ 2 * Λ * Y ^ M := by
      have hYp : Y * Y ^ (M - 1) = Y ^ M := by rw [← pow_succ']; congr 1; omega
      have h3a : Y * (∫ θ : ℝ, ‖iteratedDeriv (M - 1) (lamG ℓ ζ w η) θ‖) ≤ Λ * Y ^ M := by
        calc Y * (∫ θ : ℝ, ‖iteratedDeriv (M - 1) (lamG ℓ ζ w η) θ‖) ≤ Y * (Λ * Y ^ (M - 1)) := by gcongr
          _ = Λ * Y ^ M := by rw [← hYp]; ring
      have h3b : (∫ θ : ℝ, ‖iteratedDeriv (M - 1 + 2) (lamG ℓ ζ w η) θ‖) / Y ≤ Λ * Y ^ M := by
        rw [div_le_iff₀ hY]
        calc (∫ θ : ℝ, ‖iteratedDeriv (M - 1 + 2) (lamG ℓ ζ w η) θ‖) ≤ Λ * Y ^ (M + 1) := hAp1
          _ = Λ * Y ^ M * Y := by rw [pow_succ]; ring
      linarith
    have h4 : (ℓ ^ (M - 1))⁻¹ * (Y * (∫ θ : ℝ, ‖iteratedDeriv (M - 1) (lamG ℓ ζ w η) θ‖) +
        (∫ θ : ℝ, ‖iteratedDeriv (M - 1 + 2) (lamG ℓ ζ w η) θ‖) / Y) ≤ 2 * Real.exp 1 * Λ * ℓ :=
      (mul_le_mul_of_nonneg_left h3 (by positivity)).trans hYM
    have e1 : R₀ * Λ = Λ * ℓ + 4 * (Λ * Aχ) := by simp only [R₀]; ring
    have e2 : Λ * R₀ = Λ * ℓ + 4 * (Λ * Aχ) := by simp only [R₀]; ring
    have e3 : 8 * Λ * (ℓ + Aχ) = 8 * (Λ * ℓ) + 8 * (Λ * Aχ) := by ring
    have h5 : 2 * Real.exp 1 * Λ * ℓ ≤ 6 * (Λ * ℓ) := by
      have hΛℓ : 0 ≤ Λ * ℓ := mul_nonneg hΛ hℓ.le
      calc 2 * Real.exp 1 * Λ * ℓ = 2 * Real.exp 1 * (Λ * ℓ) := by ring
        _ ≤ 2 * 3 * (Λ * ℓ) := by gcongr
        _ = 6 * (Λ * ℓ) := by ring
    rw [e3]; linarith
  calc ‖r.Lam 0 (lamG ℓ ζ w η)‖ ≤ κ ^ (M - 1) * K * 2 ^ (M - 1) * (R₀ * (∫ θ : ℝ, ‖lamG ℓ ζ w η θ‖) +
        (∫ θ : ℝ, ‖iteratedDeriv 2 (lamG ℓ ζ w η) θ‖) / R₀ +
          (ℓ ^ (M - 1))⁻¹ * (Y * (∫ θ : ℝ, ‖iteratedDeriv (M - 1) (lamG ℓ ζ w η) θ‖) +
            (∫ θ : ℝ, ‖iteratedDeriv (M - 1 + 2) (lamG ℓ ζ w η) θ‖) / Y)) := hmain
    _ ≤ κ ^ (M - 1) * K * 2 ^ (M - 1) * (8 * Λ * (ℓ + Aχ)) := by gcongr
    _ = (2 * κ) ^ (M - 1) * (16 * (w + η) * ((B₀ * B₀) * (B₀ * B₀)) * ℓ * (ℓ + 4 * A₀ / η) * E) * K := by
        rw [mul_pow, hAχeq]; simp only [Λ]; ring

end Riemann.Deep.LeanNBRate
