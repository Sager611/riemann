import Riemann.Deep.Hardy.Kernel

/-!
# Hardy's theorem, part 2: Λ on the critical line is a Fourier transform of the kernel

Λ(1/2 + i t) = (1/2) · 𝓕 Hk (t / (4π)), i.e. Λ(1/2+it) = (1/2) ∫ Hk(v) e^{-ivt/2} dv.

Route: in the strip 0 < Re s < 1, Λ(s) = (1/2) mellin (θ - 1 - x^{-1/2}) (s/2)
(Mathlib's Λ₀ is the Mellin transform of a piecewise-modified theta kernel; subtracting the
Mellin transforms of 1_{(0,1]} and x^{-1/2} 1_{(1,∞)} removes the two pole terms), and
mellin_eq_fourier turns the Mellin transform into a Fourier transform.
Also: Λ(1/2+it) is even in t and real.
-/

open Complex Real MeasureTheory Set Filter Topology HurwitzZeta
open scoped ComplexConjugate FourierTransform

namespace Fayed.Hardy

lemma hasMellin_Ioi_cpow {w : ℂ} (hw : w.re < 1 / 2) :
    HasMellin (indicator (Ioi 1) (fun x : ℝ => (x : ℂ) ^ (-(1 / 2 : ℂ)))) w (1 / (1 / 2 - w)) := by
  have ha : (w - 3 / 2).re < -1 := by
    simp only [sub_re]; norm_num; linarith
  have key : ∀ x ∈ Ioi (1 : ℝ),
      (x : ℂ) ^ (w - 1) • (x : ℂ) ^ (-(1 / 2 : ℂ)) = (x : ℂ) ^ (w - 3 / 2) := by
    intro x hx
    have hx0 : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr (by linarith [mem_Ioi.mp hx])
    rw [smul_eq_mul, ← cpow_add _ _ hx0]
    ring_nf
  have hint : IntegrableOn (fun x : ℝ => (x : ℂ) ^ (w - 3 / 2)) (Ioi 1) :=
    integrableOn_Ioi_cpow_of_lt ha one_pos
  simp_rw [HasMellin, mellin, MellinConvergent, ← indicator_smul, IntegrableOn,
    integrable_indicator_iff measurableSet_Ioi, integral_indicator measurableSet_Ioi,
    Measure.restrict_restrict_of_subset (Ioi_subset_Ioi zero_le_one)]
  refine ⟨?_, ?_⟩
  · rw [IntegrableOn, Measure.restrict_restrict_of_subset (Ioi_subset_Ioi zero_le_one)]
    exact hint.congr_fun (fun x hx => (key x hx).symm) measurableSet_Ioi
  · rw [setIntegral_congr_fun measurableSet_Ioi key, integral_Ioi_cpow_of_lt ha one_pos]
    have hw' : w - 3 / 2 + 1 ≠ 0 := by
      intro h; have := congrArg Complex.re h; simp at this; linarith
    have hw'' : (1 / 2 : ℂ) - w ≠ 0 := by
      intro h; have := congrArg Complex.re h; simp at this; linarith
    rw [ofReal_one, one_cpow]
    field_simp
    ring_nf
    rw [show (1 : ℂ) - w * 2 = -(-1 + w * 2) by ring, inv_neg]

/-- The modified theta kernel whose Mellin transform is 2 Λ(2w) in 0 < Re w < 1/2. -/
noncomputable def Gm (x : ℝ) : ℂ :=
  (hurwitzEvenFEPair 0).f_modif x - indicator (Ioc 0 1) (fun _ => (1 : ℂ)) x -
    indicator (Ioi 1) (fun x : ℝ => (x : ℂ) ^ (-(1 / 2 : ℂ))) x

lemma completedRiemannZeta_eq_mellin {s : ℂ} (h0 : 0 < s.re) (h1 : s.re < 1) :
    completedRiemannZeta s = mellin Gm (s / 2) / 2 := by
  set P := hurwitzEvenFEPair 0
  have hw0 : 0 < (s / 2).re := by simp [div_ofNat_re]; linarith
  have hw1 : (s / 2).re < 1 / 2 := by simp [div_ofNat_re]; linarith
  have hc1 : MellinConvergent P.f_modif (s / 2) := (P.isStrongFEPair_toStrongFEPair.hasMellin _).1
  have hm2 := hasMellin_one_Ioc hw0
  have hm3 := hasMellin_Ioi_cpow hw1
  have hA := hasMellin_sub hc1 hm2.1
  have hB := hasMellin_sub hA.1 hm3.1
  have hG : mellin Gm (s / 2) = mellin P.f_modif (s / 2) - 1 / (s / 2) - 1 / (1 / 2 - s / 2) := by
    have : Gm = fun t => (P.f_modif t - indicator (Ioc 0 1) (fun _ => (1 : ℂ)) t) -
        indicator (Ioi 1) (fun x : ℝ => (x : ℂ) ^ (-(1 / 2 : ℂ))) t := rfl
    rw [this, hB.2, hA.2, hm2.2, hm3.2]
  have hs0 : s ≠ 0 := by rintro rfl; simp at h0
  have hs1 : (1 : ℂ) - s ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this; linarith
  have h0' : completedRiemannZeta₀ s = mellin P.f_modif (s / 2) / 2 := rfl
  rw [completedRiemannZeta_eq, h0', hG]
  field_simp

lemma f_modif_lt_one {x : ℝ} (hx0 : 0 < x) (hx : x < 1) :
    (hurwitzEvenFEPair 0).f_modif x = (evenKernel 0 x : ℂ) - ((x ^ (-(1 / 2 : ℝ)) : ℝ) : ℂ) := by
  simp only [WeakFEPair.f_modif, Pi.add_apply]
  rw [indicator_of_notMem (by simp; linarith), indicator_of_mem (mem_Ioo.mpr ⟨hx0, hx⟩)]
  simp [hurwitzEvenFEPair]

lemma f_modif_gt_one {x : ℝ} (hx : 1 < x) :
    (hurwitzEvenFEPair 0).f_modif x = (evenKernel 0 x : ℂ) - 1 := by
  simp only [WeakFEPair.f_modif, Pi.add_apply]
  rw [indicator_of_mem (mem_Ioi.mpr hx), indicator_of_notMem (by simp; intro; linarith)]
  simp [hurwitzEvenFEPair]

/-- Pointwise identification of the Fourier integrand with the kernel Hk. -/
lemma Gm_exp_eq {v : ℝ} (hv : v ≠ 0) :
    rexp (-(1 / 4) * v) • Gm (rexp (-v)) = Hk v := by
  have hθ : ((evenKernel 0 (rexp (-v)) : ℝ) : ℂ) = theta (cexp (-(v : ℂ))) := by
    have := evenKernel_def 0 (rexp (-v))
    simp only [QuotientAddGroup.mk_zero, ofReal_zero] at this
    rw [this, theta]; simp [ofReal_exp]
  have hpow : ((rexp (-v) : ℝ) : ℂ) ^ (-(1 / 2 : ℂ)) = cexp ((v : ℂ) / 2) := by
    rw [ofReal_exp, cpow_def_of_ne_zero (Complex.exp_ne_zero _), Complex.log_exp (by simp [pi_pos])
      (by simp [pi_pos.le])]
    congr 1; push_cast; ring
  have hrpow : ((rexp (-v) ^ (-(1 / 2 : ℝ)) : ℝ) : ℂ) = cexp ((v : ℂ) / 2) := by
    rw [ofReal_cpow (Real.exp_pos _).le]; push_cast at hpow ⊢; exact hpow
  have hG : Gm (rexp (-v)) = theta (cexp (-(v : ℂ))) - 1 - cexp ((v : ℂ) / 2) := by
    rcases lt_or_gt_of_ne hv with h | h
    · have hx : 1 < rexp (-v) := by rw [Real.one_lt_exp_iff]; linarith
      rw [Gm, f_modif_gt_one hx, indicator_of_notMem (by simp; intro; linarith),
        indicator_of_mem (mem_Ioi.mpr hx), hpow, hθ]
      ring
    · have hx : rexp (-v) < 1 := by
        have := Real.exp_lt_exp.mpr (show -v < 0 by linarith); simpa using this
      have hx0 : 0 < rexp (-v) := Real.exp_pos _
      rw [Gm, f_modif_lt_one hx0 hx, indicator_of_mem (mem_Ioc.mpr ⟨hx0, hx.le⟩),
        indicator_of_notMem (by simp; linarith), hrpow, hθ]
      ring
  rw [hG, Complex.real_smul]
  simp only [Hk, Rk, Ck, Complex.two_cosh]
  push_cast
  rw [show cexp (-(1 / 4) * (v : ℂ)) = cexp (-(v : ℂ) / 4) by congr 1; ring]
  have : cexp (-(v : ℂ) / 4) * cexp ((v : ℂ) / 2) = cexp ((v : ℂ) / 4) := by
    rw [← Complex.exp_add]; congr 1; ring
  rw [show -(v : ℂ) / 4 = -((v : ℂ) / 4) by ring] at this ⊢
  linear_combination (-1 : ℂ) * this

/-- **Λ on the critical line is a Fourier transform of Hk.** -/
theorem Lam_eq_fourier (t : ℝ) :
    completedRiemannZeta (1 / 2 + t * I) = (1 / 2) * 𝓕 (fun v : ℝ => Hk v) (t / (4 * π)) := by
  have hre : ((1 / 2 + t * I : ℂ) / 2).re = 1 / 4 := by simp [div_ofNat_re]; norm_num
  have him : ((1 / 2 + t * I : ℂ) / 2).im = t / 2 := by simp [div_ofNat_im]
  rw [completedRiemannZeta_eq_mellin (by simp) (by simp; norm_num), mellin_eq_fourier, hre, him,
    show t / 2 / (2 * π) = t / (4 * π) by ring, div_eq_inv_mul, ← one_div]
  congr 1
  rw [Real.fourier_real_eq, Real.fourier_real_eq]
  refine integral_congr_ae ?_
  have hae : ∀ᵐ v : ℝ, v ≠ 0 := by
    rw [ae_iff]; simp [measure_singleton]
  filter_upwards [hae] with v hv
  rw [Gm_exp_eq hv]

/-- Λ(1/2 + i t) is even in t. -/
lemma Lam_neg (t : ℝ) :
    completedRiemannZeta (1 / 2 + (-t : ℝ) * I) = completedRiemannZeta (1 / 2 + t * I) := by
  rw [← completedRiemannZeta_one_sub]; congr 1; push_cast; ring

/-- Λ(1/2 + i t) is real. -/
lemma Lam_conj (t : ℝ) :
    conj (completedRiemannZeta (1 / 2 + t * I)) = completedRiemannZeta (1 / 2 + t * I) := by
  conv_rhs => rw [← Lam_neg t]
  rw [Lam_eq_fourier, Lam_eq_fourier, map_mul, Real.fourier_real_eq_integral_exp_smul,
    Real.fourier_real_eq_integral_exp_smul, ← integral_conj]
  congr 1
  · rw [map_div₀, map_one, map_ofNat]
  · congr 1; ext v
    rw [smul_eq_mul, smul_eq_mul, map_mul, conj_Hk_ofReal, ← exp_conj]
    congr 2
    simp only [map_mul, conj_ofReal, conj_I]
    push_cast; ring

end Fayed.Hardy
