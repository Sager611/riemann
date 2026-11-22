import Riemann.Deep.LeanNBRate.LemmaLambda
import Riemann.Deep.LeanNBRate.LemmaB
import Mathlib.Algebra.Order.Chebyshev

/-!
# [N] Lemma D: almost orthogonality of the groups

Crude-count version: with G groups whose windows [c_g - w_g, c_g + w_g] are s₀-separated,
(1 - ((1+η⁻¹)G + 1)·2π·δ) Σ_g ∫|q_g ψ|² ≤ (1+η) ∫|Σ_g q_g ψ|², where δ bounds the share of
∫|kern(q_g ψ)|² outside distance s₀/2 of the window (Lemma Λ). Plancherel: Mathlib's
SchwartzMap.integral_norm_sq_fourier.
-/

open MeasureTheory Filter Topology Set Complex FourierTransform

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem IsTest.fourier_eq_toS {G : ℝ → ℂ} (hG : IsTest G) :
    𝓕 G = ((𝓕 hG.toS : SchwartzMap ℝ ℂ) : ℝ → ℂ) := by
  rw [SchwartzMap.fourier_coe]; rfl

theorem IsTest.integral_norm_sq_fourier {G : ℝ → ℂ} (hG : IsTest G) :
    ∫ ξ : ℝ, ‖𝓕 G ξ‖ ^ 2 = ∫ u : ℝ, ‖G u‖ ^ 2 := by
  rw [hG.fourier_eq_toS, SchwartzMap.integral_norm_sq_fourier hG.toS]; rfl

theorem IsTest.integrable_norm_sq {G : ℝ → ℂ} (hG : IsTest G) :
    Integrable (fun u : ℝ => ‖G u‖ ^ 2) :=
  (hG.toS.memLp 2 volume).norm.integrable_sq

theorem IsTest.integrable_norm_sq_fourier {G : ℝ → ℂ} (hG : IsTest G) :
    Integrable (fun ξ : ℝ => ‖𝓕 G ξ‖ ^ 2) := by
  rw [hG.fourier_eq_toS]
  exact ((𝓕 hG.toS).memLp 2 volume).norm.integrable_sq

theorem norm_two_pi_inv : ‖(2 * Real.pi : ℂ)⁻¹‖ = (2 * Real.pi)⁻¹ := by
  rw [norm_inv]; congr 1
  rw [show (2 * Real.pi : ℂ) = ((2 * Real.pi : ℝ) : ℂ) by push_cast; ring, Complex.norm_real,
    Real.norm_eq_abs, abs_of_pos (by positivity)]

theorem norm_sq_kern (G : ℝ → ℂ) (v : ℝ) :
    ‖kern G v‖ ^ 2 = (2 * Real.pi)⁻¹ ^ 2 * ‖𝓕 G (v / (2 * Real.pi))‖ ^ 2 := by
  unfold kern; rw [norm_mul, norm_two_pi_inv, mul_pow]

/-- Plancherel for kernels of test functions. -/
theorem IsTest.integral_norm_sq_kern {G : ℝ → ℂ} (hG : IsTest G) :
    ∫ v : ℝ, ‖kern G v‖ ^ 2 = (2 * Real.pi)⁻¹ * ∫ u : ℝ, ‖G u‖ ^ 2 := by
  simp_rw [norm_sq_kern]
  rw [integral_const_mul, Measure.integral_comp_div (fun ξ => ‖𝓕 G ξ‖ ^ 2) (2 * Real.pi),
    hG.integral_norm_sq_fourier, smul_eq_mul, abs_of_pos (by positivity : (0:ℝ) < 2 * Real.pi)]
  have hπ : (2 * Real.pi) ≠ 0 := by positivity
  field_simp

theorem IsTest.integrable_norm_sq_kern {G : ℝ → ℂ} (hG : IsTest G) :
    Integrable (fun v : ℝ => ‖kern G v‖ ^ 2) := by
  simp_rw [norm_sq_kern]
  exact (hG.integrable_norm_sq_fourier.comp_div (by positivity : (2 * Real.pi) ≠ 0)).const_mul _

theorem integrable_exp_mul {Q : ℝ → ℂ} (hQ : Integrable Q) (v : ℝ) :
    Integrable (fun θ : ℝ => cexp (-(I * θ * v)) * Q θ) := by
  refine (hQ.mul_bdd (g := fun θ : ℝ => cexp (-(I * θ * v))) (c := 1) ?_ ?_).congr
    (Eventually.of_forall fun θ => by simp only [mul_comm])
  · exact (by fun_prop : Continuous fun θ : ℝ => cexp (-(I * θ * v))).aestronglyMeasurable
  · refine Eventually.of_forall fun θ => le_of_eq ?_
    rw [show -(I * (θ : ℂ) * v) = ((-(θ * v) : ℝ) : ℂ) * I by push_cast; ring,
      Complex.norm_exp_ofReal_mul_I]

theorem kern_finsetSum {ι : Type*} (s : Finset ι) (Q : ι → ℝ → ℂ)
    (hQ : ∀ i ∈ s, Integrable (Q i)) (v : ℝ) :
    kern (fun u => ∑ i ∈ s, Q i u) v = ∑ i ∈ s, kern (Q i) v := by
  simp only [kern_eq]
  rw [← Finset.mul_sum, ← integral_finsetSum s fun i hi => integrable_exp_mul (hQ i hi) v]
  congr 1
  refine integral_congr_ae (Eventually.of_forall fun θ => ?_)
  simp only [Finset.mul_sum]

theorem ExpSum.contDiff_eval (r : ExpSum) : ContDiff ℝ (⊤ : ℕ∞) r.eval := by
  unfold ExpSum.eval
  refine ContDiff.sum fun j _ => ContDiff.sum fun k _ => ?_
  refine (contDiff_const.mul ((contDiff_const.mul ofRealCLM.contDiff).pow k)).mul ?_
  exact (contDiff_const.mul ofRealCLM.contDiff).cexp

/-- The windowed derivative of one group: Q(u) = e^{icu} (r.dmod c)(u) ψ(u). -/
def grpQ (r : ExpSum) (c ℓ : ℝ) (u : ℝ) : ℂ :=
  cexp (I * c * u) * (r.dmod c).eval u * psiC ℓ u

theorem grpQ_isTest (r : ExpSum) (c : ℝ) {ℓ : ℝ} (hℓ : 0 < ℓ) : IsTest (grpQ r c ℓ) := by
  refine ⟨?_, ?_⟩
  · unfold grpQ
    exact (((contDiff_const.mul ofRealCLM.contDiff).cexp).mul (r.dmod c).contDiff_eval).mul
      (psiC_isTest hℓ).smooth
  · unfold grpQ
    exact (psiC_isTest hℓ).supp.mul_left

theorem norm_cexp_I_mul (x : ℝ) : ‖cexp (I * x)‖ = 1 := by
  rw [mul_comm, Complex.norm_exp_ofReal_mul_I]

theorem norm_grpQ (r : ExpSum) (c ℓ u : ℝ) :
    ‖grpQ r c ℓ u‖ = ‖(r.dmod c).eval u‖ * twin ℓ u := by
  unfold grpQ psiC
  rw [norm_mul, norm_mul, show (I * (c : ℂ) * u) = I * ((c * u : ℝ) : ℂ) by push_cast; ring,
    norm_cexp_I_mul, one_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (twin_nonneg _ _)]

/-- kern Q (ξ) = (2π)⁻¹ ∫ (r.dmod c)(u) ψ(u) e^{-i(ξ-c)u} du. -/
theorem kern_grpQ (r : ExpSum) (c ℓ ξ : ℝ) :
    kern (grpQ r c ℓ) ξ = (2 * Real.pi : ℂ)⁻¹ *
      ∫ u : ℝ, (r.dmod c).eval u * modPsi ℓ (ξ - c) u := by
  rw [kern_eq]
  congr 1
  refine integral_congr_ae (Eventually.of_forall fun u => ?_)
  unfold grpQ modPsi
  have : cexp (-(I * u * ξ)) * cexp (I * c * u) = cexp (-(I * ((ξ - c : ℝ) : ℂ) * u)) := by
    rw [← Complex.exp_add]; push_cast; ring_nf
  calc cexp (-(I * u * ξ)) * (cexp (I * c * u) * (r.dmod c).eval u * psiC ℓ u)
      = (cexp (-(I * u * ξ)) * cexp (I * c * u)) * (r.dmod c).eval u * psiC ℓ u := by ring
    _ = _ := by rw [this]; ring

/-- Sup over [0, ℓ] by the mean square over the middle half (Turán + Nikolskii). -/
theorem sup_sq_le (hT : Nazarov1993_Turan) : ∃ κ : ℝ, 0 < κ ∧
    ∀ (r : ExpSum) (ℓ : ℝ), 0 < ℓ → 1 ≤ r.M → ∀ u ∈ Icc 0 ℓ,
      ‖r.eval u‖ ^ 2 ≤ 4 * (2 * κ) ^ (4 * (r.M - 1)) / ℓ *
        ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖r.eval x‖ ^ 2 := by
  obtain ⟨κ₁, hκ₁, hTm⟩ := turan_middle hT
  obtain ⟨κ₂, hκ₂, hN⟩ := nikolskii hT
  set κ := max κ₁ κ₂
  have hκ : 0 < κ := lt_max_of_lt_left hκ₁
  refine ⟨κ, hκ, fun r ℓ hℓ hM u hu => ?_⟩
  set I₀ := ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖r.eval x‖ ^ 2
  have hI₀ : 0 ≤ I₀ := setIntegral_nonneg measurableSet_Icc fun x _ => by positivity
  set X := Real.sqrt (2 * (2 * κ₂) ^ (2 * (r.M - 1)) * ((3 * ℓ / 4 - ℓ / 4)⁻¹ * I₀))
  have hX : ∀ v ∈ Icc (ℓ / 4) (3 * ℓ / 4), ‖r.eval v‖ ≤ X := by
    intro v hv
    apply Real.le_sqrt_of_sq_le
    exact hN r (ℓ / 4) (3 * ℓ / 4) hM (by linarith) v hv
  have h1 := hTm r ℓ X hℓ hM hX u hu
  have hX0 : 0 ≤ X := Real.sqrt_nonneg _
  have hba : 3 * ℓ / 4 - ℓ / 4 = ℓ / 2 := by ring
  have hXsq : X ^ 2 = 2 * (2 * κ₂) ^ (2 * (r.M - 1)) * ((3 * ℓ / 4 - ℓ / 4)⁻¹ * I₀) :=
    Real.sq_sqrt (by rw [hba]; positivity)
  have hk1 : (2 * κ₁) ^ (r.M - 1) ≤ (2 * κ) ^ (r.M - 1) :=
    pow_le_pow_left₀ (by positivity) (by linarith [le_max_left κ₁ κ₂]) _
  have hk2 : (2 * κ₂) ^ (2 * (r.M - 1)) ≤ (2 * κ) ^ (2 * (r.M - 1)) :=
    pow_le_pow_left₀ (by positivity) (by linarith [le_max_right κ₁ κ₂]) _
  calc ‖r.eval u‖ ^ 2 ≤ ((2 * κ₁) ^ (r.M - 1) * X) ^ 2 :=
        pow_le_pow_left₀ (norm_nonneg _) h1 2
    _ = ((2 * κ₁) ^ (r.M - 1)) ^ 2 * X ^ 2 := by ring
    _ ≤ ((2 * κ) ^ (r.M - 1)) ^ 2 * (2 * (2 * κ) ^ (2 * (r.M - 1)) * ((3 * ℓ / 4 - ℓ / 4)⁻¹ * I₀)) := by
        rw [hXsq]
        have : 0 ≤ (3 * ℓ / 4 - ℓ / 4)⁻¹ * I₀ := by rw [hba]; positivity
        gcongr
    _ = 4 * (2 * κ) ^ (4 * (r.M - 1)) / ℓ * I₀ := by
        rw [show 3 * ℓ / 4 - ℓ / 4 = ℓ / 2 by ring, ← pow_mul,
          show 4 * (r.M - 1) = (r.M - 1) * 2 + 2 * (r.M - 1) by ring, pow_add]
        field_simp
        ring

/-- The pointwise almost-orthogonality inequality: at most one index outside A. -/
theorem sum_norm_sq_le_of_near {ι : Type*} [Fintype ι] [DecidableEq ι] (x : ι → ℂ)
    (A : Finset ι) (hA : ∀ g g', g ∉ A → g' ∉ A → g = g') {η : ℝ} (hη : 0 < η) :
    ∑ g, ‖x g‖ ^ 2 ≤ (1 + η) * ‖∑ g, x g‖ ^ 2 +
      ((1 + η⁻¹) * Fintype.card ι + 1) * ∑ g ∈ A, ‖x g‖ ^ 2 := by
  have hSA : 0 ≤ ∑ g ∈ A, ‖x g‖ ^ 2 := Finset.sum_nonneg fun g _ => by positivity
  have hc : 0 ≤ (1 + η⁻¹) * (Fintype.card ι : ℝ) := by positivity
  by_cases hex : ∃ g₀, g₀ ∉ A
  · obtain ⟨g₀, hg₀⟩ := hex
    have hAe : A = Finset.univ.erase g₀ := by
      ext g
      simp only [Finset.mem_erase, Finset.mem_univ, and_true]
      constructor
      · intro hg h; exact hg₀ (h ▸ hg)
      · intro hne; by_contra hg; exact hne (hA g g₀ hg hg₀)
    have hsplit : ∀ f : ι → ℂ, ∑ g, f g = f g₀ + ∑ g ∈ A, f g := by
      intro f; rw [hAe, Finset.add_sum_erase _ _ (Finset.mem_univ g₀)]
    have hsplit' : ∑ g, ‖x g‖ ^ 2 = ‖x g₀‖ ^ 2 + ∑ g ∈ A, ‖x g‖ ^ 2 := by
      rw [hAe, Finset.add_sum_erase _ (fun g => ‖x g‖ ^ 2) (Finset.mem_univ g₀)]
    have hx0 : x g₀ = (∑ g, x g) - ∑ g ∈ A, x g := by rw [hsplit x]; ring
    set a := ∑ g, x g
    set b := ∑ g ∈ A, x g
    have hab : ‖a - b‖ ^ 2 ≤ (1 + η) * ‖a‖ ^ 2 + (1 + η⁻¹) * ‖b‖ ^ 2 := by
      have h1 : ‖a - b‖ ≤ ‖a‖ + ‖b‖ := norm_sub_le a b
      have h2 : 2 * ‖a‖ * ‖b‖ ≤ η * ‖a‖ ^ 2 + η⁻¹ * ‖b‖ ^ 2 := by
        have h0 := sq_nonneg (η * ‖a‖ - ‖b‖)
        have hη' : η⁻¹ * (η * ‖a‖ - ‖b‖) ^ 2 = η * ‖a‖ ^ 2 - 2 * ‖a‖ * ‖b‖ + η⁻¹ * ‖b‖ ^ 2 := by
          field_simp; ring
        have := mul_nonneg (inv_nonneg.mpr hη.le) h0
        linarith
      calc ‖a - b‖ ^ 2 ≤ (‖a‖ + ‖b‖) ^ 2 := pow_le_pow_left₀ (norm_nonneg _) h1 2
        _ = ‖a‖ ^ 2 + 2 * ‖a‖ * ‖b‖ + ‖b‖ ^ 2 := by ring
        _ ≤ _ := by linarith
    have hb : ‖b‖ ^ 2 ≤ (Fintype.card ι : ℝ) * ∑ g ∈ A, ‖x g‖ ^ 2 := by
      have h1 : ‖b‖ ≤ ∑ g ∈ A, ‖x g‖ := norm_sum_le _ _
      have h2 : (∑ g ∈ A, ‖x g‖) ^ 2 ≤ (A.card : ℝ) * ∑ g ∈ A, ‖x g‖ ^ 2 :=
        sq_sum_le_card_mul_sum_sq
      have h3 : (A.card : ℝ) ≤ Fintype.card ι := by exact_mod_cast Finset.card_le_univ A
      calc ‖b‖ ^ 2 ≤ (∑ g ∈ A, ‖x g‖) ^ 2 := pow_le_pow_left₀ (norm_nonneg _) h1 2
        _ ≤ (A.card : ℝ) * ∑ g ∈ A, ‖x g‖ ^ 2 := h2
        _ ≤ _ := by gcongr
    rw [hsplit', hx0]
    have h4 : (1 + η⁻¹) * ‖b‖ ^ 2 ≤ (1 + η⁻¹) * ((Fintype.card ι : ℝ) * ∑ g ∈ A, ‖x g‖ ^ 2) :=
      mul_le_mul_of_nonneg_left hb (by positivity)
    have e5 : ((1 + η⁻¹) * Fintype.card ι + 1) * ∑ g ∈ A, ‖x g‖ ^ 2 =
        (1 + η⁻¹) * ((Fintype.card ι : ℝ) * ∑ g ∈ A, ‖x g‖ ^ 2) + ∑ g ∈ A, ‖x g‖ ^ 2 := by ring
    rw [e5]; linarith
  · push Not at hex
    have hAu : A = Finset.univ := Finset.eq_univ_of_forall hex
    rw [hAu]
    rw [hAu] at hSA
    have h0 : 0 ≤ (1 + η) * ‖∑ g, x g‖ ^ 2 := by positivity
    have e5 : ((1 + η⁻¹) * Fintype.card ι + 1) * ∑ g, ‖x g‖ ^ 2 =
        (1 + η⁻¹) * Fintype.card ι * ∑ g, ‖x g‖ ^ 2 + ∑ g, ‖x g‖ ^ 2 := by ring
    have h1 : 0 ≤ (1 + η⁻¹) * Fintype.card ι * ∑ g, ‖x g‖ ^ 2 := mul_nonneg hc hSA
    rw [e5]; linarith

/-- Distance at least s₀/2 from the window [c - w, c + w]. -/
def farSet (c w s₀ : ℝ) : Set ℝ := {ξ | w + s₀ / 2 ≤ |ξ - c|}

theorem measurableSet_farSet (c w s₀ : ℝ) : MeasurableSet (farSet c w s₀) :=
  (isClosed_le continuous_const (continuous_abs.comp (continuous_id.sub continuous_const))).measurableSet

/-- The Lemma Λ constant with η = s₀/4, the sup-to-mean-square factor and the tail integral. -/
def tailC (κ : ℝ) (M : ℕ) (w ℓ s₀ : ℝ) (N : ℕ) : ℝ :=
  (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * (16 * (w + s₀ / 4) * ((B₀ * B₀) * (B₀ * B₀)) * ℓ *
    (ℓ + 4 * A₀ / (s₀ / 4)) * (2 * ((M : ℝ) + 1 + 8 * A₀) * (N : ℝ) ^ 2 / ℓ) ^ N)) ^ 2 *
  (4 * (2 * κ) ^ (4 * (M - 1)) / ℓ) * ((s₀ / 2) ^ (2 * N - 2))⁻¹ * (8 * Real.pi / s₀)

theorem integral_inv_sq_add_sq_sub {R : ℝ} (hR : 0 < R) (a : ℝ) :
    ∫ v : ℝ, (R ^ 2 + (v - a) ^ 2)⁻¹ = Real.pi / R := by
  rw [integral_sub_right_eq_self (fun v => (R ^ 2 + v ^ 2)⁻¹) a]
  exact integral_inv_sq_add_sq hR

theorem integrable_inv_sq_add_sq_sub {R : ℝ} (hR : 0 < R) (a : ℝ) :
    Integrable fun v : ℝ => (R ^ 2 + (v - a) ^ 2)⁻¹ :=
  (integrable_inv_sq_add_sq hR).comp_sub_right a

theorem far_decay_le {R d a b : ℝ} {N : ℕ} (hR : 0 < R) (hd : R ≤ d) (hN : 1 ≤ N)
    (hab : d ^ 2 = a ^ 2 ∨ d ^ 2 = b ^ 2) :
    ((d ^ N) ^ 2)⁻¹ ≤ (R ^ (2 * N - 2))⁻¹ * (2 * ((R ^ 2 + a ^ 2)⁻¹ + (R ^ 2 + b ^ 2)⁻¹)) := by
  have hd0 : 0 < d := lt_of_lt_of_le hR hd
  have hdN : ((d ^ N) ^ 2)⁻¹ ≤ (R ^ (2 * N - 2))⁻¹ * (d ^ 2)⁻¹ := by
    rw [← mul_inv]
    apply inv_anti₀ (by positivity)
    have e : (d ^ N) ^ 2 = d ^ (2 * N - 2) * d ^ 2 := by
      rw [← pow_mul, ← pow_add]; congr 1; omega
    rw [e]
    gcongr
  have hd2 : (d ^ 2)⁻¹ ≤ 2 * (R ^ 2 + d ^ 2)⁻¹ := by
    have hRd : R ^ 2 + d ^ 2 ≤ 2 * d ^ 2 := by nlinarith
    have hd2p : 0 < d ^ 2 := by positivity
    calc (d ^ 2)⁻¹ = 2 * (2 * d ^ 2)⁻¹ := by field_simp
      _ ≤ 2 * (R ^ 2 + d ^ 2)⁻¹ := by gcongr
  have hdsq : (R ^ 2 + d ^ 2)⁻¹ ≤ (R ^ 2 + a ^ 2)⁻¹ + (R ^ 2 + b ^ 2)⁻¹ := by
    have h1 : 0 ≤ (R ^ 2 + a ^ 2)⁻¹ := by positivity
    have h2 : 0 ≤ (R ^ 2 + b ^ 2)⁻¹ := by positivity
    rcases hab with h | h <;> rw [h] <;> linarith
  calc ((d ^ N) ^ 2)⁻¹ ≤ (R ^ (2 * N - 2))⁻¹ * (d ^ 2)⁻¹ := hdN
    _ ≤ (R ^ (2 * N - 2))⁻¹ * (2 * (R ^ 2 + d ^ 2)⁻¹) := by gcongr
    _ ≤ _ := by gcongr

set_option maxHeartbeats 1000000 in
/-- One group's frequency tail: the share of ∫|kern Q|² at distance ≥ s₀/2 from its window. -/
theorem tail_bound (hT : Nazarov1993_Turan) : ∃ κ : ℝ, 0 < κ ∧
    ∀ (r : ExpSum) (c w ℓ s₀ : ℝ) (N : ℕ), 0 ≤ w → 0 < s₀ → 0 < ℓ → 1 ≤ r.M → 1 ≤ N →
      (∀ j, r.node j ∈ Icc (-w) w) →
      (r.M : ℝ) * (4 * A₀ / (s₀ / 4)) * ((r.M : ℝ) + 1) ^ 2 ≤ ℓ →
      ∫ ξ : ℝ, (farSet c w s₀).indicator (fun ξ => ‖kern (grpQ r c ℓ) ξ‖ ^ 2) ξ ≤
        tailC κ r.M w ℓ s₀ N * ∫ u : ℝ, ‖grpQ r c ℓ u‖ ^ 2 := by
  obtain ⟨κ₁, hκ₁, hΛ⟩ := lemmaLambda hT
  obtain ⟨κ₂, hκ₂, hS⟩ := sup_sq_le hT
  set κ := max κ₁ κ₂
  have hκ : 0 < κ := lt_max_of_lt_left hκ₁
  refine ⟨κ, hκ, fun r c w ℓ s₀ N hw hs₀ hℓ hM hN hnodes hℓbig => ?_⟩
  set rq := r.dmod c
  have hMq : rq.M = r.M := r.dmod_M c
  have hnq : ∀ j, rq.node j ∈ Icc (-w) w := fun j => by rw [r.dmod_node c]; exact hnodes j
  set M := r.M
  set Q := grpQ r c ℓ
  have hQ := grpQ_isTest r c hℓ
  have hA0 := A₀_pos
  have hB0 := B₀_pos
  set η := s₀ / 4
  have hη : 0 < η := by positivity
  -- mean square of rq over the middle half ≤ ∫ |Q|²
  set Imid := ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖rq.eval x‖ ^ 2
  set IQ := ∫ u : ℝ, ‖Q u‖ ^ 2
  have hImid : Imid ≤ IQ := by
    have e : Imid = ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖Q x‖ ^ 2 := by
      refine setIntegral_congr_fun measurableSet_Icc fun x hx => ?_
      simp only [Q, rq]
      rw [norm_grpQ, twin_eq_one hℓ hx, mul_one]
    rw [e]
    exact setIntegral_le_integral hQ.integrable_norm_sq (Eventually.of_forall fun x => by positivity)
  have hImid0 : 0 ≤ Imid := setIntegral_nonneg measurableSet_Icc fun x _ => by positivity
  set K := Real.sqrt (4 * (2 * κ₂) ^ (4 * (M - 1)) / ℓ * Imid)
  have hK : ∀ u ∈ Icc 0 ℓ, ‖rq.eval u‖ ≤ K := by
    intro u hu
    apply Real.le_sqrt_of_sq_le
    have := hS rq ℓ hℓ (by rw [hMq]; exact hM) u hu
    rwa [hMq] at this
  have hK0 : 0 ≤ K := Real.sqrt_nonneg _
  have hKsq : K ^ 2 ≤ 4 * (2 * κ) ^ (4 * (M - 1)) / ℓ * IQ := by
    rw [Real.sq_sqrt (by positivity)]
    have h1 : (2 * κ₂) ^ (4 * (M - 1)) ≤ (2 * κ) ^ (4 * (M - 1)) :=
      pow_le_pow_left₀ (by positivity) (by linarith [le_max_right κ₁ κ₂]) _
    have h2 : 4 * (2 * κ₂) ^ (4 * (M - 1)) / ℓ ≤ 4 * (2 * κ) ^ (4 * (M - 1)) / ℓ := by gcongr
    calc 4 * (2 * κ₂) ^ (4 * (M - 1)) / ℓ * Imid ≤ 4 * (2 * κ) ^ (4 * (M - 1)) / ℓ * Imid := by gcongr
      _ ≤ _ := by gcongr
  -- the constants
  set C₁ := 16 * (w + η) * ((B₀ * B₀) * (B₀ * B₀)) * ℓ * (ℓ + 4 * A₀ / η)
  set β := 2 * ((M : ℝ) + 1 + 8 * A₀) * (N : ℝ) ^ 2 / ℓ
  have hC₁ : 0 ≤ C₁ := by positivity
  have hβ : 0 ≤ β := by positivity
  set R := s₀ / 2
  have hR : 0 < R := by positivity
  set P := (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * K ^ 2 * (R ^ (2 * N - 2))⁻¹ * 2
  have hP : 0 ≤ P := by positivity
  set g : ℝ → ℝ := fun ξ => P * ((R ^ 2 + (ξ - (c + w)) ^ 2)⁻¹ + (R ^ 2 + (ξ - (c - w)) ^ 2)⁻¹)
  have hg_int : Integrable g :=
    ((integrable_inv_sq_add_sq_sub hR (c + w)).add (integrable_inv_sq_add_sq_sub hR (c - w))).const_mul P
  have hg_val : ∫ ξ : ℝ, g ξ = P * (4 * Real.pi / s₀) := by
    simp only [g]
    rw [integral_const_mul, integral_add (integrable_inv_sq_add_sq_sub hR (c + w))
      (integrable_inv_sq_add_sq_sub hR (c - w)), integral_inv_sq_add_sq_sub hR,
      integral_inv_sq_add_sq_sub hR]
    simp only [R]; field_simp; ring
  -- pointwise domination
  have hpt : ∀ ξ, (farSet c w s₀).indicator (fun ξ => ‖kern Q ξ‖ ^ 2) ξ ≤ g ξ := by
    intro ξ
    by_cases hξ : ξ ∈ farSet c w s₀
    · rw [Set.indicator_of_mem hξ]
      have hξ' : w + s₀ / 2 ≤ |ξ - c| := hξ
      set d := |ξ - c| - w
      have hd : R ≤ d := by simp only [d, R]; linarith
      have hd0 : 0 < d := lt_of_lt_of_le hR hd
      have hΛξ := hΛ rq w η ℓ d (ξ - c) K N hw hη hℓ (by rw [hMq]; exact hM) hnq
        (by simp only [η, d]; linarith) (by simp only [d]; linarith)
        (by rw [hMq]; exact hℓbig) hK
      rw [hMq] at hΛξ
      have hk1 : (2 * κ₁) ^ (M - 1) ≤ (2 * κ) ^ (M - 1) :=
        pow_le_pow_left₀ (by positivity) (by linarith [le_max_left κ₁ κ₂]) _
      have hβd : 2 * ((M : ℝ) + 1 + 8 * A₀) * (N : ℝ) ^ 2 / (ℓ * d) = β / d := by
        simp only [β]; rw [div_div]
      rw [hβd] at hΛξ
      have hkern : ‖kern Q ξ‖ ≤ (2 * Real.pi)⁻¹ * ((2 * κ) ^ (M - 1) * C₁ * β ^ N * K / d ^ N) := by
        simp only [Q]
        rw [kern_grpQ, norm_mul, norm_two_pi_inv]
        gcongr
        refine hΛξ.trans ?_
        rw [div_pow]
        calc (2 * κ₁) ^ (M - 1) * (C₁ * (β ^ N / d ^ N)) * K
            ≤ (2 * κ) ^ (M - 1) * (C₁ * (β ^ N / d ^ N)) * K := by gcongr
          _ = (2 * κ) ^ (M - 1) * C₁ * β ^ N * K / d ^ N := by ring
      have hsq : ‖kern Q ξ‖ ^ 2 ≤ (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * K ^ 2 *
          ((d ^ N) ^ 2)⁻¹ := by
        calc ‖kern Q ξ‖ ^ 2 ≤ ((2 * Real.pi)⁻¹ * ((2 * κ) ^ (M - 1) * C₁ * β ^ N * K / d ^ N)) ^ 2 :=
              pow_le_pow_left₀ (norm_nonneg _) hkern 2
          _ = _ := by rw [div_eq_mul_inv]; ring
      have hab : d ^ 2 = (ξ - (c + w)) ^ 2 ∨ d ^ 2 = (ξ - (c - w)) ^ 2 := by
        rcases le_or_gt 0 (ξ - c) with h | h
        · left; simp only [d]; rw [abs_of_nonneg h]; ring
        · right; simp only [d]; rw [abs_of_neg h]; ring
      have hfar := far_decay_le hR hd hN hab
      have hX : 0 ≤ (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * K ^ 2 := by positivity
      calc ‖kern Q ξ‖ ^ 2 ≤ (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * K ^ 2 *
            ((d ^ N) ^ 2)⁻¹ := hsq
        _ ≤ (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * K ^ 2 *
            ((R ^ (2 * N - 2))⁻¹ * (2 * ((R ^ 2 + (ξ - (c + w)) ^ 2)⁻¹ +
              (R ^ 2 + (ξ - (c - w)) ^ 2)⁻¹))) := by
            exact mul_le_mul_of_nonneg_left hfar hX
        _ = g ξ := by simp only [g, P]; ring
    · rw [Set.indicator_of_notMem hξ]
      simp only [g]; positivity
  have hind_int : Integrable ((farSet c w s₀).indicator (fun ξ => ‖kern Q ξ‖ ^ 2)) :=
    hQ.integrable_norm_sq_kern.indicator (measurableSet_farSet c w s₀)
  calc ∫ ξ : ℝ, (farSet c w s₀).indicator (fun ξ => ‖kern Q ξ‖ ^ 2) ξ ≤ ∫ ξ : ℝ, g ξ :=
        integral_mono hind_int hg_int hpt
    _ = P * (4 * Real.pi / s₀) := hg_val
    _ ≤ tailC κ M w ℓ s₀ N * IQ := by
        simp only [P, tailC]
        have hIQ : 0 ≤ IQ := integral_nonneg fun u => by positivity
        have e1 : (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * K ^ 2 *
            (R ^ (2 * N - 2))⁻¹ * 2 * (4 * Real.pi / s₀) =
            ((2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * (R ^ (2 * N - 2))⁻¹ *
              (8 * Real.pi / s₀)) * K ^ 2 := by ring
        rw [e1]
        have hY : 0 ≤ (2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * (R ^ (2 * N - 2))⁻¹ *
            (8 * Real.pi / s₀) := by positivity
        calc ((2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * (R ^ (2 * N - 2))⁻¹ *
              (8 * Real.pi / s₀)) * K ^ 2
            ≤ ((2 * Real.pi)⁻¹ ^ 2 * ((2 * κ) ^ (M - 1) * C₁ * β ^ N) ^ 2 * (R ^ (2 * N - 2))⁻¹ *
              (8 * Real.pi / s₀)) * (4 * (2 * κ) ^ (4 * (M - 1)) / ℓ * IQ) :=
              mul_le_mul_of_nonneg_left hKsq hY
          _ = _ := by simp only [C₁, β, R, η]; ring

theorem IsTest.zero : IsTest (fun _ : ℝ => (0 : ℂ)) :=
  ⟨contDiff_const, by exact HasCompactSupport.zero⟩

theorem IsTest.add {F G : ℝ → ℂ} (hF : IsTest F) (hG : IsTest G) : IsTest (fun u => F u + G u) :=
  ⟨hF.smooth.add hG.smooth, hF.supp.add hG.supp⟩

theorem IsTest.finsetSum {ι : Type*} [DecidableEq ι] (s : Finset ι) {Q : ι → ℝ → ℂ}
    (hQ : ∀ i ∈ s, IsTest (Q i)) : IsTest (fun u => ∑ i ∈ s, Q i u) := by
  induction s using Finset.induction_on with
  | empty => simpa using IsTest.zero
  | insert a s ha ih =>
    have e : (fun u => ∑ i ∈ insert a s, Q i u) = fun u => Q a u + ∑ i ∈ s, Q i u := by
      funext u; rw [Finset.sum_insert ha]
    rw [e]
    exact (hQ a (Finset.mem_insert_self a s)).add (ih fun i hi => hQ i (Finset.mem_insert_of_mem hi))

/-- **[N] Lemma D** (almost orthogonality, crude-count form): for s₀-separated windows and
frequency tails bounded by δ, (1 - ((1+η⁻¹)G + 1)·2π·δ) Σ_g ∫|Q_g|² ≤ (1+η) ∫|Σ_g Q_g|². -/
theorem lemmaD {ι : Type*} [Fintype ι] [DecidableEq ι] (Q : ι → ℝ → ℂ) (hQ : ∀ g, IsTest (Q g))
    (c w : ι → ℝ) {s₀ : ℝ} (hs₀ : 0 < s₀)
    (hsep : ∀ g g', g ≠ g' → w g + w g' + s₀ ≤ |c g - c g'|) {δ : ℝ}
    (hδ : ∀ g, ∫ ξ : ℝ, (farSet (c g) (w g) s₀).indicator (fun ξ => ‖kern (Q g) ξ‖ ^ 2) ξ ≤
      δ * ∫ u : ℝ, ‖Q g u‖ ^ 2) {η : ℝ} (hη : 0 < η) :
    (1 - ((1 + η⁻¹) * Fintype.card ι + 1) * (2 * Real.pi) * δ) * ∑ g, ∫ u : ℝ, ‖Q g u‖ ^ 2 ≤
      (1 + η) * ∫ u : ℝ, ‖∑ g, Q g u‖ ^ 2 := by
  classical
  set C := (1 + η⁻¹) * (Fintype.card ι : ℝ) + 1
  have hC : 0 ≤ C := by positivity
  have hS := IsTest.finsetSum Finset.univ fun g _ => hQ g
  have hpt : ∀ ξ : ℝ, ∑ g, ‖kern (Q g) ξ‖ ^ 2 ≤ (1 + η) * ‖kern (fun u => ∑ g, Q g u) ξ‖ ^ 2 +
      C * ∑ g, (farSet (c g) (w g) s₀).indicator (fun ξ => ‖kern (Q g) ξ‖ ^ 2) ξ := by
    intro ξ
    set A := Finset.univ.filter fun g => ξ ∈ farSet (c g) (w g) s₀
    have hA : ∀ g g', g ∉ A → g' ∉ A → g = g' := by
      intro g g' hg hg'
      by_contra hne
      simp only [A, Finset.mem_filter, Finset.mem_univ, true_and, farSet, Set.mem_setOf_eq,
        not_le] at hg hg'
      have h := hsep g g' hne
      have : |c g - c g'| ≤ |ξ - c g| + |ξ - c g'| := by
        calc |c g - c g'| = |(ξ - c g') - (ξ - c g)| := by ring_nf
          _ ≤ |ξ - c g'| + |ξ - c g| := abs_sub _ _
          _ = _ := by ring
      linarith
    have h := sum_norm_sq_le_of_near (fun g => kern (Q g) ξ) A hA hη
    rw [kern_finsetSum Finset.univ Q (fun g _ => (hQ g).integrable)]
    have hind : ∑ g ∈ A, ‖kern (Q g) ξ‖ ^ 2 =
        ∑ g, (farSet (c g) (w g) s₀).indicator (fun ξ => ‖kern (Q g) ξ‖ ^ 2) ξ := by
      rw [Finset.sum_filter]
      refine Finset.sum_congr rfl fun g _ => ?_
      by_cases hg : ξ ∈ farSet (c g) (w g) s₀
      · rw [if_pos hg, Set.indicator_of_mem hg]
      · rw [if_neg hg, Set.indicator_of_notMem hg]
    rw [← hind]
    exact h
  have hind_int : ∀ g, Integrable ((farSet (c g) (w g) s₀).indicator (fun ξ => ‖kern (Q g) ξ‖ ^ 2)) :=
    fun g => (hQ g).integrable_norm_sq_kern.indicator (measurableSet_farSet _ _ _)
  have hL : Integrable fun ξ : ℝ => ∑ g, ‖kern (Q g) ξ‖ ^ 2 :=
    integrable_finsetSum _ fun g _ => (hQ g).integrable_norm_sq_kern
  have hR1 : Integrable fun ξ : ℝ => (1 + η) * ‖kern (fun u => ∑ g, Q g u) ξ‖ ^ 2 :=
    hS.integrable_norm_sq_kern.const_mul _
  have hR2 : Integrable fun ξ : ℝ =>
      C * ∑ g, (farSet (c g) (w g) s₀).indicator (fun ξ => ‖kern (Q g) ξ‖ ^ 2) ξ :=
    (integrable_finsetSum _ fun g _ => hind_int g).const_mul _
  have hint : ∫ ξ : ℝ, ∑ g, ‖kern (Q g) ξ‖ ^ 2 ≤ ∫ ξ : ℝ, ((1 + η) * ‖kern (fun u => ∑ g, Q g u) ξ‖ ^ 2 +
      C * ∑ g, (farSet (c g) (w g) s₀).indicator (fun ξ => ‖kern (Q g) ξ‖ ^ 2) ξ) :=
    integral_mono hL (hR1.add hR2) hpt
  rw [integral_finsetSum _ fun g _ => (hQ g).integrable_norm_sq_kern, integral_add hR1 hR2,
    integral_const_mul, integral_const_mul, integral_finsetSum _ fun g _ => hind_int g,
    hS.integral_norm_sq_kern] at hint
  simp_rw [IsTest.integral_norm_sq_kern (hQ _)] at hint
  have hsum : ∑ g, ∫ ξ : ℝ, (farSet (c g) (w g) s₀).indicator (fun ξ => ‖kern (Q g) ξ‖ ^ 2) ξ ≤
      δ * ∑ g, ∫ u : ℝ, ‖Q g u‖ ^ 2 := by
    rw [Finset.mul_sum]; exact Finset.sum_le_sum fun g _ => hδ g
  have hπ : 0 < 2 * Real.pi := by positivity
  set T := ∑ g, ∫ u : ℝ, ‖Q g u‖ ^ 2
  set U := ∫ u : ℝ, ‖∑ g, Q g u‖ ^ 2
  have e1 : ∑ g, (2 * Real.pi)⁻¹ * ∫ u : ℝ, ‖Q g u‖ ^ 2 = (2 * Real.pi)⁻¹ * T := by
    rw [← Finset.mul_sum]
  rw [e1] at hint
  have h2 : (2 * Real.pi)⁻¹ * T ≤ (1 + η) * ((2 * Real.pi)⁻¹ * U) + C * (δ * T) := by
    have := mul_le_mul_of_nonneg_left hsum hC
    linarith
  have h3 : T ≤ (1 + η) * U + (2 * Real.pi) * C * δ * T := by
    have := mul_le_mul_of_nonneg_left h2 hπ.le
    have e2 : 2 * Real.pi * ((2 * Real.pi)⁻¹ * T) = T := by field_simp
    have e3 : 2 * Real.pi * ((1 + η) * ((2 * Real.pi)⁻¹ * U) + C * (δ * T)) =
        (1 + η) * U + (2 * Real.pi) * C * δ * T := by field_simp
    linarith
  have e4 : (1 - C * (2 * Real.pi) * δ) * T = T - (2 * Real.pi) * C * δ * T := by ring
  rw [e4]; linarith

end Riemann.Deep.LeanNBRate



