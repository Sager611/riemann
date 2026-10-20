import Riemann.Deep.PNTError.Basic
import Riemann.Deep.WeilEF.Z23.FromPNTPlus.StrongPNTPrefix
import Riemann.Deep.Convexity.Partial

/-!
# Power-saving PNT, part 2: ζ'/ζ ≪ log|t| to the right of a zero-free half-plane

If ζ has no zero in Re s > σ₀ (1/2 ≤ σ₀ < 1), then for every η > 0 there is C with
‖ζ'(s)/ζ(s)‖ ≤ C log|t| for Re s ≥ σ₀ + η, |t| ≥ 3 (`logDeriv_bound`), and for |t| ≥ 2
(`logDeriv_bound_two`).

Proof (Titchmarsh §14.2 style, Borel–Carathéodory): on the disc |z| ≤ R = 2 - σ₀ - η/2 put
B(z) = ζ(2 + it + z), which has no zero there; J = log B - log B(0) is analytic with
Re J = log|B| - log|B(0)| ≤ log(15(|t|+2)) + log 2, using the crude bound ‖ζ(w)‖ ≤ 15|Im w|
(`norm_zeta_le_crude`, from Convexity/Partial.lean) and ‖ζ(2+it)‖ ≥ 1/2 (`norm_zeta_ge_half`,
from L(1,s)·L(μ,s) = 1). Borel–Carathéodory for J' = ζ'/ζ (`BorelCaratheodoryDeriv`, ported PNT+
code) on |z| ≤ r = 2 - σ₀ - η gives the bound for σ₀ + η ≤ Re s ≤ 2. For Re s ≥ 2,
‖ζ'/ζ(s)‖ ≤ Σ Λ(n)/n² (`norm_logDeriv_le_K2`).
-/

open Complex Set Filter Topology
open scoped Real ComplexConjugate LSeries.notation ArithmeticFunction.vonMangoldt
  ArithmeticFunction.Moebius

namespace Fayed.PNTError

noncomputable section

/-- Crude growth bound: ‖ζ(w)‖ ≤ 15|Im w| for 1/2 ≤ Re w ≤ 4, |Im w| ≥ 1. -/
lemma norm_zeta_le_crude {w : ℂ} (h1 : 1 / 2 ≤ w.re) (h4 : w.re ≤ 4) (him : 1 ≤ |w.im|) :
    ‖riemannZeta w‖ ≤ 15 * |w.im| := by
  have key : ∀ v : ℂ, 1 / 2 ≤ v.re → v.re ≤ 4 → 1 ≤ v.im → ‖riemannZeta v‖ ≤ 15 * v.im := by
    intro v hv1 hv4 hv
    have hvim : 0 < v.im := by linarith
    have hvim' : v.im ≠ 0 := hvim.ne'
    have hU : v ∈ Fayed.Convexity.U := ⟨show 0 < v.re by linarith, show 0 < v.im from hvim⟩
    have h := Fayed.Convexity.norm_zeta_le_on_U hU
    have hn : ‖v‖ ≤ 5 * v.im := by
      calc ‖v‖ ≤ |v.re| + |v.im| := Complex.norm_le_abs_re_add_abs_im v
        _ = v.re + v.im := by rw [abs_of_pos (by linarith), abs_of_pos hvim]
        _ ≤ 5 * v.im := by linarith
    have hn1 : v.im ≤ ‖v - 1‖ := by
      have := Complex.abs_im_le_norm (v - 1)
      rw [sub_im, one_im, sub_zero] at this
      linarith [le_abs_self v.im]
    calc ‖riemannZeta v‖ ≤ ‖v‖ / ‖v - 1‖ + ‖v‖ * (1 / v.re) := h
      _ ≤ 5 * v.im / v.im + 5 * v.im * 2 := by
          apply add_le_add
          · exact div_le_div₀ (by positivity) hn hvim hn1
          · apply mul_le_mul hn _ (by positivity) (by positivity)
            rw [div_le_iff₀ (by linarith)]
            linarith
      _ = 5 + 10 * v.im := by rw [mul_div_cancel_right₀ _ hvim']; ring
      _ ≤ 15 * v.im := by linarith
  rcases le_or_gt 0 w.im with hpos | hneg
  · rw [abs_of_nonneg hpos] at him ⊢
    exact key w h1 h4 him
  · have hc : riemannZeta w = conj (riemannZeta (conj w)) := by
      rw [← riemannZeta_conj, Complex.conj_conj]
    rw [abs_of_neg hneg] at him ⊢
    rw [hc, Complex.norm_conj]
    have := key (conj w) (by simpa using h1) (by simpa using h4) (by simpa using him)
    simpa using this

/-- An L-series with coefficients of norm ≤ 1 is ≤ π²/6 in norm on Re s ≥ 2. -/
lemma norm_LSeries_le_of_norm_le_one {f : ℕ → ℂ} (hf : ∀ n, ‖f n‖ ≤ 1) {s : ℂ}
    (hs : 2 ≤ s.re) : ‖LSeries f s‖ ≤ π ^ 2 / 6 := by
  have hb : ∀ n : ℕ, ‖LSeries.term f s n‖ ≤ 1 / (n : ℝ) ^ 2 := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · simp
    · simp only [LSeries.norm_term_eq, hn, ↓reduceIte]
      have hn1 : (1 : ℝ) ≤ n := Nat.one_le_cast.mpr (Nat.one_le_iff_ne_zero.mpr hn)
      have hpow : (n : ℝ) ^ (2 : ℕ) ≤ (n : ℝ) ^ s.re := by
        rw [← Real.rpow_natCast]
        exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num; linarith)
      calc ‖f n‖ / (n : ℝ) ^ s.re ≤ 1 / (n : ℝ) ^ s.re :=
            div_le_div_of_nonneg_right (hf n) (by positivity)
        _ ≤ 1 / (n : ℝ) ^ 2 := one_div_le_one_div_of_le (by positivity) hpow
  have hsum : Summable (fun n : ℕ => ‖LSeries.term f s n‖) :=
    Summable.of_nonneg_of_le (fun _ => norm_nonneg _) hb
      (Real.summable_one_div_nat_pow.mpr (by norm_num))
  calc ‖LSeries f s‖ = ‖∑' n, LSeries.term f s n‖ := rfl
    _ ≤ ∑' n, ‖LSeries.term f s n‖ := norm_tsum_le_tsum_norm hsum
    _ ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ 2 := hsum.tsum_le_tsum hb hasSum_zeta_two.summable
    _ = π ^ 2 / 6 := hasSum_zeta_two.tsum_eq

/-- ‖ζ(s)‖ ≥ 1/2 on Re s ≥ 2 (from ζ(s) · Σ μ(n) n^{-s} = 1). -/
lemma norm_zeta_ge_half {s : ℂ} (hs : 2 ≤ s.re) : 1 / 2 ≤ ‖riemannZeta s‖ := by
  have hs1 : 1 < s.re := by linarith
  have h1 := LSeries_one_mul_Lseries_moebius hs1
  rw [LSeries_one_eq_riemannZeta hs1] at h1
  have hb : ‖LSeries ↗μ s‖ ≤ π ^ 2 / 6 := by
    refine norm_LSeries_le_of_norm_le_one (fun n => ?_) hs
    rw [Complex.norm_intCast]
    exact_mod_cast ArithmeticFunction.abs_moebius_le_one
  have hπ : π ^ 2 / 6 < 2 := by nlinarith [Real.pi_lt_d2, Real.pi_pos]
  have hn := congrArg norm h1
  rw [norm_mul, norm_one] at hn
  nlinarith [norm_nonneg (riemannZeta s), norm_nonneg (LSeries ↗μ s)]

/-- Σ Λ(n)/n², the bound for ζ'/ζ on Re s ≥ 2. -/
def K2 : ℝ := ∑' n : ℕ, ‖LSeries.term ↗Λ 2 n‖

lemma K2_nonneg : 0 ≤ K2 := tsum_nonneg (fun _ => norm_nonneg _)

lemma norm_logDeriv_le_K2 {s : ℂ} (hs : 2 ≤ s.re) :
    ‖deriv riemannZeta s / riemannZeta s‖ ≤ K2 := by
  have hs1 : 1 < s.re := by linarith
  have e : ‖deriv riemannZeta s / riemannZeta s‖ = ‖LSeries ↗Λ s‖ := by
    rw [ArithmeticFunction.LSeries_vonMangoldt_eq_deriv_riemannZeta_div hs1, neg_div, norm_neg]
  rw [e]
  have hsum2 : Summable (fun n => ‖LSeries.term ↗Λ 2 n‖) :=
    summable_norm_iff.mpr (ArithmeticFunction.LSeriesSummable_vonMangoldt (by norm_num))
  have hle : ∀ n, ‖LSeries.term ↗Λ s n‖ ≤ ‖LSeries.term ↗Λ 2 n‖ :=
    LSeries.norm_term_le_of_re_le_re _ (by simpa using hs)
  have hsums : Summable (fun n => ‖LSeries.term ↗Λ s n‖) :=
    Summable.of_nonneg_of_le (fun _ => norm_nonneg _) hle hsum2
  calc ‖LSeries ↗Λ s‖ = ‖∑' n, LSeries.term ↗Λ s n‖ := rfl
    _ ≤ ∑' n, ‖LSeries.term ↗Λ s n‖ := norm_tsum_le_tsum_norm hsums
    _ ≤ K2 := hsums.tsum_le_tsum hle hsum2

/-- Borel–Carathéodory: ‖ζ'/ζ(s)‖ ≤ 1728/η³ · log|t| for σ₀ + η ≤ Re s ≤ 2, |t| ≥ 3. -/
theorem logDeriv_bound_strip {σ₀ : ℝ} (hZ : ZeroFree σ₀) (h12 : 1 / 2 ≤ σ₀) (h1 : σ₀ < 1)
    {η : ℝ} (hη0 : 0 < η) (hη : η ≤ 1 / 2) {s : ℂ} (hs1 : σ₀ + η ≤ s.re) (hs2 : s.re ≤ 2)
    (ht : 3 ≤ |s.im|) :
    ‖deriv riemannZeta s / riemannZeta s‖ ≤ 1728 / η ^ 3 * Real.log |s.im| := by
  obtain ⟨c, hc⟩ : ∃ c : ℂ, c = (2 : ℂ) + (s.im : ℂ) * I := ⟨_, rfl⟩
  have hcre : c.re = 2 := by rw [hc]; simp
  have hcim : c.im = s.im := by rw [hc]; simp
  obtain ⟨R, hR⟩ : ∃ R : ℝ, R = 2 - σ₀ - η / 2 := ⟨_, rfl⟩
  obtain ⟨r, hr⟩ : ∃ r : ℝ, r = 2 - σ₀ - η := ⟨_, rfl⟩
  have hr0 : 0 < r := by rw [hr]; linarith
  have hrR : r < R := by rw [hr, hR]; linarith
  have hR32 : R ≤ 3 / 2 := by rw [hR]; linarith
  have hR0 : 0 < R := by linarith
  have hgeom : ∀ z ∈ Metric.closedBall (0 : ℂ) R,
      σ₀ + η / 2 ≤ (c + z).re ∧ (c + z).re ≤ 4 ∧ |s.im| - 3 / 2 ≤ |(c + z).im| ∧
        |(c + z).im| ≤ |s.im| + 3 / 2 := by
    intro z hz
    rw [mem_closedBall_zero_iff] at hz
    have hre := Complex.abs_re_le_norm z
    have him := Complex.abs_im_le_norm z
    rw [add_re, add_im, hcre, hcim]
    have hre' := abs_le.mp (hre.trans hz)
    have him' : |z.im| ≤ 3 / 2 := (him.trans hz).trans hR32
    refine ⟨by linarith [hre'.1], by linarith [hre'.2], ?_, ?_⟩
    · have h := abs_add_le (s.im + z.im) (-z.im)
      rw [add_neg_cancel_right, abs_neg] at h
      linarith
    · have h := abs_add_le s.im z.im
      linarith
  have hne1 : ∀ z ∈ Metric.closedBall (0 : ℂ) R, c + z ≠ 1 := by
    intro z hz h
    have h' := (hgeom z hz).2.2.1
    rw [h, one_im, abs_zero] at h'
    linarith
  have hanal : AnalyticOnNhd ℂ (fun z => riemannZeta (c + z)) (Metric.closedBall 0 R) := by
    intro z hz
    exact (analyticAt_zeta (hne1 z hz)).comp (analyticAt_const.add analyticAt_id)
  have hB0 : ∀ z ∈ Metric.closedBall (0 : ℂ) R, riemannZeta (c + z) ≠ 0 := fun z hz =>
    hZ _ (by linarith [(hgeom z hz).1])
  obtain ⟨J, hJan, hJ0, hJd, hJre⟩ := LogOfAnalyticFunction hr0 hrR hanal hB0
  obtain ⟨M, hM⟩ : ∃ M : ℝ, M = Real.log (15 * (|s.im| + 2)) + Real.log 2 := ⟨_, rfl⟩
  have hlow : 1 / 2 ≤ ‖riemannZeta c‖ := norm_zeta_ge_half (by rw [hcre])
  have hM0 : 0 < M := by
    have l1 : 0 < Real.log (15 * (|s.im| + 2)) := Real.log_pos (by linarith [abs_nonneg s.im])
    have l2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
    rw [hM]; linarith
  have hJM : ∀ z ∈ Metric.ball (0 : ℂ) R, (J z).re ≤ M := by
    intro z hz
    have hz' := Metric.ball_subset_closedBall hz
    rw [← hJre z hz]
    obtain ⟨g1, g2, g3, g4⟩ := hgeom z hz'
    have hBz : 0 < ‖riemannZeta (c + z)‖ := norm_pos_iff.mpr (hB0 z hz')
    have hup : ‖riemannZeta (c + z)‖ ≤ 15 * (|s.im| + 2) := by
      have := norm_zeta_le_crude (w := c + z) (by linarith) g2 (by linarith)
      linarith
    have l1 := Real.log_le_log hBz hup
    have l2 : Real.log (1 / 2) ≤ Real.log ‖riemannZeta c‖ := Real.log_le_log (by norm_num) hlow
    rw [one_div, Real.log_inv] at l2
    show Real.log ‖riemannZeta (c + z)‖ - Real.log ‖riemannZeta (c + 0)‖ ≤ M
    rw [add_zero, hM]
    linarith
  have hsc : s - c ∈ Metric.closedBall (0 : ℂ) r := by
    rw [mem_closedBall_zero_iff]
    have e : s - c = ((s.re - 2 : ℝ) : ℂ) := by
      apply Complex.ext <;> simp [hcre, hcim]
    rw [e, Complex.norm_real, Real.norm_eq_abs, abs_of_nonpos (by linarith), hr]
    linarith
  have hbc := BorelCaratheodoryDeriv hM0 hr0 hrR hJan.analyticOn hJ0 hJM hsc
  rw [hJd _ hsc] at hbc
  have hd : deriv (fun z => riemannZeta (c + z)) (s - c) = deriv riemannZeta s := by
    rw [deriv_comp_const_add]
    congr 1
    ring
  have hBs : riemannZeta (c + (s - c)) = riemannZeta s := by
    congr 1
    ring
  simp only [hd, hBs] at hbc
  refine hbc.trans ?_
  have hRr : R - r = η / 2 := by rw [hR, hr]; ring
  rw [hRr]
  have hlogt : Real.log 3 ≤ Real.log |s.im| := Real.log_le_log (by norm_num) ht
  have hlog3 : 0 < Real.log 3 := Real.log_pos (by norm_num)
  have hM6 : M ≤ 6 * Real.log |s.im| := by
    have ht0 : 0 < |s.im| := by linarith
    have h81 : (81 : ℝ) ≤ |s.im| ^ 4 := by
      have := pow_le_pow_left₀ (by norm_num : (0 : ℝ) ≤ 3) ht 4
      norm_num at this
      linarith
    have hA : 15 * (|s.im| + 2) ≤ |s.im| ^ 5 := by
      have e5 : |s.im| ^ 5 = |s.im| * |s.im| ^ 4 := by ring
      rw [e5]
      nlinarith
    have l1 : Real.log (15 * (|s.im| + 2)) ≤ 5 * Real.log |s.im| := by
      calc Real.log (15 * (|s.im| + 2)) ≤ Real.log (|s.im| ^ 5) :=
            Real.log_le_log (by positivity) hA
        _ = 5 * Real.log |s.im| := by rw [Real.log_pow]; norm_num
    have l2 : Real.log 2 ≤ Real.log |s.im| := Real.log_le_log (by norm_num) (by linarith)
    rw [hM]
    linarith
  have hR2 : R ^ 2 ≤ 9 / 4 := by nlinarith
  have hη0' : η ≠ 0 := hη0.ne'
  calc 16 * M * R ^ 2 / (η / 2) ^ 3 = 128 * (M * R ^ 2) / η ^ 3 := by
        field_simp
        ring
    _ ≤ 128 * (6 * Real.log |s.im| * (9 / 4)) / η ^ 3 := by
        apply div_le_div_of_nonneg_right _ (by positivity)
        apply mul_le_mul_of_nonneg_left _ (by norm_num)
        exact mul_le_mul hM6 hR2 (by positivity) (by linarith)
    _ = 1728 / η ^ 3 * Real.log |s.im| := by ring

/-- **Rung 1.** If ζ has no zero in Re s > σ₀ (1/2 ≤ σ₀ < 1), then ζ'/ζ ≪_η log|t| on
Re s ≥ σ₀ + η, |t| ≥ 3. -/
theorem logDeriv_bound {σ₀ : ℝ} (hZ : ZeroFree σ₀) (h12 : 1 / 2 ≤ σ₀) (h1 : σ₀ < 1)
    {η : ℝ} (hη0 : 0 < η) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ s : ℂ, σ₀ + η ≤ s.re → 3 ≤ |s.im| →
      ‖deriv riemannZeta s / riemannZeta s‖ ≤ C * Real.log |s.im| := by
  obtain ⟨η', hη'⟩ : ∃ η' : ℝ, η' = min η (1 / 2) := ⟨_, rfl⟩
  have hη'0 : 0 < η' := by rw [hη']; exact lt_min hη0 (by norm_num)
  have hη'1 : η' ≤ 1 / 2 := by rw [hη']; exact min_le_right _ _
  have hη'η : η' ≤ η := by rw [hη']; exact min_le_left _ _
  refine ⟨max (1728 / η' ^ 3) (K2 / Real.log 3),
    le_max_of_le_left (div_nonneg (by norm_num) (pow_pos hη'0 3).le), ?_⟩
  intro s hs ht
  have hlog3 : 0 < Real.log 3 := Real.log_pos (by norm_num)
  have hlogt : Real.log 3 ≤ Real.log |s.im| := Real.log_le_log (by norm_num) ht
  have hlt0 : 0 ≤ Real.log |s.im| := by linarith
  rcases le_or_gt s.re 2 with hs2 | hs2
  · calc ‖deriv riemannZeta s / riemannZeta s‖ ≤ 1728 / η' ^ 3 * Real.log |s.im| :=
          logDeriv_bound_strip hZ h12 h1 hη'0 hη'1 (by linarith) hs2 ht
      _ ≤ _ := mul_le_mul_of_nonneg_right (le_max_left _ _) hlt0
  · calc ‖deriv riemannZeta s / riemannZeta s‖ ≤ K2 := norm_logDeriv_le_K2 hs2.le
      _ = K2 / Real.log 3 * Real.log 3 := by field_simp
      _ ≤ K2 / Real.log 3 * Real.log |s.im| :=
          mul_le_mul_of_nonneg_left hlogt (div_nonneg K2_nonneg hlog3.le)
      _ ≤ _ := mul_le_mul_of_nonneg_right (le_max_right _ _) hlt0

/-- ζ'/ζ is bounded on compact subsets of the zero-free region that avoid s = 1. -/
lemma logDeriv_bdd_compact {σ₀ : ℝ} (hZ : ZeroFree σ₀) {K : Set ℂ} (hK : IsCompact K)
    (hK1 : ∀ s ∈ K, σ₀ < s.re ∧ s ≠ 1) :
    ∃ B : ℝ, ∀ s ∈ K, ‖deriv riemannZeta s / riemannZeta s‖ ≤ B := by
  have hc : ContinuousOn (fun s => deriv riemannZeta s / riemannZeta s) K := by
    intro s hs
    obtain ⟨h1, h2⟩ := hK1 s hs
    exact ((differentiableAt_deriv_zeta h2).continuousAt.div
      (differentiableAt_riemannZeta h2).continuousAt (hZ s h1)).continuousWithinAt
  exact hK.exists_bound_of_continuousOn hc

/-- **Rung 1, as stated in the brief**: ζ'/ζ ≪_η log|t| on Re s ≥ σ₀ + η, |t| ≥ 2. -/
theorem logDeriv_bound_two {σ₀ : ℝ} (hZ : ZeroFree σ₀) (h12 : 1 / 2 ≤ σ₀) (h1 : σ₀ < 1)
    {η : ℝ} (hη0 : 0 < η) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ s : ℂ, σ₀ + η ≤ s.re → 2 ≤ |s.im| →
      ‖deriv riemannZeta s / riemannZeta s‖ ≤ C * Real.log |s.im| := by
  obtain ⟨C₁, hC₁0, hC₁⟩ := logDeriv_bound hZ h12 h1 hη0
  set K : Set ℂ := {s | σ₀ + η ≤ s.re ∧ s.re ≤ 2 ∧ 2 ≤ |s.im| ∧ |s.im| ≤ 3} with hK
  have hKc : IsCompact K := by
    apply Metric.isCompact_of_isClosed_isBounded
    · refine (isClosed_le continuous_const continuous_re).inter
        ((isClosed_le continuous_re continuous_const).inter
          ((isClosed_le continuous_const (continuous_abs.comp continuous_im)).inter
            (isClosed_le (continuous_abs.comp continuous_im) continuous_const)))
    · refine (Metric.isBounded_closedBall (x := (0 : ℂ)) (r := 6)).subset (fun s hs => ?_)
      obtain ⟨a1, a2, _, a4⟩ := hs
      rw [mem_closedBall_zero_iff]
      have hre : |s.re| ≤ 2 := abs_le.mpr ⟨by linarith, a2⟩
      linarith [Complex.norm_le_abs_re_add_abs_im s]
  obtain ⟨B, hB⟩ := logDeriv_bdd_compact hZ hKc (fun s hs => by
    obtain ⟨a1, _, a3, _⟩ := hs
    refine ⟨by linarith, fun h => ?_⟩
    rw [h, one_im, abs_zero] at a3
    linarith)
  have hlog2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hD0 : 0 ≤ max B 0 + K2 := add_nonneg (le_max_right _ _) K2_nonneg
  refine ⟨max C₁ ((max B 0 + K2) / Real.log 2), le_max_of_le_left hC₁0, ?_⟩
  intro s hs ht
  have hlogt : Real.log 2 ≤ Real.log |s.im| := Real.log_le_log (by norm_num) ht
  have hlt0 : 0 ≤ Real.log |s.im| := by linarith
  by_cases h3 : 3 ≤ |s.im|
  · exact (hC₁ s hs h3).trans (mul_le_mul_of_nonneg_right (le_max_left _ _) hlt0)
  · rw [not_le] at h3
    have hD : ‖deriv riemannZeta s / riemannZeta s‖ ≤ max B 0 + K2 := by
      by_cases hs2 : s.re ≤ 2
      · have hb := hB s ⟨hs, hs2, ht, h3.le⟩
        linarith [le_max_left B 0, K2_nonneg]
      · rw [not_le] at hs2
        have hK2 := norm_logDeriv_le_K2 hs2.le
        linarith [le_max_right B 0]
    calc ‖deriv riemannZeta s / riemannZeta s‖ ≤ max B 0 + K2 := hD
      _ = (max B 0 + K2) / Real.log 2 * Real.log 2 := by field_simp
      _ ≤ (max B 0 + K2) / Real.log 2 * Real.log |s.im| :=
          mul_le_mul_of_nonneg_left hlogt (div_nonneg hD0 hlog2.le)
      _ ≤ _ := mul_le_mul_of_nonneg_right (le_max_right _ _) hlt0

end

end Fayed.PNTError
