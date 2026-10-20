/-
Riemann/Deep/Yoshida/Hybrid.lean — Weil positivity from a one-variable cell certificate.

Below `log 2`: `Re Q(f) = (1/2π) ∫ |F(r)|² W(r) dr` with `W = 2π μ + J` (`J` = pole kernel,
`PoleKernel.lean`). Layer cake with `|F|² ≤ 2aN`, `∫|F|² = 2πN` (`N = ∫|f|²`): if
`W(r) ≥ s − t(r)` for a step function `t ≥ 0` (cells `lo ≤ |r| < hi`, height `τ`), then
`∫ |F|² W ≥ 2πN s − 2aN ∫ t`, so `Q ≥ 0` once `2a ∫ t ≤ 2π s`.
-/
import Riemann.Deep.Yoshida.PoleKernel
import Riemann.Deep.Yoshida.MuBounds
import Riemann.Deep.Yoshida.Taylor

noncomputable section

open Complex MeasureTheory Set Filter
open scoped ComplexConjugate

namespace Fayed.Yoshida

open Riemann.Deep.Weil Riemann.Deep.WeilEF Zeta23 Zeta23.EF

/-- The Fourier-side weight `W(r) = 2π μ(r) + J(r)`. -/
def W (a r : ℝ) : ℝ := 2 * Real.pi * mu r + J a r

lemma abs_J_le {a : ℝ} (ha : 0 ≤ a) (r : ℝ) : |J a r| ≤ 8 * a * Real.cosh a := by
  rw [← integral_J]
  have hb : ∀ x ∈ Set.uIoc (-(2 * a)) (2 * a),
      ‖2 * Real.cosh (x / 2) * Real.cos (r * x)‖ ≤ 2 * Real.cosh a := by
    intro x hx
    rw [Set.uIoc_of_le (by linarith)] at hx
    have hxa : |x / 2| ≤ a := by rw [abs_le]; constructor <;> linarith [hx.1, hx.2]
    have hch : Real.cosh (x / 2) ≤ Real.cosh a := by
      rw [← Real.cosh_abs (x / 2)]
      exact Real.cosh_le_cosh.2 (by rw [abs_abs, abs_of_nonneg ha]; exact hxa)
    rw [Real.norm_eq_abs, abs_mul, abs_of_pos (by positivity : 0 < 2 * Real.cosh (x / 2))]
    have := Real.abs_cos_le_one (r * x)
    nlinarith [Real.cosh_pos (x / 2)]
  have h := intervalIntegral.norm_integral_le_of_norm_le_const hb
  rw [Real.norm_eq_abs] at h
  refine h.trans (le_of_eq ?_)
  rw [show 2 * a - -(2 * a) = 4 * a by ring, abs_of_nonneg (by linarith)]
  ring

/-- **Prime-free Fourier form.** Below `log 2`: `Re Q(f) = (1/2π) ∫ |F|² W`. -/
theorem weilQ_re_eq_W {a : ℝ} (ha : 0 ≤ a) (hlog : 2 * a < Real.log 2) {f : ℝ → ℂ}
    (hf : TestFn a f) :
    Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * W a r) ∧
    (weilQ f).re = 1 / (2 * Real.pi) * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W a r := by
  have hk := hf.autocorr_contDiff
  have hkc := hf.autocorr_hasCompactSupport
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  have hφi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2) := by
    refine (integrable_paperFT_ofReal hFk).re.congr (Eventually.of_forall fun r => ?_)
    simp only [paperFT_autocorr_real hf, RCLike.re_to_complex, Complex.ofReal_re]
  obtain ⟨hμi, harch⟩ := arch_re hf
  have hJi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * J a r) :=
    hφi.mul_bdd (continuous_J a).aestronglyMeasurable
      (Eventually.of_forall fun r => by rw [Real.norm_eq_abs]; exact abs_J_le ha r)
  have hμi' : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r)) := by
    refine (hμi.const_mul (2 * Real.pi)).congr (Eventually.of_forall fun r => ?_)
    simp only
    ring
  have hsum : (fun r : ℝ => ‖paperFT f r‖ ^ 2 * W a r)
      = fun r : ℝ => ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r) + ‖paperFT f r‖ ^ 2 * J a r := by
    funext r
    unfold W
    ring
  refine ⟨by rw [hsum]; exact hμi'.add hJi, ?_⟩
  rw [weilQ_eq_fix]
  unfold weilQFix
  simp only [Complex.add_re, Complex.sub_re]
  rw [primeTerm_eq_zero hf.2 hlog, harch, pole_re_eq_J ha hf, Complex.zero_re, sub_zero, hsum,
    integral_add hμi' hJi]
  have hπ := Real.pi_pos
  have e : ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * (2 * Real.pi * mu r)
      = 2 * Real.pi * ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * mu r := by
    rw [← integral_const_mul]
    congr 1
    funext r
    ring
  rw [e]
  field_simp
  ring

/-! ### Layer cake with a step function -/

lemma bath {φ w t : ℝ → ℝ} {c s T : ℝ} (hφ0 : ∀ r, 0 ≤ φ r) (hφc : ∀ r, φ r ≤ c)
    (hφi : Integrable φ) (hφw : Integrable (fun r => φ r * w r)) (hT : ∫ r, φ r = T)
    (ht0 : ∀ r, 0 ≤ t r) (hti : Integrable t) (hw : ∀ r, s - t r ≤ w r) :
    s * T - c * ∫ r, t r ≤ ∫ r, φ r * w r := by
  have hpt : ∀ r, s * φ r - c * t r ≤ φ r * w r := by
    intro r
    have h1 := mul_le_mul_of_nonneg_left (hw r) (hφ0 r)
    have h2 := mul_le_mul_of_nonneg_right (hφc r) (ht0 r)
    nlinarith
  have hint : Integrable (fun r => s * φ r - c * t r) := (hφi.const_mul s).sub (hti.const_mul c)
  have h := integral_mono hint hφw hpt
  rw [integral_sub (hφi.const_mul s) (hti.const_mul c), integral_const_mul, integral_const_mul,
    hT] at h
  exact h

/-- One cell `(lo, hi, τ)`: height `τ` on `lo ≤ |r| < hi` (both signs; `r = 0` may count twice). -/
def cellFn (c : ℝ × ℝ × ℝ) (r : ℝ) : ℝ :=
  (Ico c.1 c.2.1).indicator (fun _ => c.2.2) r + (Ioc (-c.2.1) (-c.1)).indicator (fun _ => c.2.2) r

def stepFn : List (ℝ × ℝ × ℝ) → ℝ → ℝ
  | [], _ => 0
  | c :: cs, r => cellFn c r + stepFn cs r

/-- Total integral of the step function. -/
def cost : List (ℝ × ℝ × ℝ) → ℝ
  | [] => 0
  | c :: cs => c.2.2 * (2 * (c.2.1 - c.1)) + cost cs

/-- Consecutive cells from `x` to `R`. -/
def Chain : ℝ → List (ℝ × ℝ × ℝ) → ℝ → Prop
  | x, [], R => x = R
  | x, c :: cs, R => c.1 = x ∧ Chain c.2.1 cs R

/-- A certified cell: `W ≥ s − τ` on `lo ≤ |r| < hi`. -/
def CellOK (a s : ℝ) (c : ℝ × ℝ × ℝ) : Prop :=
  0 ≤ c.1 ∧ c.1 ≤ c.2.1 ∧ 0 ≤ c.2.2 ∧ ∀ r : ℝ, c.1 ≤ |r| → |r| < c.2.1 → s - c.2.2 ≤ W a r

lemma cellFn_nonneg {c : ℝ × ℝ × ℝ} (h : 0 ≤ c.2.2) (r : ℝ) : 0 ≤ cellFn c r := by
  unfold cellFn
  have h1 : 0 ≤ (Ico c.1 c.2.1).indicator (fun _ => c.2.2) r := Set.indicator_nonneg (fun _ _ => h) r
  have h2 : 0 ≤ (Ioc (-c.2.1) (-c.1)).indicator (fun _ => c.2.2) r :=
    Set.indicator_nonneg (fun _ _ => h) r
  linarith

lemma integrable_cellFn (c : ℝ × ℝ × ℝ) : Integrable (cellFn c) := by
  unfold cellFn
  refine Integrable.add ?_ ?_
  · exact (integrableOn_const (C := c.2.2) (by simp)).integrable_indicator measurableSet_Ico
  · exact (integrableOn_const (C := c.2.2) (by simp)).integrable_indicator measurableSet_Ioc

lemma integral_cellFn {c : ℝ × ℝ × ℝ} (h : c.1 ≤ c.2.1) :
    ∫ r, cellFn c r = c.2.2 * (2 * (c.2.1 - c.1)) := by
  unfold cellFn
  rw [integral_add ((integrableOn_const (C := c.2.2) (by simp)).integrable_indicator measurableSet_Ico)
    ((integrableOn_const (C := c.2.2) (by simp)).integrable_indicator measurableSet_Ioc),
    integral_indicator_const _ measurableSet_Ico, integral_indicator_const _ measurableSet_Ioc,
    Real.volume_real_Ico_of_le h, Real.volume_real_Ioc_of_le (by linarith), smul_eq_mul, smul_eq_mul]
  ring

lemma le_cellFn {c : ℝ × ℝ × ℝ} (hτ : 0 ≤ c.2.2) {r : ℝ} (h1 : c.1 ≤ |r|) (h2 : |r| < c.2.1) :
    c.2.2 ≤ cellFn c r := by
  unfold cellFn
  have i1 : 0 ≤ (Ico c.1 c.2.1).indicator (fun _ => c.2.2) r := Set.indicator_nonneg (fun _ _ => hτ) r
  have i2 : 0 ≤ (Ioc (-c.2.1) (-c.1)).indicator (fun _ => c.2.2) r :=
    Set.indicator_nonneg (fun _ _ => hτ) r
  rcases le_total 0 r with hr | hr
  · rw [abs_of_nonneg hr] at h1 h2
    rw [Set.indicator_of_mem (show r ∈ Ico c.1 c.2.1 from ⟨h1, h2⟩)]
    linarith
  · rw [abs_of_nonpos hr] at h1 h2
    rw [Set.indicator_of_mem (show r ∈ Ioc (-c.2.1) (-c.1) from ⟨by linarith, by linarith⟩)]
    linarith

lemma stepFn_nonneg : ∀ (cs : List (ℝ × ℝ × ℝ)), (∀ c ∈ cs, 0 ≤ c.2.2) → ∀ r, 0 ≤ stepFn cs r
  | [], _, _ => le_refl _
  | c :: cs, h, r => by
    simp only [List.mem_cons, forall_eq_or_imp] at h
    exact add_nonneg (cellFn_nonneg h.1 r) (stepFn_nonneg cs h.2 r)

lemma integrable_stepFn : ∀ (cs : List (ℝ × ℝ × ℝ)), Integrable (stepFn cs)
  | [] => by
    have : stepFn [] = fun _ => (0 : ℝ) := rfl
    rw [this]
    exact integrable_zero _ _ _
  | c :: cs => (integrable_cellFn c).add (integrable_stepFn cs)

lemma integral_stepFn : ∀ (cs : List (ℝ × ℝ × ℝ)), (∀ c ∈ cs, c.1 ≤ c.2.1) →
    ∫ r, stepFn cs r = cost cs
  | [], _ => by simp [stepFn, cost]
  | c :: cs, h => by
    simp only [List.mem_cons, forall_eq_or_imp] at h
    have e : (fun r => stepFn (c :: cs) r) = fun r => cellFn c r + stepFn cs r := rfl
    rw [e, integral_add (integrable_cellFn c) (integrable_stepFn cs), integral_cellFn h.1,
      integral_stepFn cs h.2]
    rfl

lemma le_stepFn : ∀ (cs : List (ℝ × ℝ × ℝ)), (∀ c ∈ cs, 0 ≤ c.2.2) → ∀ c ∈ cs, ∀ r : ℝ,
    c.1 ≤ |r| → |r| < c.2.1 → c.2.2 ≤ stepFn cs r
  | [], _, c, hc, _, _, _ => by simp at hc
  | d :: cs, h, c, hc, r, h1, h2 => by
    simp only [List.mem_cons, forall_eq_or_imp] at h
    have hs : stepFn (d :: cs) r = cellFn d r + stepFn cs r := rfl
    rw [hs]
    rcases List.mem_cons.1 hc with rfl | hc'
    · have := stepFn_nonneg cs h.2 r
      linarith [le_cellFn h.1 h1 h2]
    · have := cellFn_nonneg h.1 r
      linarith [le_stepFn cs h.2 c hc' r h1 h2]

lemma chain_cover : ∀ (cs : List (ℝ × ℝ × ℝ)) (x R : ℝ), Chain x cs R → ∀ y, x ≤ y → y < R →
    ∃ c ∈ cs, c.1 ≤ y ∧ y < c.2.1
  | [], x, R, h, y, h1, h2 => by
    simp only [Chain] at h
    linarith
  | c :: cs, x, R, h, y, h1, h2 => by
    simp only [Chain] at h
    rcases lt_or_ge y c.2.1 with hy | hy
    · exact ⟨c, List.mem_cons_self .., by rw [h.1]; exact h1, hy⟩
    · obtain ⟨d, hd, hd'⟩ := chain_cover cs c.2.1 R h.2 y hy h2
      exact ⟨d, List.mem_cons_of_mem _ hd, hd'⟩

/-- **Certificate ⇒ Weil positivity** on the window `2a`. -/
theorem weilPos_of_cert {a s R : ℝ} (ha : 0 < a) (hlog : 2 * a < Real.log 2)
    (cs : List (ℝ × ℝ × ℝ)) (hchain : Chain 0 cs R) (hok : ∀ c ∈ cs, CellOK a s c)
    (htail : ∀ r : ℝ, R ≤ |r| → s ≤ W a r) (hnum : 2 * a * cost cs ≤ 2 * Real.pi * s) :
    WeilPos (2 * a) := by
  intro f hf
  have hf' : TestFn a f := by
    have e : 2 * a / 2 = a := by ring
    rwa [e] at hf
  obtain ⟨hφw, hre⟩ := weilQ_re_eq_W ha.le hlog hf'
  have hk := hf'.autocorr_contDiff
  have hkc := hf'.autocorr_hasCompactSupport
  have hFk := integrable_fourier_of_contDiff_two hk hkc
  set N := ∫ u, ‖f u‖ ^ 2 with hN
  have hN0 : 0 ≤ N := integral_nonneg fun u => by positivity
  have hφi : Integrable (fun r : ℝ => ‖paperFT f r‖ ^ 2) := by
    refine (integrable_paperFT_ofReal hFk).re.congr (Eventually.of_forall fun r => ?_)
    simp only [paperFT_autocorr_real hf', RCLike.re_to_complex, Complex.ofReal_re]
  have hφc : ∀ r : ℝ, ‖paperFT f r‖ ^ 2 ≤ 2 * a * N := fun r =>
    (pow_le_pow_left₀ (norm_nonneg _) (norm_paperFT_le f r) 2).trans (sq_integral_norm_le ha hf')
  have hτ : ∀ c ∈ cs, 0 ≤ c.2.2 := fun c hc => (hok c hc).2.2.1
  have hle : ∀ c ∈ cs, c.1 ≤ c.2.1 := fun c hc => (hok c hc).2.1
  have hw : ∀ r : ℝ, s - stepFn cs r ≤ W a r := by
    intro r
    rcases le_or_gt R |r| with hR | hR
    · have := stepFn_nonneg cs hτ r
      linarith [htail r hR]
    · obtain ⟨c, hc, h1, h2⟩ := chain_cover cs 0 R hchain |r| (abs_nonneg r) hR
      have := le_stepFn cs hτ c hc r h1 h2
      linarith [(hok c hc).2.2.2 r h1 h2]
  have hb := bath (φ := fun r : ℝ => ‖paperFT f r‖ ^ 2) (w := W a) (t := stepFn cs)
    (fun r => by positivity) hφc hφi hφw (plancherel hf') (stepFn_nonneg cs hτ)
    (integrable_stepFn cs) hw
  rw [integral_stepFn cs hle] at hb
  rw [hre]
  have hπ := Real.pi_pos
  have h2 : 0 ≤ s * (2 * Real.pi * N) - 2 * a * N * cost cs := by
    have := mul_le_mul_of_nonneg_left hnum hN0
    nlinarith
  have h3 : 0 ≤ ∫ r : ℝ, ‖paperFT f r‖ ^ 2 * W a r := le_trans h2 hb
  positivity

/-! ### Cell and tail lemmas used by generated certificates -/

/-- A cell from a `μ` bound at `lo`, a `J` bound at `p`, and the Lipschitz bound of `J`. -/
lemma cellOK_of {a s lo hi τ p d m j Λ : ℝ} (ha : 0 ≤ a) (hlo : 0 ≤ lo) (hlh : lo ≤ hi) (hτ : 0 ≤ τ)
    (hd1 : p - lo ≤ d) (hd2 : hi - p ≤ d) (hm : m ≤ 2 * Real.pi * mu lo) (hj : j ≤ J a p)
    (hΛ : 8 * a ^ 2 * Real.cosh a ≤ Λ) (hnum : s - τ ≤ m + j - Λ * d) :
    CellOK a s (lo, hi, τ) := by
  refine ⟨hlo, hlh, hτ, fun r h1 h2 => ?_⟩
  simp only at h1 h2
  have hmu : mu lo ≤ mu r := mu_abs_mono hlo h1
  have hJr : J a r = J a |r| := by
    rcases le_total 0 r with h | h
    · rw [abs_of_nonneg h]
    · rw [abs_of_nonpos h, J_even]
  have hL := abs_J_sub_le ha |r| p
  have hdist : |(|r| - p)| ≤ d := by
    rw [abs_le]; constructor <;> linarith
  have hd0 : 0 ≤ d := le_trans (abs_nonneg _) hdist
  have hc0 : 0 ≤ 8 * a ^ 2 * Real.cosh a := by positivity
  have k1 : 8 * a ^ 2 * Real.cosh a * |(|r| - p)| ≤ Λ * d :=
    mul_le_mul hΛ hdist (abs_nonneg _) (le_trans hc0 hΛ)
  have k2 := (abs_le.1 hL).1
  unfold W
  rw [hJr]
  have hπ := Real.pi_pos
  nlinarith

/-- `J` lower bound from lower bounds for `cos`, `sin` at `2ap`. -/
lemma J_ge_of {a p cl sl X Y : ℝ} (ha : 0 ≤ a) (hp : 0 ≤ p) (hc : cl ≤ Real.cos (2 * a * p))
    (hs : sl ≤ Real.sin (2 * a * p)) (hX : X ≤ Real.sinh a * cl) (hY : Y ≤ Real.cosh a * sl) :
    4 * (1 / 2 * X + p * Y) / (1 / 4 + p ^ 2) ≤ J a p := by
  unfold J
  have hsh : 0 ≤ Real.sinh a := Real.sinh_nonneg_iff.2 ha
  have hch : 0 ≤ Real.cosh a := (Real.cosh_pos a).le
  have k1 := mul_le_mul_of_nonneg_left hc hsh
  have k2 := mul_le_mul_of_nonneg_left hs hch
  have hden : 0 < 1 / 4 + p ^ 2 := by positivity
  apply div_le_div_of_nonneg_right _ hden.le
  have : p * Y ≤ p * (Real.cosh a * Real.sin (2 * a * p)) :=
    mul_le_mul_of_nonneg_left (by linarith) hp
  nlinarith

/-- The tail: `W(r) ≥ m − 4(sinh a/(2R²) + cosh a/R)` for `|r| ≥ R ≥ 1`. -/
lemma tail_of {a s R m sh ch : ℝ} (ha : 0 ≤ a) (hR : 1 ≤ R) (hm : m ≤ 2 * Real.pi * mu R)
    (hsh : Real.sinh a ≤ sh) (hch : Real.cosh a ≤ ch)
    (hnum : s ≤ m - 4 * (sh / (2 * R ^ 2) + ch / R)) : ∀ r : ℝ, R ≤ |r| → s ≤ W a r := by
  intro r hr
  have hmu : mu R ≤ mu r := mu_abs_mono (by linarith) hr
  have hJ := J_ge_tail ha hR hr
  have hR0 : 0 < R := by linarith
  have k : Real.sinh a / (2 * R ^ 2) + Real.cosh a / R ≤ sh / (2 * R ^ 2) + ch / R := by
    gcongr
  unfold W
  have hπ := Real.pi_pos
  nlinarith

end Fayed.Yoshida
