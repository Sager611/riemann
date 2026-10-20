/-
Riemann/Deep/WeilEF/Suff.lean — sufficiency half of Weil's criterion, sorry-free:
positivity of the Weil functional on every window implies RH.

Route (Weil/Bombieri, made elementary): fix a test function `f` and put
`D τ = Σ_ρ w_ρ e^{(ρ−1/2)τ}`, `w_ρ = m_ρ F(γ_ρ) conj F(conj γ_ρ)`.  The explicit formula (ported zeta23,
`Bridge.lean`) applied to `f + c·f(· − τ)` gives `Q = (1+|c|²) D 0 + c D τ + conj c · D (−τ)`, and the
zero reflection `ρ ↦ 1 − conj ρ` gives `D (−τ) = conj (D τ)`; positivity for all `c` forces
`‖D τ‖ ≤ Re D 0`.  A bounded `D` cannot carry weight at a zero with `Re ρ > 1/2` (`Core.lean`:
Laplace transform + identity theorem), while a bump times `e^{−i t₀ u}` gives `w_{ρ₀} ≠ 0`.
-/
import Riemann.Deep.WeilEF.Bridge
import Riemann.Deep.WeilEF.Core
import Riemann.Tree.R007
import Mathlib.Analysis.Calculus.BumpFunction.InnerProduct

noncomputable section

open Complex MeasureTheory Set Filter
open scoped ComplexConjugate

namespace Fayed.WeilSuff

open Riemann.Deep.Weil Riemann.Deep.WeilEF Zeta23 Zeta23.EF

/-! ### Zeros: countability, exponents, reflection -/

instance countable_carrier : Countable zetaZeroConfig.carrier := by
  have h : zetaZeroConfig.carrier ⊆
      ⋃ n : ℤ, zetaZeroConfig.carrier ∩ {ρ | (n : ℝ) < ρ.im ∧ ρ.im ≤ (n : ℝ) + 1} := by
    intro ρ hρ
    refine mem_iUnion.mpr ⟨⌈ρ.im⌉ - 1, hρ, ?_, ?_⟩
    · push_cast; linarith [Int.ceil_lt_add_one ρ.im]
    · push_cast; linarith [Int.le_ceil ρ.im]
  exact ((countable_iUnion fun n : ℤ =>
    (zetaZeroConfig.finite_window (n : ℝ) ((n : ℝ) + 1)).countable).mono h).to_subtype

lemma strip_of_mem {ρ : ℂ} (h : ρ ∈ zetaZeroConfig.carrier) : 0 < ρ.re ∧ ρ.re < 1 := h.2

/-- Exponent `ρ − 1/2 = i γ_ρ`. -/
def pe (ρ : zetaZeroConfig.carrier) : ℂ := (ρ : ℂ) - 1 / 2

lemma pe_re (ρ : zetaZeroConfig.carrier) : (pe ρ).re = (ρ : ℂ).re - 1 / 2 := by
  simp [pe]

lemma pe_im (ρ : zetaZeroConfig.carrier) : (pe ρ).im = (ρ : ℂ).im := by simp [pe]

lemma abs_pe_re_le (ρ : zetaZeroConfig.carrier) : |(pe ρ).re| ≤ 1 / 2 := by
  have := strip_of_mem ρ.2
  rw [pe_re, abs_le]; constructor <;> linarith

lemma I_mul_gammaOf (z : ℂ) : I * gammaOf z = z - 1 / 2 := by
  unfold gammaOf; field_simp

/-- Weight of a zero for the test function `f`. -/
def wt (f : ℝ → ℂ) (ρ : zetaZeroConfig.carrier) : ℂ :=
  (zetaZeroConfig.mult ρ : ℂ) * (paperFT f (gammaOf ρ) * conj (paperFT f (conj (gammaOf ρ))))

/-- The reflection `ρ ↦ 1 − conj ρ` of the zero set. -/
def reflEquiv : zetaZeroConfig.carrier ≃ zetaZeroConfig.carrier :=
  Function.Involutive.toPerm (fun ρ => ⟨reflect ρ, zetaZeroConfig.reflect_mem ρ ρ.2⟩)
    (fun ρ => Subtype.ext (by simp [reflect]))

lemma reflEquiv_val (ρ : zetaZeroConfig.carrier) : ((reflEquiv ρ : zetaZeroConfig.carrier) : ℂ)
    = reflect ρ := rfl

lemma gammaOf_reflect (z : ℂ) : gammaOf (reflect z) = conj (gammaOf z) := by
  unfold gammaOf reflect
  rw [map_div₀, Complex.conj_I, map_sub, map_div₀, map_one, map_ofNat]
  field_simp
  ring

lemma pe_refl (ρ : zetaZeroConfig.carrier) : pe (reflEquiv ρ) = -conj (pe ρ) := by
  simp only [pe, reflEquiv_val, reflect, map_sub, map_div₀, map_one, map_ofNat]
  ring

lemma wt_refl (f : ℝ → ℂ) (ρ : zetaZeroConfig.carrier) : wt f (reflEquiv ρ) = conj (wt f ρ) := by
  unfold wt
  rw [reflEquiv_val, gammaOf_reflect, zetaZeroConfig.mult_reflect ρ ρ.2]
  simp only [map_mul, map_natCast, Complex.conj_conj]
  ring

/-! ### Shifted combinations `f + c · f(· − τ)` -/

/-- `f + c · f(· − τ)`. -/
def shiftComb (f : ℝ → ℂ) (c : ℂ) (τ : ℝ) : ℝ → ℂ := fun u => f u + c * f (u - τ)

lemma testFn_shiftComb {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (c : ℂ) (τ : ℝ) :
    TestFn (a + |τ|) (shiftComb f c τ) := by
  have h1 := hf.1
  refine ⟨by unfold shiftComb; fun_prop, ?_⟩
  intro u hu
  simp only [Function.mem_support, shiftComb] at hu
  have h2 := le_abs_self τ
  have h3 := neg_abs_le τ
  by_cases h : f u = 0
  · rw [h, zero_add] at hu
    have hf2 : f (u - τ) ≠ 0 := fun h' => hu (by rw [h', mul_zero])
    have := hf.2 (Function.mem_support.mpr hf2)
    simp only [mem_Icc] at this ⊢
    constructor <;> linarith [this.1, this.2]
  · have := hf.2 (Function.mem_support.mpr h)
    simp only [mem_Icc] at this ⊢
    constructor <;> linarith [this.1, this.2, abs_nonneg τ]

lemma integrable_mul_cexp {g : ℝ → ℂ} (hg : Continuous g) (hgs : HasCompactSupport g) (z : ℂ) :
    Integrable (fun u : ℝ => g u * cexp (I * z * u)) :=
  (hg.mul (by fun_prop)).integrable_of_hasCompactSupport hgs.mul_right

lemma paperFT_shiftComb {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (c : ℂ) (τ : ℝ) (z : ℂ) :
    paperFT (shiftComb f c τ) z = paperFT f z * (1 + c * cexp (I * z * τ)) := by
  have hc := hf.continuous
  have hcs := hf.hasCompactSupport
  have hcs' : HasCompactSupport (fun u : ℝ => f (u - τ)) :=
    hcs.comp_homeomorph (Homeomorph.subRight τ)
  have hc' : Continuous (fun u : ℝ => f (u - τ)) := hc.comp (continuous_id.sub continuous_const)
  unfold paperFT shiftComb
  have e1 : (fun u : ℝ => (f u + c * f (u - τ)) * cexp (I * z * u))
      = fun u => f u * cexp (I * z * u) + c * (f (u - τ) * cexp (I * z * u)) := by
    funext u; ring
  rw [e1, integral_add (integrable_mul_cexp hc hcs z) ((integrable_mul_cexp hc' hcs' z).const_mul c),
    integral_const_mul]
  have e2 : ∫ u : ℝ, f (u - τ) * cexp (I * z * u)
      = (∫ u : ℝ, f u * cexp (I * z * u)) * cexp (I * z * τ) := by
    rw [← integral_mul_const,
      ← integral_sub_right_eq_self (fun u : ℝ => f u * cexp (I * z * u) * cexp (I * z * τ)) τ]
    congr 1; funext u
    rw [mul_assoc (f (u - τ)), ← Complex.exp_add]
    congr 2; push_cast; ring
  rw [e2]; ring

/-- Termwise expansion of the zero sum of `f + c·f(· − τ)`. -/
lemma term_eq {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (c : ℂ) (τ : ℝ)
    (ρ : zetaZeroConfig.carrier) :
    (zetaZeroConfig.mult ρ : ℂ) * (paperFT (shiftComb f c τ) (gammaOf ρ) *
        conj (paperFT (shiftComb f c τ) (conj (gammaOf ρ))))
      = (1 + c * conj c) * wt f ρ + c * (wt f ρ * cexp (pe ρ * τ))
        + conj c * (wt f ρ * cexp (pe ρ * ((-τ : ℝ) : ℂ))) := by
  rw [paperFT_shiftComb hf, paperFT_shiftComb hf]
  have hp : I * gammaOf ρ = pe ρ := I_mul_gammaOf _
  have hc : conj (cexp (I * conj (gammaOf ρ) * τ)) = cexp (pe ρ * ((-τ : ℝ) : ℂ)) := by
    rw [← Complex.exp_conj]; congr 1
    rw [← hp]
    simp only [map_mul, Complex.conj_I, Complex.conj_conj, Complex.conj_ofReal, ofReal_neg]
    ring
  have hE : cexp (pe ρ * τ) * cexp (pe ρ * ((-τ : ℝ) : ℂ)) = 1 := by
    rw [← Complex.exp_add, ofReal_neg, mul_neg, add_neg_cancel, Complex.exp_zero]
  rw [map_mul, map_add, map_one, map_mul, hc, show I * gammaOf ρ * τ = pe ρ * τ by rw [hp]]
  unfold wt
  linear_combination ((zetaZeroConfig.mult ρ : ℂ) * paperFT f (gammaOf ρ) *
    conj (paperFT f (conj (gammaOf ρ))) * c * conj c) * hE

/-! ### The exponential sum and its symmetry -/

/-- `D τ = Σ_ρ w_ρ e^{(ρ − 1/2) τ}`. -/
abbrev D (f : ℝ → ℂ) (τ : ℝ) : ℂ := Core.expSum pe (wt f) τ

lemma summable_wt {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) : Summable (wt f) :=
  (weilQFix_eq_zeroSum hf).1

lemma summable_norm_wt {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) : Summable (fun ρ => ‖wt f ρ‖) :=
  summable_norm_iff.mpr (summable_wt hf)

lemma summable_wt_exp {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (τ : ℝ) :
    Summable (fun ρ => wt f ρ * cexp (pe ρ * τ)) :=
  Summable.of_norm_bounded ((summable_norm_wt hf).mul_right (Real.exp |τ|))
    (fun ρ => Core.term_norm_le abs_pe_re_le ρ τ)

lemma D_neg (f : ℝ → ℂ) (τ : ℝ) : D f (-τ) = conj (D f τ) := by
  unfold D Core.expSum
  rw [conj_tsum, ← Equiv.tsum_eq reflEquiv (fun ρ => wt f ρ * cexp (pe ρ * ((-τ : ℝ) : ℂ)))]
  congr 1; funext ρ
  rw [wt_refl, pe_refl, map_mul, ← Complex.exp_conj]
  congr 2
  simp only [map_mul, Complex.conj_ofReal, ofReal_neg]
  ring

lemma D_zero_conj (f : ℝ → ℂ) : conj (D f 0) = D f 0 := by
  rw [← D_neg, neg_zero]

/-- Zero-sum expansion of `weilQFix (f + c·f(· − τ))`. -/
lemma weilQFix_shiftComb {a : ℝ} {f : ℝ → ℂ} (hf : TestFn a f) (c : ℂ) (τ : ℝ) :
    weilQFix (shiftComb f c τ) = (1 + c * conj c) * D f 0 + c * D f τ + conj c * conj (D f τ) := by
  rw [(weilQFix_eq_zeroSum (testFn_shiftComb hf c τ)).2, tsum_congr (term_eq hf c τ), ← D_neg]
  have s1 := (summable_wt hf).mul_left (1 + c * conj c)
  have s2 := (summable_wt_exp hf τ).mul_left c
  have s3 := (summable_wt_exp hf (-τ)).mul_left (conj c)
  rw [Summable.tsum_add (s1.add s2) s3, Summable.tsum_add s1 s2, tsum_mul_left, tsum_mul_left,
    tsum_mul_left]
  unfold D Core.expSum
  simp only [ofReal_zero, mul_zero, Complex.exp_zero, mul_one]

/-- **Positivity bounds the exponential sum**: `‖D τ‖ ≤ Re D 0`. -/
lemma norm_D_le {a : ℝ} (ha : 0 < a) {f : ℝ → ℂ} (hf : TestFn a f)
    (hW : ∀ L : ℝ, 0 < L → WeilPosFix L) (τ : ℝ) : ‖D f τ‖ ≤ (D f 0).re := by
  have hpos : ∀ c : ℂ, 0 ≤ (weilQFix (shiftComb f c τ)).re := by
    intro c
    refine hW (2 * (a + |τ|)) (by positivity) _ ?_
    rw [show 2 * (a + |τ|) / 2 = a + |τ| by ring]
    exact testFn_shiftComb hf c τ
  by_cases hD : D f τ = 0
  · have := hpos 0
    rw [weilQFix_shiftComb hf] at this
    simp only [map_zero, mul_zero, add_zero, zero_mul, one_mul] at this
    rw [hD, norm_zero]; exact this
  · set Dt := D f τ with hDt
    set n : ℝ := ‖Dt‖ with hn
    have hn0 : 0 < n := norm_pos_iff.mpr hD
    have hnc : (n : ℂ) ≠ 0 := ofReal_ne_zero.mpr hn0.ne'
    have hnn : conj Dt * Dt = (n : ℂ) ^ 2 := by
      rw [mul_comm, Complex.mul_conj, Complex.normSq_eq_norm_sq]; push_cast; rfl
    set c : ℂ := -conj Dt / (n : ℂ) with hc
    have hcc : conj c = -Dt / (n : ℂ) := by
      rw [hc, map_div₀, map_neg, Complex.conj_conj, Complex.conj_ofReal]
    have h1 : c * conj c = 1 := by
      rw [hcc, hc]; field_simp; try rw [hnn]
    have h2 : c * Dt = -(n : ℂ) := by
      rw [hc]; field_simp; try rw [hnn]; try ring
    have h3 : conj c * conj Dt = -(n : ℂ) := by
      rw [← map_mul, h2, map_neg, Complex.conj_ofReal]
    have := hpos c
    rw [weilQFix_shiftComb hf, h1, h2, h3] at this
    have hre : ((1 + 1) * D f 0 + -(n : ℂ) + -(n : ℂ)).re = 2 * (D f 0).re - 2 * n := by
      simp; ring
    rw [hre] at this
    linarith

/-! ### The test function: bump × `e^{−i t₀ u}` -/

/-- A fixed smooth bump on `ℝ`, `1` on `[−1/2, 1/2]`, supported in `(−1, 1)`. -/
def bump : ContDiffBump (0 : ℝ) := ⟨1 / 2, 1, by norm_num, by norm_num⟩

/-- `f₁(u) = bump(u) e^{−i t₀ u}`. -/
def f1 (t0 : ℝ) : ℝ → ℂ := fun u => (bump u : ℂ) * cexp (-(I * t0 * u))

lemma testFn_f1 (t0 : ℝ) : TestFn 1 (f1 t0) := by
  have h1 : ContDiff ℝ ((⊤ : ℕ∞) : WithTop ℕ∞) (fun u : ℝ => (bump u : ℂ)) :=
    ofRealCLM.contDiff.comp bump.contDiff
  have h2 : ContDiff ℝ ((⊤ : ℕ∞) : WithTop ℕ∞) (fun u : ℝ => -(I * (t0 : ℂ) * (u : ℂ))) :=
    (contDiff_const.mul ofRealCLM.contDiff).neg
  refine ⟨h1.mul h2.cexp, ?_⟩
  intro u hu
  have hb : bump u ≠ 0 := by
    intro h0; apply hu; simp [f1, h0]
  have : u ∈ Metric.ball (0 : ℝ) 1 := by
    have := bump.support_eq ▸ Function.mem_support.mpr hb
    simpa [bump] using this
  rw [Metric.mem_ball, dist_zero_right, Real.norm_eq_abs] at this
  exact ⟨by linarith [neg_abs_le u], by linarith [le_abs_self u]⟩

lemma paperFT_f1 (t0 b : ℝ) :
    paperFT (f1 t0) ((t0 : ℂ) - I * b) = ((∫ u : ℝ, bump u * Real.exp (b * u) : ℝ) : ℂ) := by
  rw [← integral_complex_ofReal]
  unfold paperFT f1
  congr 1; funext u
  push_cast
  rw [mul_assoc, ← Complex.exp_add]
  congr 2
  ring_nf
  rw [I_sq]
  ring

lemma integral_bump_exp_pos (b : ℝ) : 0 < ∫ u : ℝ, bump u * Real.exp (b * u) := by
  refine Continuous.integral_pos_of_hasCompactSupport_nonneg_nonzero (x := 0)
    (bump.continuous.mul (by fun_prop)) bump.hasCompactSupport.mul_right
    (fun u => mul_nonneg bump.nonneg (Real.exp_pos _).le) ?_
  rw [bump.one_of_mem_closedBall (by simp [bump])]
  simp

/-! ### Choosing a right-isolated off-line zero -/

lemma exists_sep (s : Set ℂ) (hs : s.Finite) (t0 : ℝ) (h : ∀ z ∈ s, z.im ≠ t0) :
    ∃ η : ℝ, 0 < η ∧ η ≤ 1 ∧ ∀ z ∈ s, η ≤ |z.im - t0| := by
  rcases s.eq_empty_or_nonempty with he | hne
  · exact ⟨1, one_pos, le_rfl, by simp [he]⟩
  · obtain ⟨z0, hz0, hmin⟩ := Set.exists_min_image s (fun z => |z.im - t0|) hs hne
    have hpos : 0 < |z0.im - t0| := abs_pos.mpr (sub_ne_zero.mpr (h z0 hz0))
    exact ⟨min 1 |z0.im - t0|, lt_min one_pos hpos, min_le_left _ _,
      fun z hz => (min_le_right _ _).trans (hmin z hz)⟩

/-- **Weil's criterion, sufficiency** (intended functional `weilQFix`): positivity on every window
implies the Riemann hypothesis. -/
theorem RH_of_weilPosFix (h : ∀ L : ℝ, 0 < L → WeilPosFix L) : RiemannHypothesis := by
  rw [Fayed.R7.RH_iff_QuasiRH_half]
  intro s hs hs1 hz
  have hs' : s.re < 1 := by
    by_contra hc
    exact riemannZeta_ne_zero_of_one_le_re (not_lt.mp hc) hz
  have hρ1 : s ∈ zetaZeroConfig.carrier := ⟨hz, by linarith, hs'⟩
  -- ρ₀: largest real part in the ordinate window around s, among zeros with Re ≥ Re s
  set W := zetaZeroConfig.carrier ∩ {ρ : ℂ | s.im - 1 < ρ.im ∧ ρ.im ≤ s.im + 1} ∩
    {ρ : ℂ | s.re ≤ ρ.re} with hW
  have hWfin : W.Finite := (zetaZeroConfig.finite_window _ _).subset inter_subset_left
  obtain ⟨ρ0, hρ0W, hmax⟩ := Set.exists_max_image W (fun ρ => ρ.re) hWfin
    ⟨s, ⟨hρ1, by linarith, by linarith⟩, (le_rfl : s.re ≤ s.re)⟩
  have hρ0 : ρ0 ∈ zetaZeroConfig.carrier := hρ0W.1.1
  set t0 := ρ0.im with ht0
  have hb : 0 < ρ0.re - 1 / 2 := by have := hρ0W.2; simp only [mem_ofPred_eq] at this; linarith
  -- separation η
  set B := zetaZeroConfig.carrier ∩ {ρ : ℂ | t0 - 1 < ρ.im ∧ ρ.im ≤ t0 + 1} ∩
    {ρ : ℂ | ρ0.re < ρ.re} with hB
  have hBfin : B.Finite := (zetaZeroConfig.finite_window _ _).subset inter_subset_left
  have hBne : ∀ z ∈ B, z.im ≠ t0 := by
    intro z hz hzt
    have hzW : z ∈ W := by
      refine ⟨⟨hz.1.1, ?_, ?_⟩, ?_⟩
      · rw [hzt]; exact hρ0W.1.2.1
      · rw [hzt]; exact hρ0W.1.2.2
      · have h1 : ρ0.re < z.re := hz.2
        have h2 : s.re ≤ ρ0.re := hρ0W.2
        show s.re ≤ z.re
        linarith
    have := hmax z hzW
    have := hz.2
    simp only [mem_ofPred_eq] at *
    linarith
  obtain ⟨η, hη0, hη1, hηB⟩ := exists_sep B hBfin t0 hBne
  have hsep : ∀ ρ ∈ zetaZeroConfig.carrier, ρ0.re < ρ.re → η ≤ |ρ.im - t0| := by
    intro ρ hρ hre
    by_cases hc : |ρ.im - t0| < 1
    · obtain ⟨h1, h2⟩ := abs_lt.mp hc
      exact hηB ρ ⟨⟨hρ, by linarith, by linarith⟩, hre⟩
    · push Not at hc; linarith
  -- the core lemma
  have hf := testFn_f1 t0
  have hw0 := Core.weight_eq_zero (summable_norm_wt hf) abs_pe_re_le
    (M := (D (f1 t0) 0).re) (fun τ _ => norm_D_le one_pos hf h τ)
    (⟨ρ0, hρ0⟩ : zetaZeroConfig.carrier)
    (by rw [pe_re]; exact hb) hη0
    (fun ρ hρ => by
      rw [pe_re, pe_re] at hρ
      rw [pe_im, pe_im]
      exact hsep ρ ρ.2 (by simpa using hρ))
    (fun ρ hρ => Subtype.ext (by simpa [pe] using hρ))
  -- but the weight at ρ₀ is nonzero
  have hγ : gammaOf ρ0 = (t0 : ℂ) - I * ((ρ0.re - 1 / 2 : ℝ) : ℂ) := by
    apply Complex.ext
    · rw [WeilEF.gammaOf_re]; simp [ht0]
    · rw [WeilEF.gammaOf_im]; simp
  have hγc : conj (gammaOf ρ0) = (t0 : ℂ) - I * ((-(ρ0.re - 1 / 2) : ℝ) : ℂ) := by
    rw [hγ]; apply Complex.ext <;> simp
  have hm : (1 : ℕ) ≤ zetaZeroConfig.mult ρ0 := zetaZeroConfig.one_le_mult ρ0 hρ0
  have hw0' : wt (f1 t0) ⟨ρ0, hρ0⟩ ≠ 0 := by
    unfold wt
    simp only
    rw [hγc, hγ, paperFT_f1, paperFT_f1, Complex.conj_ofReal]
    refine mul_ne_zero (Nat.cast_ne_zero.mpr (by omega)) (mul_ne_zero ?_ ?_)
    · exact ofReal_ne_zero.mpr (integral_bump_exp_pos _).ne'
    · exact ofReal_ne_zero.mpr (integral_bump_exp_pos _).ne'
  exact hw0' hw0

/-- **Weil's criterion, necessity**: under RH every zero has real `γ_ρ`, so each term of the zero
sum is `m_ρ |F(γ_ρ)|² ≥ 0`. -/
theorem weilPosFix_of_RH (hRH : RiemannHypothesis) (L : ℝ) : WeilPosFix L := by
  intro f hf
  obtain ⟨hs, he⟩ := weilQFix_eq_zeroSum hf
  rw [he, Complex.re_tsum hs]
  refine tsum_nonneg fun ρ => ?_
  have hre : (ρ : ℂ).re = 1 / 2 := RH_implies_on_line hRH ρ.2
  have hγ : conj (gammaOf ρ) = gammaOf ρ :=
    Complex.ext (by simp) (by rw [conj_im, WeilEF.gammaOf_im, hre]; norm_num)
  rw [hγ, Complex.mul_conj]
  simp only [mul_re, natCast_re, ofReal_re, natCast_im, ofReal_im, mul_zero, sub_zero]
  exact mul_nonneg (Nat.cast_nonneg _) (Complex.normSq_nonneg _)

/-- **Weil's criterion** (both directions, sorry-free, standard axioms). -/
theorem weil_criterion : RiemannHypothesis ↔ ∀ L : ℝ, 0 < L → WeilPosFix L :=
  ⟨fun h L _ => weilPosFix_of_RH h L, RH_of_weilPosFix⟩

/-- The brief's statement, for `Riemann/Deep/Weil/Defs.lean`'s (fixed) `WeilPos`. -/
theorem RH_of_weilPos (h : ∀ L : ℝ, 0 < L → WeilPos L) : RiemannHypothesis :=
  RH_of_weilPosFix fun L hL => (weilPos_iff_fix L).mp (h L hL)

/-- **Weil's criterion** for `Riemann/Deep/Weil/Defs.lean`'s `WeilPos`. -/
theorem weil_criterion_weilPos : RiemannHypothesis ↔ ∀ L : ℝ, 0 < L → WeilPos L :=
  weil_criterion

/-- `Defs.WeilCriterion` (all windows `L : ℝ`) holds outright, sorry-free. -/
theorem weilCriterion_holds : WeilCriterion :=
  ⟨fun h L => weilPosFix_of_RH h L, fun h => RH_of_weilPos fun L _ => h L⟩

end Fayed.WeilSuff
