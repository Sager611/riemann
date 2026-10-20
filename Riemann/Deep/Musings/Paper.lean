import Mathlib

/-!
# Nadim, arXiv 2605.30767 ("Musings on the Riemann Hypothesis"), audited

The paper claims NO proof of RH. Its chain (`z = s - 1/2`, `ξ = f + i g`):
* §2  `f` even, `g` odd in each variable; Taylor monomials of `z^{2n}`.  (one typo: `Im z⁴`)
* §3  Lemma: zero contours `f = 0` / `g = 0` emanating from the critical line never meet
      off it ("mixed Dirichlet/Neumann region ⇒ f ≡ 0").                 CORRECT — `no_enclosed_region`
* §4  Strip model `f = Σ bₙ sin(ny) cosh(nx)`; examples with off-axis zeros.  (fine for finite sums)
* §5  Taylor criterion: `(-1)^k F^{(2k)}(y) ≥ 0 ∀k ⇒ f(·,y) ≠ 0`.         CORRECT — `taylor_criterion`
      but VACUOUS for `ξ` wherever `Re ξ(x+iy)` has a zero — `criterion_fails`; and the "different
      harmonic function" `Σ sin((2k+1)y)cosh((2k+1)x)/(2k+1)³` diverges for every `x > 0` —
      `sec5_series_not_summable`.
-/

namespace Riemann.Deep.Musings

open Complex Filter Topology

/-! ## §2: `Im z⁴` -/

theorem im_z4 (z : ℂ) : (z ^ 4).im = 4 * z.re ^ 3 * z.im - 4 * z.re * z.im ^ 3 := by
  simp only [pow_succ, pow_zero, one_mul, mul_im, mul_re]; ring

/-- The paper's displayed `Im z⁴ = x³y − xy³` is false (typo; factor 4 missing). -/
theorem paper_z4_false : ¬ ∀ z : ℂ, (z ^ 4).im = z.re ^ 3 * z.im - z.re * z.im ^ 3 := by
  intro h
  have h1 := im_z4 (2 + I)
  have h2 := h (2 + I)
  simp at h1 h2
  rw [h1] at h2; norm_num at h2

/-! ## §3: the non-intersection lemma is correct

Sharper form: if a bounded open `U` has `∂U ⊆ {Re f · Im f = 0}` for entire `f`, then
`Im f² = 2 Re f Im f` vanishes on `∂U`, hence (max principle for `exp(± i f²)`) on `U`,
hence `f²` is real on an open set, hence constant. Two axis-emanating contours of `ξ` meeting
off the axis would, with a critical-line segment (`g = 0` there), bound such a `U`. -/

theorem im_sq (w : ℂ) : (w ^ 2).im = 2 * w.re * w.im := by
  simp only [sq, mul_im]; ring

theorem im_sq_zero_on_closure {f : ℂ → ℂ} (hf : Differentiable ℂ f) {U : Set ℂ}
    (hU : Bornology.IsBounded U) (hfr : ∀ z ∈ frontier U, (f z).re * (f z).im = 0) :
    ∀ z ∈ closure U, (f z ^ 2).im = 0 := by
  have key : ∀ c : ℂ, (c = I ∨ c = -I) → ∀ z ∈ closure U, ‖cexp (c * f z ^ 2)‖ ≤ 1 := by
    intro c hc
    have hd : Differentiable ℂ (fun z => cexp (c * f z ^ 2)) :=
      ((hf.pow 2).const_mul c).cexp
    intro z hz
    refine Complex.norm_le_of_forall_mem_frontier_norm_le hU hd.diffContOnCl ?_ hz
    intro z hz
    rw [Complex.norm_exp]
    have h0 : (f z ^ 2).im = 0 := by rw [im_sq]; linarith [hfr z hz]
    rcases hc with rfl | rfl <;> simp [mul_re, h0]
  intro z hz
  have h1 := key I (Or.inl rfl) z hz
  have h2 := key (-I) (Or.inr rfl) z hz
  rw [Complex.norm_exp, Real.exp_le_one_iff] at h1 h2
  simp [mul_re] at h1 h2
  linarith

/-- The real line is not a neighbourhood of any point. -/
theorem setOf_im_eq_zero_not_mem_nhds (w : ℂ) : {u : ℂ | u.im = 0} ∉ 𝓝 w := by
  intro h
  obtain ⟨ε, hε, hb⟩ := Metric.mem_nhds_iff.mp h
  have hw : w.im = 0 := hb (Metric.mem_ball_self hε)
  have hp : w + ((ε / 2 : ℝ) : ℂ) * I ∈ Metric.ball w ε := by
    rw [Metric.mem_ball, dist_eq_norm, add_sub_cancel_left]
    simp [abs_of_pos hε]
    linarith
  have h1 : (w + ((ε / 2 : ℝ) : ℂ) * I).im = 0 := hb hp
  simp only [add_im, mul_im, ofReal_re, ofReal_im, I_re, I_im] at h1
  linarith

/-- §3 lemma (correct): an entire `f` admits no nonempty bounded open region whose boundary
lies on `{Re f = 0} ∪ {Im f = 0}` unless `f²` is constant. -/
theorem no_enclosed_region {f : ℂ → ℂ} (hf : Differentiable ℂ f) {U : Set ℂ} (hUo : IsOpen U)
    (hU : Bornology.IsBounded U) {z₀ : ℂ} (hz₀ : z₀ ∈ U)
    (hfr : ∀ z ∈ frontier U, (f z).re * (f z).im = 0) : ∀ z, f z ^ 2 = f z₀ ^ 2 := by
  have him := im_sq_zero_on_closure hf hU hfr
  have hg : Differentiable ℂ (fun z => f z ^ 2) := hf.pow 2
  rcases (hg.analyticAt z₀).eventually_constant_or_nhds_le_map_nhds with hc | hn
  · have := (show AnalyticOnNhd ℂ (fun z => f z ^ 2) Set.univ from
        fun z _ => hg.analyticAt z).eqOn_of_preconnected_of_eventuallyEq
      (analyticOnNhd_const (v := f z₀ ^ 2)) isPreconnected_univ (Set.mem_univ z₀) hc
    exact fun z => this (Set.mem_univ z)
  · exfalso
    have hmem : {u : ℂ | u.im = 0} ∈ map (fun z => f z ^ 2) (𝓝 z₀) :=
      Filter.mem_map.mpr (Filter.mem_of_superset (hUo.mem_nhds hz₀)
        fun z hz => show (f z ^ 2).im = 0 from him z (subset_closure hz))
    exact setOf_im_eq_zero_not_mem_nhds _ (hn hmem)

/-- The paper's §3 chain as Props. `Encloses`: two axis-emanating contours meeting off the
axis bound a nonempty bounded open region with boundary on `{f g = 0}` (Jordan topology). -/
def Encloses (ξ : ℂ → ℂ) (P : Prop) : Prop :=
  P → ∃ U : Set ℂ, IsOpen U ∧ Bornology.IsBounded U ∧ U.Nonempty ∧
    ∀ z ∈ frontier U, (ξ z).re * (ξ z).im = 0

theorem sec3_chain {ξ : ℂ → ℂ} (hξ : Differentiable ℂ ξ) (hnc : ∃ a b, ξ a ^ 2 ≠ ξ b ^ 2)
    {Meet : Prop} (henc : Encloses ξ Meet) : ¬ Meet := by
  intro hm
  obtain ⟨U, hUo, hUb, ⟨z₀, hz₀⟩, hfr⟩ := henc hm
  obtain ⟨a, b, hab⟩ := hnc
  have h := no_enclosed_region hξ hUo hUb hz₀ hfr
  exact hab ((h a).trans (h b).symm)

/-! ## §5: Taylor criterion — correct, but vacuous for `ξ` -/

/-- `f(x,y) = Σ_k x^{2k}/(2k)! · c_k` with `c_k = (-1)^k F^{(2k)}(y)` (valid for entire `ξ`).
If every `c_k ≥ 0` then `f(x,y) ≥ F(y)`. -/
theorem taylor_criterion (c : ℕ → ℝ) (x : ℝ) (hc : ∀ k, 0 ≤ c k)
    (hs : Summable fun k => x ^ (2 * k) / (2 * k).factorial * c k) :
    c 0 ≤ ∑' k, x ^ (2 * k) / (2 * k).factorial * c k := by
  have h := hs.le_tsum 0 (fun k _ => by
    have : 0 ≤ x ^ (2 * k) := by rw [pow_mul]; positivity
    exact mul_nonneg (div_nonneg this (by positivity)) (hc k))
  simpa using h

/-- Hence wherever `Re ξ(x+iy) = 0` for some `x` while `F(y) > 0`, the criterion fails at `y`.
For `ξ`, `arg Γ(s/2)` grows like `(y/2) log x` as `x → ∞`, so `Re ξ(·+iy)` changes sign infinitely
often on EVERY horizontal line `y ≠ 0` (numerics `.scratch/musings/signs.py`): the criterion holds
at no `y` with `F(y) ≠ 0`, so "if `F` had this property everywhere, RH would hold" is vacuous. -/
theorem criterion_fails (c : ℕ → ℝ) (x : ℝ) (hF : 0 < c 0)
    (hs : Summable fun k => x ^ (2 * k) / (2 * k).factorial * c k)
    (hzero : ∑' k, x ^ (2 * k) / (2 * k).factorial * c k = 0) : ¬ ∀ k, 0 ≤ c k := by
  intro hc
  have := taylor_criterion c x hc hs
  linarith

/-- §5: the paper's "different harmonic function" `Σ sin((2k+1)y) cosh((2k+1)x)/(2k+1)³`
diverges at `y = π/2` for every `x > 0` (terms are unbounded), so it defines no function: the
Cauchy problem `f = F, ∂ₓf = 0` on `x = 0` with `F = (π/8)y(π−y)` has no solution. -/
theorem sec5_series_not_summable {x : ℝ} (hx : 0 < x) :
    ¬ Summable (fun k : ℕ => Real.sin ((2 * k + 1) * (Real.pi / 2)) *
      Real.cosh ((2 * k + 1) * x) / (2 * k + 1) ^ 3) := by
  intro h
  have ht := h.tendsto_atTop_zero.norm
  simp only [norm_zero] at ht
  have hsin : ∀ k : ℕ, |Real.sin ((2 * k + 1) * (Real.pi / 2))| = 1 := by
    intro k
    have : (2 * (k : ℝ) + 1) * (Real.pi / 2) = k * Real.pi + Real.pi / 2 := by ring
    rw [this, Real.sin_add_pi_div_two, Real.cos_nat_mul_pi]; simp
  -- lower bound by `x³/2 · exp t / t³` with `t = (2k+1)x → ∞`
  have hlin : Tendsto (fun k : ℕ => (2 * (k : ℝ) + 1) * x) atTop atTop := by
    apply Tendsto.atTop_mul_const hx
    exact tendsto_atTop_add_const_right _ 1
      (Tendsto.const_mul_atTop (by norm_num) tendsto_natCast_atTop_atTop)
  have hbig : Tendsto (fun k : ℕ => x ^ 3 / 2 * (Real.exp ((2 * (k : ℝ) + 1) * x) /
      ((2 * (k : ℝ) + 1) * x) ^ 3)) atTop atTop :=
    ((Real.tendsto_exp_div_pow_atTop 3).comp hlin).const_mul_atTop (by positivity)
  have hge : ∀ k : ℕ, x ^ 3 / 2 * (Real.exp ((2 * (k : ℝ) + 1) * x) /
      ((2 * (k : ℝ) + 1) * x) ^ 3) ≤ ‖Real.sin ((2 * k + 1) * (Real.pi / 2)) *
      Real.cosh ((2 * k + 1) * x) / (2 * k + 1) ^ 3‖ := by
    intro k
    have hn : (0 : ℝ) < 2 * k + 1 := by positivity
    have hx3 : 0 < x ^ 3 := by positivity
    have heq : x ^ 3 / 2 * (Real.exp ((2 * k + 1) * x) / ((2 * k + 1) * x) ^ 3) =
        Real.exp ((2 * k + 1) * x) / 2 / (2 * k + 1) ^ 3 := by
      rw [mul_pow]; field_simp
    rw [heq, Real.norm_eq_abs, abs_div, abs_mul, hsin, one_mul, abs_of_pos (Real.cosh_pos _),
      abs_of_pos (by positivity)]
    have hc : Real.exp ((2 * k + 1) * x) / 2 ≤ Real.cosh ((2 * k + 1) * x) := by
      rw [Real.cosh_eq]; have := Real.exp_pos (-((2 * k + 1) * x)); linarith
    gcongr
  exact not_tendsto_nhds_of_tendsto_atTop (tendsto_atTop_mono hge hbig) 0 ht

end Riemann.Deep.Musings
