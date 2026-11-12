import Riemann.Deep.NBLp.Lean.Basic
import Mathlib.MeasureTheory.Function.LpSeminorm.CompareExp
import Mathlib.NumberTheory.LSeries.Nonvanishing

/-!
# NBLp in Lean, part 2: PROOF.md §0 definitions and Proposition N

* rho k x = {1/(kx)} (Báez-Duarte's ρ_k on (0,∞)), fk k x = {1/(kx)} − k⁻¹{1/x} (f_k on (0,1)),
  fk_eq_rhoNB: f_k = ρ_{1/k} in the page's notation ρ_a(x) = {a/x} − a{1/x}.
* ChiInW p: χ = 1_(0,1) lies in W_p, the L^p(0,1)-closure of V = span_ℂ{f_k : k ≥ 2}
  (every ε > 0 admits N and c_2..c_N ∈ ℂ with ‖χ − Σ_{k=2}^N c_k f_k‖_p < ε, i.e. d_{N,p} → 0).
* propN (PROOF.md §3, Proposition N): ζ(w₀) = 0, Re w₀ > 1/p, 1 < p ⇒ ¬ ChiInW p
  (Hölder against x^{w₀−1} ∈ L^q(0,1)).
-/

open Complex MeasureTheory Set Filter Topology
open scoped ENNReal

namespace NBLp

/-- PROOF.md §0: ρ_k(x) = {1/(kx)}, x > 0. -/
noncomputable def rho (k : ℕ) (x : ℝ) : ℝ := Int.fract (1 / (k * x))

/-- PROOF.md §0: f_k(x) = {1/(kx)} − k⁻¹{1/x}, 0 < x < 1. -/
noncomputable def fk (k : ℕ) (x : ℝ) : ℝ := Int.fract (1 / (k * x)) - (k : ℝ)⁻¹ * Int.fract (1 / x)

/-- Lebesgue measure on (0,1] (L^p(0,1); the endpoint is null). -/
noncomputable abbrev vol01 : Measure ℝ := volume.restrict (Ioc (0:ℝ) 1)

/-- χ − Σ_{k=2}^N c_k f_k on (0,1). -/
noncomputable def approx (N : ℕ) (c : ℕ → ℂ) (x : ℝ) : ℂ :=
  1 - ∑ k ∈ Finset.Icc 2 N, c k * (fk k x : ℂ)

/-- PROOF.md §1 Theorem A (ii): χ ∈ W_p, i.e. d_{N,p} → 0. -/
def ChiInW (p : ℝ) : Prop :=
  ∀ ε : ℝ, 0 < ε → ∃ (N : ℕ) (c : ℕ → ℂ),
    eLpNorm (approx N c) (ENNReal.ofReal p) vol01 < ENNReal.ofReal ε

/-- PROOF.md §1 Theorem A (i) with θ = 1/p: ζ(s) ≠ 0 for Re s > θ. -/
def ZeroFree (θ : ℝ) : Prop := ∀ s : ℂ, θ < s.re → riemannZeta s ≠ 0

/-- The page's ρ_a(x) = {a/x} − a{1/x} at a = 1/k is PROOF.md's f_k. -/
lemma fk_eq_rhoNB (k : ℕ) (x : ℝ) : fk k x = rhoNB (1 / k) x := by
  unfold fk rhoNB
  rw [div_div]
  simp only [one_div]

lemma measurable_fk (k : ℕ) : Measurable (fk k) := by
  have : fk k = rhoNB (1 / k) := funext (fk_eq_rhoNB k)
  rw [this]; exact measurable_rhoNB _

lemma measurable_approx (N : ℕ) (c : ℕ → ℂ) : Measurable (approx N c) :=
  measurable_const.sub (Finset.measurable_sum _ fun k _ =>
    measurable_const.mul (Complex.measurable_ofReal.comp (measurable_fk k)))

/-- PROOF.md §3, Proposition N (necessity, (ii) ⇒ (i)). -/
theorem propN {p : ℝ} (hp : 1 < p) {w : ℂ} (hz : riemannZeta w = 0) (hw : 1 / p < w.re) :
    ¬ ChiInW p := by
  intro hW
  have hp0 : 0 < p := by linarith
  have hp1 : p - 1 ≠ 0 := by intro h; linarith
  have hw0 : 0 < w.re := lt_trans (by positivity) hw
  have hw1 : w.re < 1 := by
    by_contra h; exact riemannZeta_ne_zero_of_one_le_re (not_lt.mp h) hz
  have hwne : w ≠ 0 := by rintro rfl; simp at hw0
  set q : ℝ := p / (p - 1) with hqdef
  have hq0 : 0 < q := div_pos hp0 (by linarith)
  have hpq : p⁻¹ + q⁻¹ = 1 := by rw [hqdef, inv_div]; field_simp; ring
  have : ENNReal.HolderTriple (ENNReal.ofReal p) (ENNReal.ofReal q) 1 := ⟨by
    rw [← ENNReal.ofReal_inv_of_pos hp0, ← ENNReal.ofReal_inv_of_pos hq0,
      ← ENNReal.ofReal_add (by positivity) (by positivity), hpq, ENNReal.ofReal_one, inv_one]⟩
  set g : ℝ → ℂ := fun x => (x : ℂ) ^ (w - 1) with hg
  have hgint : IntegrableOn g (Ioc 0 1) := by
    have := intervalIntegral.intervalIntegrable_cpow' (a := 0) (b := 1) (r := w - 1)
      (by simp; linarith)
    rwa [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
  have hgI : ∫ x in Ioc (0:ℝ) 1, g x = 1 / w := by
    rw [← intervalIntegral.integral_of_le zero_le_one,
      integral_cpow (Or.inl (by simp; linarith)), sub_add_cancel, ofReal_one, ofReal_zero,
      one_cpow, zero_cpow hwne, sub_zero]
  -- x^{w-1} ∈ L^q(0,1): q(Re w - 1) > -1 ⟺ Re w > 1/p
  have hexp : -1 < (w.re - 1) * q := by
    have h1 : 1 < w.re * p := by rwa [div_lt_iff₀ hp0] at hw
    have e : (w.re - 1) * q + 1 = (w.re * p - 1) / (p - 1) := by
      rw [hqdef]; field_simp; ring
    have : 0 < (w.re * p - 1) / (p - 1) := div_pos (by linarith) (by linarith)
    linarith
  have hGfin : eLpNorm g (ENNReal.ofReal q) vol01 < ⊤ := by
    rw [eLpNorm_eq_lintegral_rpow_enorm_toReal (by simp [hq0]) ENNReal.ofReal_ne_top
      hgint.aestronglyMeasurable, ENNReal.toReal_ofReal hq0.le]
    refine ENNReal.rpow_lt_top_of_nonneg (by positivity) ?_
    have hint : IntegrableOn (fun x : ℝ => x ^ ((w.re - 1) * q)) (Ioc 0 1) := by
      have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := 1)
        (r := (w.re - 1) * q) hexp
      rwa [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
    have heq : ∀ x ∈ Ioc (0:ℝ) 1, ‖g x‖ₑ ^ q = ENNReal.ofReal (x ^ ((w.re - 1) * q)) := by
      intro x hx
      rw [← ofReal_norm, hg, norm_cpow_eq_rpow_re_of_pos hx.1,
        ENNReal.ofReal_rpow_of_nonneg (Real.rpow_nonneg hx.1.le _) hq0.le,
        ← Real.rpow_mul hx.1.le]
      simp
    rw [setLIntegral_congr_fun measurableSet_Ioc heq]
    exact (hint.lintegral_lt_top).ne
  set G := eLpNorm g (ENNReal.ofReal q) vol01 with hGdef
  have hGne : G ≠ ⊤ := hGfin.ne
  set c0 : ℝ := ‖1 / w‖ with hc0
  have hc0pos : 0 < c0 := norm_pos_iff.mpr (one_div_ne_zero hwne)
  obtain ⟨N, c, hNc⟩ := hW (c0 / (G.toReal + 1)) (by positivity)
  set Fn := approx N c with hFn
  have hrho : ∀ k ∈ Finset.Icc 2 N, Integrable (fun x : ℝ => (fk k x : ℂ) * g x) vol01 := by
    intro k hk
    have hk2 : 2 ≤ k := (Finset.mem_Icc.mp hk).1
    have hkpos : (0:ℝ) < k := by exact_mod_cast (by omega : 0 < k)
    have hθ : (0:ℝ) < 1 / k := by positivity
    have := integrableOn_rho hθ hw0 hw1
    simpa [fk_eq_rhoNB, hg, IntegrableOn] using this
  have hzero : ∀ k ∈ Finset.Icc 2 N, ∫ x in Ioc (0:ℝ) 1, (fk k x : ℂ) * g x = 0 := by
    intro k hk
    have hk2 : 2 ≤ k := (Finset.mem_Icc.mp hk).1
    have hkpos : (0:ℝ) < k := by exact_mod_cast (by omega : 0 < k)
    have hθ : (0:ℝ) < 1 / k := by positivity
    have hθ1 : (1:ℝ) / k ≤ 1 := by
      rw [div_le_one hkpos]; exact_mod_cast (by omega : 1 ≤ k)
    have := integral_rho_eq_zero hθ hθ1 hw0 hw1 hz
    simpa [fk_eq_rhoNB, hg] using this
  have hFI : ∫ x in Ioc (0:ℝ) 1, Fn x * g x = 1 / w := by
    have e : ∀ x, Fn x * g x = g x - ∑ k ∈ Finset.Icc 2 N, c k * ((fk k x : ℂ) * g x) := by
      intro x
      simp only [hFn, approx]
      rw [sub_mul, one_mul, Finset.sum_mul]
      congr 1
      exact Finset.sum_congr rfl fun k _ => by ring
    simp_rw [e]
    rw [integral_sub hgint (integrable_finsetSum _ fun k hk => (hrho k hk).const_mul _),
      integral_finsetSum _ fun k hk => (hrho k hk).const_mul _,
      Finset.sum_eq_zero fun k hk => by rw [integral_const_mul, hzero k hk, mul_zero],
      sub_zero, hgI]
  have hH : ‖∫ x in Ioc (0:ℝ) 1, Fn x * g x‖ₑ ≤ eLpNorm Fn (ENNReal.ofReal p) vol01 * G := by
    calc ‖∫ x in Ioc (0:ℝ) 1, Fn x * g x‖ₑ ≤ ∫⁻ x in Ioc (0:ℝ) 1, ‖Fn x * g x‖ₑ :=
          enorm_integral_le_lintegral_enorm _
      _ ≤ eLpNorm (Fn • g) 1 vol01 := lintegral_enorm_le_eLpNorm_one (f := Fn • g)
      _ ≤ eLpNorm Fn (ENNReal.ofReal p) vol01 * G :=
          eLpNorm_smul_le_mul_eLpNorm (measurable_approx N c).aestronglyMeasurable
            hgint.aestronglyMeasurable
  rw [hFI] at hH
  have key : ‖1 / w‖ₑ < ‖1 / w‖ₑ := by
    calc ‖1 / w‖ₑ ≤ eLpNorm Fn (ENNReal.ofReal p) vol01 * G := hH
      _ ≤ ENNReal.ofReal (c0 / (G.toReal + 1)) * G := by gcongr
      _ = ENNReal.ofReal (c0 / (G.toReal + 1) * G.toReal) := by
          rw [ENNReal.ofReal_mul (by positivity), ENNReal.ofReal_toReal hGne]
      _ < ENNReal.ofReal c0 := by
          rw [ENNReal.ofReal_lt_ofReal_iff hc0pos, div_mul_eq_mul_div,
            div_lt_iff₀ (by positivity)]
          have := ENNReal.toReal_nonneg (a := G)
          nlinarith
      _ = ‖1 / w‖ₑ := by rw [hc0, ofReal_norm]
  exact lt_irrefl _ key

/-- (ii) ⇒ (i) of Theorem A. -/
theorem zeroFree_of_chiInW {p : ℝ} (hp : 1 < p) (h : ChiInW p) : ZeroFree (1 / p) :=
  fun _ hs hz => propN hp hz hs h

end NBLp

#print axioms NBLp.propN


