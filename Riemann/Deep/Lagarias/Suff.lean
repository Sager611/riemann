import Riemann.Tree.R009
import Riemann.Tree.R007

/-!
# Lagarias positivity ⇒ RH (sufficiency direction)

Lagarias (1999): RH ⟺ Re ξ'(s)/ξ(s) > 0 on Re s > 1/2. Here: the direction ⇐.

* `xi_eq`: ξ(s) = s(s-1)/2 · Λ(s) for s ≠ 0, 1 (Mathlib's Λ = Λ₀ - 1/s - 1/(1-s)).
* `xi_zero_iff`: on Re s > 0, s ≠ 1, ξ(s) = 0 ↔ ζ(s) = 0 (ζ = Λ / Γℝ, Γℝ ≠ 0).
* `not_LagariasPositivity_of_zero`: a zero ρ of order m ≥ 1 with Re ρ > 1/2 makes
  Re ξ'/ξ(ρ - δ) ≈ -m/δ < 0 at nearby non-zeros, refuting positivity honestly.
* `RH_of_LagariasPositivity`: conclude via `Fayed.R7.RH_iff_QuasiRH_half`.
-/

open Complex

namespace Fayed.Lagarias

open Fayed.R9

/-- ξ(s) = s(s-1)/2 · Λ(s) away from the poles of Λ. -/
theorem xi_eq (s : ℂ) (h0 : s ≠ 0) (h1 : s ≠ 1) :
    xi s = s * (s - 1) / 2 * completedRiemannZeta s := by
  have h1' : (1 : ℂ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  rw [xi, completedRiemannZeta_eq]
  field_simp
  ring

/-- On Re s > 0, s ≠ 1: zeros of ξ are exactly zeros of ζ. -/
theorem xi_zero_iff {s : ℂ} (hs : 0 < s.re) (h1 : s ≠ 1) :
    xi s = 0 ↔ riemannZeta s = 0 := by
  have h0 : s ≠ 0 := by rintro rfl; simp at hs
  have h1' : s - 1 ≠ 0 := sub_ne_zero.mpr h1
  have hG : Gammaℝ s ≠ 0 := Gammaℝ_ne_zero_of_re_pos hs
  rw [xi_eq s h0 h1, riemannZeta_def_of_ne_zero h0, div_eq_zero_iff, or_iff_left hG]
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h | h
    · exfalso
      rcases div_eq_zero_iff.mp h with h | h
      · exact (mul_ne_zero h0 h1') h
      · exact two_ne_zero h
    · exact h
  · intro h; rw [h, mul_zero]

/-- ξ is entire. -/
theorem xi_differentiable : Differentiable ℂ xi := by
  unfold xi
  have := differentiable_completedZeta₀
  fun_prop

/-- ξ(2) ≠ 0, so ξ is not identically zero. -/
theorem xi_two_ne_zero : xi 2 ≠ 0 := by
  rw [Ne, xi_zero_iff (by norm_num) (by norm_num)]
  exact riemannZeta_ne_zero_of_one_lt_re (by norm_num)

/-- A zero ρ of ξ with Re ρ > 1/2 violates Lagarias positivity at genuine non-zeros of ξ:
writing ξ(s) = (s - ρ)^m g(s) with m ≥ 1, g(ρ) ≠ 0, one has ξ'/ξ = m/(s - ρ) + g'/g, and at
s = ρ - δ the term m/(s - ρ) = -m/δ dominates the bounded g'/g for small real δ > 0. -/
theorem not_LagariasPositivity_of_zero {ρ : ℂ} (hρ : 1 / 2 < ρ.re) (hxi : xi ρ = 0) :
    ¬ LagariasPositivity := by
  intro hpos
  have han : ∀ z, AnalyticAt ℂ xi z := fun z => xi_differentiable.analyticAt z
  have h2 : analyticOrderAt xi 2 ≠ ⊤ := by
    rw [(han 2).analyticOrderAt_eq_zero.mpr xi_two_ne_zero]; exact ENat.zero_ne_top
  have htop : analyticOrderAt xi ρ ≠ ⊤ :=
    AnalyticOnNhd.analyticOrderAt_ne_top_of_isPreconnected (fun z _ => han z) isPreconnected_univ
      (Set.mem_univ 2) (Set.mem_univ ρ) h2
  obtain ⟨m, hm⟩ := ENat.ne_top_iff_exists.mp htop
  have hm0 : m ≠ 0 := by
    rintro rfl
    exact ((han ρ).analyticOrderAt_ne_zero.mpr hxi) (by rw [← hm]; rfl)
  obtain ⟨g, hg, hgρ, hfac⟩ := (han ρ).analyticOrderAt_eq_natCast.mp hm.symm
  simp only [smul_eq_mul] at hfac
  have hfacEq : xi =ᶠ[nhds ρ] fun z => (z - ρ) ^ m * g z := hfac
  have hLD : logDeriv xi =ᶠ[nhds ρ] logDeriv (fun z => (z - ρ) ^ m * g z) :=
    logDeriv_congr_nhds hfacEq
  have hgne : ∀ᶠ z in nhds ρ, g z ≠ 0 := hg.continuousAt.eventually_ne hgρ
  have hgan : ∀ᶠ z in nhds ρ, AnalyticAt ℂ g z := hg.eventually_analyticAt
  have hre : ∀ᶠ z in nhds ρ, 1 / 2 < z.re :=
    Complex.continuous_re.continuousAt.eventually (lt_mem_nhds hρ)
  set B := (logDeriv g ρ).re + 1 with hBdef
  have hcont : ContinuousAt (fun z => (logDeriv g z).re) ρ := by
    have : ContinuousAt (logDeriv g) ρ := hg.deriv.continuousAt.div hg.continuousAt hgρ
    exact Complex.continuous_re.continuousAt.comp this
  have hB : ∀ᶠ z in nhds ρ, (logDeriv g z).re < B :=
    hcont.eventually (gt_mem_nhds (by rw [hBdef]; linarith))
  have hpath : Filter.Tendsto (fun δ : ℝ => ρ - (δ : ℂ)) (nhdsWithin 0 (Set.Ioi 0)) (nhds ρ) := by
    have hc : Continuous (fun δ : ℝ => ρ - (δ : ℂ)) := by fun_prop
    have := (hc.tendsto 0).mono_left (nhdsWithin_le_nhds (s := Set.Ioi (0 : ℝ)))
    simpa using this
  have hmpos : (0 : ℝ) < m := by exact_mod_cast Nat.pos_of_ne_zero hm0
  have hsmall : ∀ᶠ δ : ℝ in nhdsWithin 0 (Set.Ioi 0), δ * B < m := by
    have ht : Filter.Tendsto (fun δ : ℝ => δ * B) (nhds 0) (nhds 0) := by
      have hc : Continuous fun δ : ℝ => δ * B := by fun_prop
      simpa using hc.tendsto 0
    exact (ht.eventually (gt_mem_nhds hmpos)).filter_mono nhdsWithin_le_nhds
  obtain ⟨δ, ⟨hδ, hsδ⟩, hL, hfs, hgs, hgas, hres, hBs⟩ :=
    ((eventually_mem_nhdsWithin.and hsmall).and
      (hpath.eventually (hLD.and (hfacEq.and (hgne.and (hgan.and (hre.and hB))))))).exists
  have hδ : (0 : ℝ) < δ := hδ
  have hsub : (ρ - (δ : ℂ)) - ρ = -(δ : ℂ) := by ring
  have hδne : (δ : ℂ) ≠ 0 := by exact_mod_cast hδ.ne'
  have hpne : ((ρ - (δ : ℂ)) - ρ) ^ m ≠ 0 :=
    pow_ne_zero _ (by rw [hsub]; exact neg_ne_zero.mpr hδne)
  have hxis : xi (ρ - (δ : ℂ)) ≠ 0 := by
    rw [hfs]; exact mul_ne_zero hpne hgs
  have hcalc : logDeriv xi (ρ - (δ : ℂ)) = (m : ℂ) * (1 / -(δ : ℂ)) + logDeriv g (ρ - (δ : ℂ)) := by
    rw [hL, logDeriv_fun_mul (f := fun z => (z - ρ) ^ m) _ hpne hgs (by fun_prop)
      hgas.differentiableAt, logDeriv_fun_pow (by fun_prop)]
    rw [logDeriv_apply (fun z => z - ρ), deriv_sub_const, deriv_id'', hsub]
  have hpos' := hpos (ρ - (δ : ℂ)) hres hxis
  change 0 < (logDeriv xi (ρ - (δ : ℂ))).re at hpos'
  have hre' : ((m : ℂ) * (1 / -(δ : ℂ))).re = -((m : ℝ) / δ) := by
    rw [show (m : ℂ) * (1 / -(δ : ℂ)) = ((-((m : ℝ) / δ) : ℝ) : ℂ) by push_cast; ring]
    exact Complex.ofReal_re _
  rw [hcalc, Complex.add_re, hre'] at hpos'
  have e : δ * (-((m : ℝ) / δ) + (logDeriv g (ρ - (δ : ℂ))).re)
      = -(m : ℝ) + δ * (logDeriv g (ρ - (δ : ℂ))).re := by
    field_simp
  have h1 := mul_pos hδ hpos'
  have h3 := mul_lt_mul_of_pos_left hBs hδ
  linarith

/-- Lagarias positivity implies RH.

`LagariasPositivity` only constrains ξ'/ξ at non-zeros of ξ, so a zero ρ (Re ρ > 1/2) is
refuted through genuine values: `not_LagariasPositivity_of_zero` uses the simple pole of ξ'/ξ
with residue m = ord_ρ ξ ≥ 1 and the point s = ρ - δ, where Re ξ'/ξ ≈ -m/δ < 0. -/
theorem RH_of_LagariasPositivity (h : LagariasPositivity) : RiemannHypothesis := by
  refine Fayed.R7.RH_iff_QuasiRH_half.mpr ?_
  intro s hs hs1 hz
  have hs' : 1 / 2 < s.re := by linarith
  have hxi : xi s = 0 := (xi_zero_iff (by linarith) hs1).mpr hz
  exact not_LagariasPositivity_of_zero hs' hxi h

end Fayed.Lagarias

#print axioms Fayed.Lagarias.xi_eq
#print axioms Fayed.Lagarias.xi_zero_iff
#print axioms Fayed.Lagarias.not_LagariasPositivity_of_zero
#print axioms Fayed.Lagarias.RH_of_LagariasPositivity
