import Riemann.Deep.FirstZero.Min

/-!
# The first zeta zero is simple

Fayed.FirstZero.first_zero_simple: the certified zero 1/2 + it₁ (14.13 < t₁ < 14.14) has
ζ'(1/2 + it₁) ≠ 0, and it is the only zero of ζ in the strip with 0 < Im s ≤ 15.

Route: the folded argument principle of MinArg.winding on [1/4, 2] × [−15, 15], with the arg change
certified in Min.lean (< 3π), gives Σ ord < 4 over the zeros in the box. The zero c and its
conjugate are both in the box and have equal order, so 2·ord c < 4, i.e. ord c = 1. Order 1 of an
analytic function at a zero means nonzero derivative (Mathlib analyticOrderAt_deriv_add_one).
-/

open Complex Set

namespace Fayed.FirstZero.Simple
open Fayed.FirstZero Fayed.FirstZero.Min

/-- An analytic function with a zero of order exactly 1 has nonzero derivative there. -/
lemma deriv_ne_zero_of_order_one {f : ℂ → ℂ} {x : ℂ} (hf : AnalyticAt ℂ f x) (h0 : f x = 0)
    (h1 : analyticOrderAt f x = 1) : deriv f x ≠ 0 := by
  have e := hf.analyticOrderAt_deriv_add_one
  simp only [h0, sub_zero] at e
  rw [h1] at e
  have hd : analyticOrderAt (deriv f) x = 0 := by
    induction h : analyticOrderAt (deriv f) x using ENat.recTopCoe with
    | top => rw [h] at e; simp at e
    | coe n =>
      rw [h] at e
      norm_cast at e ⊢
      omega
  rcases analyticOrderAt_eq_zero.mp hd with h | h
  · exact absurd hf.deriv h
  · exact h

/-- Count-based simplicity: with the certified arg change < 3π, the certified zero has order 1. -/
theorem order_one_of_cert {σL T Utop Ur Uleft : ℝ} (hσ0 : 0 < σL) (hσh : σL < 1 / 2)
    (hT : 14.14 ≤ T) (htop : TopOK T σL 2 Utop) (hleft : LeftOK σL 0 T Uleft)
    (hr : arg (riemannZeta ((2 : ℝ) + T * I)) ≤ Ur) (hsum : Utop + Ur + Uleft < 3 * Real.pi)
    (t1 : ℝ) (ht1a : 14.13 < t1) (ht1b : t1 < 14.14) (hz1 : riemannZeta (1 / 2 + t1 * I) = 0) :
    analyticOrderAt riemannZeta (1 / 2 + t1 * I) = 1 := by
  classical
  have hT0 : 0 < T := by linarith
  obtain ⟨Z, hZ, heq⟩ := winding hσ0 (by linarith) hT0 hleft.1 htop.1
  have hR := right_edge (T := T) hT0.le
  have hup : ((∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ)) - 1) * Real.pi < 3 * Real.pi := by
    rw [heq, hR]; linarith [htop.2, hleft.2]
  have hcount : (∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ)) < 4 := by
    have := Real.pi_pos
    by_contra hc; push Not at hc; nlinarith
  set c : ℂ := 1 / 2 + (t1 : ℂ) * I with hc
  have hcre : c.re = 1 / 2 := by simp [c]
  have hcim : c.im = t1 := by simp [c]
  have hc1 : c ≠ 1 := by intro h; have := congrArg Complex.im h; rw [hcim] at this; simp at this; linarith
  have hcc1 : (starRingEnd ℂ) c ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp [hcim] at this; linarith
  have hrect : ∀ x : ℂ, σL ≤ x.re → x.re ≤ 2 → -T ≤ x.im → x.im ≤ T →
      x ∈ Rectangle ((σL : ℂ) - T * I) (2 + T * I) := by
    intro x h1 h2 h3 h4
    simp only [Rectangle, mem_reProdIm, sub_re, ofReal_re, mul_re, I_re, I_im, ofReal_im,
      add_re, sub_im, add_im, mul_im, re_ofNat, im_ofNat]
    rw [uIcc_of_le (by norm_num; linarith), uIcc_of_le (by norm_num; linarith)]
    constructor <;> constructor <;> simp <;> linarith
  have hcZ : c ∈ Z := hZ c (hrect c (by rw [hcre]; linarith) (by rw [hcre]; norm_num)
    (by rw [hcim]; linarith) (by rw [hcim]; linarith)) hc1 hz1
  have hcZ' : (starRingEnd ℂ) c ∈ Z := hZ _ (hrect _ (by simp [hcre]; linarith)
    (by simp [hcre]; norm_num) (by simp [hcim]; linarith) (by simp [hcim]; linarith)) hcc1
    (by rw [riemannZeta_conj, hz1, map_zero])
  have hne : c ≠ (starRingEnd ℂ) c := by
    intro h; have := congrArg Complex.im h; simp [hcim] at this; linarith
  have hord : analyticOrderNatAt riemannZeta ((starRingEnd ℂ) c) =
      analyticOrderNatAt riemannZeta c := by
    unfold analyticOrderNatAt
    rw [Zeta23.analyticOrderAt_zeta_conj hc1]
  set S : Finset ℂ := {c, (starRingEnd ℂ) c} with hS
  have hSZ : S ⊆ Z := by
    intro x hx
    simp only [S, Finset.mem_insert, Finset.mem_singleton] at hx
    rcases hx with rfl | rfl
    exacts [hcZ, hcZ']
  have hsub := Finset.sum_le_sum_of_subset_of_nonneg hSZ
    (f := fun ρ => (analyticOrderNatAt riemannZeta ρ : ℝ)) (fun _ _ _ => Nat.cast_nonneg _)
  have hsumS : ∑ ρ ∈ S, (analyticOrderNatAt riemannZeta ρ : ℝ) =
      2 * (analyticOrderNatAt riemannZeta c : ℝ) := by
    rw [hS, Finset.sum_insert (by simpa using hne), Finset.sum_singleton, hord]; ring
  have hlt : (analyticOrderNatAt riemannZeta c : ℝ) < 2 := by
    have h2 : ∑ ρ ∈ S, (analyticOrderNatAt riemannZeta ρ : ℝ) ≤
        ∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ) := hsub
    rw [hsumS] at h2; linarith
  have hge : 1 ≤ analyticOrderNatAt riemannZeta c :=
    one_le_ord_of_zero hz1 (by rw [hcre]; norm_num) (by rw [hcre]; norm_num)
  have hnat : analyticOrderNatAt riemannZeta c = 1 := by
    have : analyticOrderNatAt riemannZeta c < 2 := by exact_mod_cast hlt
    omega
  have hnetop := Zeta23.analyticOrderAt_riemannZeta_ne_top hc1
  rw [← Nat.cast_analyticOrderNatAt hnetop, hnat, Nat.cast_one]

end Fayed.FirstZero.Simple

namespace Fayed.FirstZero

/-- **The first zeta zero is simple**, and it is the only zero of ζ in the critical strip with
0 < Im s ≤ 15 (in particular the only zero 1/2 + it with 14.13 < t < 14.14). -/
theorem first_zero_simple_unique : ∃ t : ℝ, 14.13 < t ∧ t < 14.14 ∧
    riemannZeta (1 / 2 + t * I) = 0 ∧ deriv riemannZeta (1 / 2 + t * I) ≠ 0 ∧
    analyticOrderAt riemannZeta (1 / 2 + t * I) = 1 ∧
    ∀ s : ℂ, riemannZeta s = 0 → 0 < s.re → s.re < 1 → 0 < s.im → s.im ≤ 15 →
      s = 1 / 2 + t * I := by
  obtain ⟨t1, ha, hb, hz, hu⟩ := Min.Cert.unique_below
  have hord := Simple.order_one_of_cert (σL := (1 / 4 : ℝ)) (T := (15 : ℝ)) (by norm_num)
    (by norm_num) (by norm_num) Min.Cert.top_all Min.Cert.left_all Min.Cert.right_bd
    (by have := Real.pi_gt_d2; norm_num at this ⊢; linarith) t1 ha hb hz
  have hc1 : (1 / 2 + (t1 : ℂ) * I) ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp at this; linarith
  refine ⟨t1, ha, hb, hz, Simple.deriv_ne_zero_of_order_one
    (Zeta23.riemannZeta_analyticOnNhd_compl_one _ hc1) hz hord, hord, ?_⟩
  intro s h0 h1 h2 h3 h4
  exact hu s h0 h1 h2 h3 (by simpa using h4)

/-- **The first zeta zero is simple.** -/
theorem first_zero_simple : ∃ t : ℝ, 14.13 < t ∧ t < 14.14 ∧
    riemannZeta (1 / 2 + t * I) = 0 ∧ deriv riemannZeta (1 / 2 + t * I) ≠ 0 := by
  obtain ⟨t, ha, hb, hz, hd, -⟩ := first_zero_simple_unique
  exact ⟨t, ha, hb, hz, hd⟩

end Fayed.FirstZero

#print axioms Fayed.FirstZero.first_zero_simple
#print axioms Fayed.FirstZero.first_zero_simple_unique
