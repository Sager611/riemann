import Riemann.Deep.FirstZero.MinArg

/-!
# Minimality of the first zero, part 3: from certified arg bounds to "no zero below 14.13"

no_zero_below_of_cert: if the arg change of ζ along 2 → 2+iT → σL+iT → σL (σL < 1/2, T ≥ 14.14)
is certified < 3π, then ζ has no zero with 0 < Re < 1, 0 < Im < 14.13. Proof: the folded
argument principle gives Σ ord − 1 < 3, while a hypothetical zero ρ (or its reflection 1 − ρ̄),
its conjugate, the certified first zero 1/2 + it₁ and its conjugate are four distinct zeros in
the rectangle, each of order ≥ 1.
-/

open Complex Set

namespace Fayed.FirstZero.Min
open Fayed.FirstZero

lemma one_le_ord_of_zero {ρ : ℂ} (hz : riemannZeta ρ = 0) (h0 : 0 < ρ.re) (h1 : ρ.re < 1) :
    1 ≤ analyticOrderNatAt riemannZeta ρ :=
  Zeta23.ZetaSeam.one_le_mult_holds ρ ⟨hz, h0, h1⟩

lemma one_le_ord_conj {ρ : ℂ} (hz : riemannZeta ρ = 0) (h0 : 0 < ρ.re) (h1 : ρ.re < 1) :
    1 ≤ analyticOrderNatAt riemannZeta ((starRingEnd ℂ) ρ) := by
  have hρ1 : ρ ≠ 1 := by intro h; rw [h] at h1; simp at h1
  have := one_le_ord_of_zero hz h0 h1
  unfold analyticOrderNatAt at this ⊢
  rw [Zeta23.analyticOrderAt_zeta_conj hρ1]
  exact this

theorem unique_of_cert {σL T Utop Ur Uleft : ℝ} (hσ0 : 0 < σL) (hσh : σL < 1 / 2)
    (hT : 14.14 ≤ T) (htop : TopOK T σL 2 Utop) (hleft : LeftOK σL 0 T Uleft)
    (hr : arg (riemannZeta ((2 : ℝ) + T * I)) ≤ Ur) (hsum : Utop + Ur + Uleft < 3 * Real.pi) :
    ∃ t1 : ℝ, 14.13 < t1 ∧ t1 < 14.14 ∧ riemannZeta (1 / 2 + t1 * I) = 0 ∧
      ∀ s : ℂ, riemannZeta s = 0 → 0 < s.re → s.re < 1 → 0 < s.im → s.im ≤ T →
        s = 1 / 2 + t1 * I := by
  classical
  obtain ⟨t1, ht1a, ht1b, hz1⟩ := first_zero
  refine ⟨t1, ht1a, ht1b, hz1, fun s hz h0 h1 hi0 hi1 => ?_⟩
  by_contra hne
  have hT0 : 0 < T := by linarith
  obtain ⟨Z, hZ, heq⟩ := winding hσ0 (by linarith) hT0 hleft.1 htop.1
  -- upper bound on the count
  have hR := right_edge (T := T) hT0.le
  have hup : ((∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ)) - 1) * Real.pi < 3 * Real.pi := by
    rw [heq, hR]; linarith [htop.2, hleft.2]
  have hcount : (∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ)) < 4 := by
    have := Real.pi_pos
    by_contra hc; push Not at hc; nlinarith
  -- four zeros
  set c : ℂ := 1 / 2 + (t1 : ℂ) * I with hc
  have hcre : c.re = 1 / 2 := by simp [c]
  have hcim : c.im = t1 := by simp [c]
  set a : ℂ := if σL ≤ s.re then s else Zeta23.reflect s with ha
  have hare : σL ≤ a.re ∧ a.re < 1 ∧ 0 < a.re := by
    simp only [a]; split_ifs with h
    · exact ⟨h, h1, h0⟩
    · simp only [Zeta23.reflect, sub_re, one_re, conj_re]; push Not at h
      exact ⟨by linarith, by linarith, by linarith⟩
  have haim : a.im = s.im := by
    simp only [a]; split_ifs
    · rfl
    · simp [Zeta23.reflect]
  have haz : riemannZeta a = 0 := by
    simp only [a]; split_ifs
    · exact hz
    · exact (Zeta23.zeta_reflect_zero s ⟨hz, h0, h1⟩).1
  have hrect : ∀ x : ℂ, σL ≤ x.re → x.re ≤ 2 → -T ≤ x.im → x.im ≤ T →
      x ∈ Rectangle ((σL : ℂ) - T * I) (2 + T * I) := by
    intro x h1 h2 h3 h4
    simp only [Rectangle, mem_reProdIm, sub_re, ofReal_re, mul_re, I_re, I_im, ofReal_im,
      add_re, sub_im, add_im, mul_im, re_ofNat, im_ofNat]
    rw [uIcc_of_le (by norm_num; linarith), uIcc_of_le (by norm_num; linarith)]
    constructor <;> constructor <;> simp <;> linarith
  have hne1 : ∀ x : ℂ, x.im ≠ 0 → x ≠ 1 := by
    intro x hx h; rw [h] at hx; simp at hx
  have hzc : ∀ x : ℂ, riemannZeta x = 0 → riemannZeta ((starRingEnd ℂ) x) = 0 := by
    intro x hx; rw [riemannZeta_conj, hx, map_zero]
  have mem : ∀ x : ℂ, riemannZeta x = 0 → σL ≤ x.re → x.re ≤ 2 → 0 < |x.im| → |x.im| ≤ T →
      x ∈ Z ∧ (starRingEnd ℂ) x ∈ Z := by
    intro x hx hr1 hr2 hi hi2
    rw [abs_le] at hi2
    have hx0 : x.im ≠ 0 := by intro h; rw [h] at hi; simp at hi
    refine ⟨hZ x (hrect x hr1 hr2 hi2.1 hi2.2) (hne1 x hx0) hx, hZ _ ?_ (hne1 _ ?_) (hzc x hx)⟩
    · exact hrect _ (by simpa using hr1) (by simpa using hr2) (by simp; linarith [hi2.2])
        (by simp; linarith [hi2.1])
    · simpa using hx0
  have hA := mem a haz hare.1 (by linarith [hare.2.1]) (by rw [haim, abs_of_pos hi0]; exact hi0)
    (by rw [haim, abs_of_pos hi0]; exact hi1)
  have hac : a ≠ c := by
    simp only [a]; split_ifs with h
    · exact hne
    · intro e
      have := congrArg Complex.re e
      simp only [Zeta23.reflect, sub_re, one_re, conj_re] at this
      rw [hcre] at this; push Not at h; linarith
  have hC := mem c hz1 (by rw [hcre]; linarith) (by rw [hcre]; norm_num)
    (by rw [hcim, abs_of_pos (by linarith)]; linarith) (by rw [hcim, abs_of_pos (by linarith)]; linarith)
  -- distinctness by imaginary parts
  have hb_im : ((starRingEnd ℂ) a).im = -s.im := by simp [haim]
  have hd_im : ((starRingEnd ℂ) c).im = -t1 := by simp [hcim]
  set S : Finset ℂ := {a, (starRingEnd ℂ) a, c, (starRingEnd ℂ) c} with hS
  have hSZ : S ⊆ Z := by
    intro x hx
    simp only [S, Finset.mem_insert, Finset.mem_singleton] at hx
    rcases hx with rfl | rfl | rfl | rfl
    exacts [hA.1, hA.2, hC.1, hC.2]
  have hne' : ∀ x y : ℂ, x.im ≠ y.im → x ≠ y := fun x y h e => h (by rw [e])
  have n1 : a ∉ ({(starRingEnd ℂ) a, c, (starRingEnd ℂ) c} : Finset ℂ) := by
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or]
    refine ⟨hne' _ _ ?_, hac, hne' _ _ ?_⟩
    · rw [haim, hb_im]; linarith
    · rw [haim, hd_im]; linarith
  have n2 : (starRingEnd ℂ) a ∉ ({c, (starRingEnd ℂ) c} : Finset ℂ) := by
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or]
    refine ⟨hne' _ _ ?_, fun e => hac ?_⟩
    · rw [hb_im, hcim]; linarith
    · have := congrArg (starRingEnd ℂ) e
      simpa using this
  have n3 : c ∉ ({(starRingEnd ℂ) c} : Finset ℂ) := by
    simp only [Finset.mem_singleton]
    exact hne' _ _ (by rw [hcim, hd_im]; linarith)
  have hcre0 : 0 < c.re := by rw [hcre]; norm_num
  have hcre1 : c.re < 1 := by rw [hcre]; norm_num
  have hsumS : (4 : ℝ) ≤ ∑ ρ ∈ S, (analyticOrderNatAt riemannZeta ρ : ℝ) := by
    rw [hS, Finset.sum_insert n1, Finset.sum_insert n2, Finset.sum_insert n3,
      Finset.sum_singleton]
    have e1 : (1 : ℝ) ≤ analyticOrderNatAt riemannZeta a := by
      exact_mod_cast one_le_ord_of_zero haz hare.2.2 hare.2.1
    have e2 : (1 : ℝ) ≤ analyticOrderNatAt riemannZeta ((starRingEnd ℂ) a) := by
      exact_mod_cast one_le_ord_conj haz hare.2.2 hare.2.1
    have e3 : (1 : ℝ) ≤ analyticOrderNatAt riemannZeta c := by
      exact_mod_cast one_le_ord_of_zero hz1 hcre0 hcre1
    have e4 : (1 : ℝ) ≤ analyticOrderNatAt riemannZeta ((starRingEnd ℂ) c) := by
      exact_mod_cast one_le_ord_conj hz1 hcre0 hcre1
    linarith
  have hsub := Finset.sum_le_sum_of_subset_of_nonneg hSZ
    (f := fun ρ => (analyticOrderNatAt riemannZeta ρ : ℝ)) (fun _ _ _ => Nat.cast_nonneg _)
  linarith

theorem no_zero_below_of_cert {σL T Utop Ur Uleft : ℝ} (hσ0 : 0 < σL) (hσh : σL < 1 / 2)
    (hT : 14.14 ≤ T) (htop : TopOK T σL 2 Utop) (hleft : LeftOK σL 0 T Uleft)
    (hr : arg (riemannZeta ((2 : ℝ) + T * I)) ≤ Ur) (hsum : Utop + Ur + Uleft < 3 * Real.pi)
    (s : ℂ) (hz : riemannZeta s = 0) (h0 : 0 < s.re) (h1 : s.re < 1) (hi0 : 0 < s.im)
    (hi1 : s.im < 14.13) : False := by
  obtain ⟨t1, ht1a, _, _, hu⟩ := unique_of_cert hσ0 hσh hT htop hleft hr hsum
  have := congrArg Complex.im (hu s hz h0 h1 hi0 (by linarith))
  simp at this
  linarith

end Fayed.FirstZero.Min
