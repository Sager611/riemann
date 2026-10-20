import Riemann.Deep.DH.ZeroCertA
import Riemann.Deep.DH.ZeroCertB
import Riemann.Tree.R081

/-!
# Davenport–Heilbronn zero, part 6: the zero (main theorem)

offLineZero_holds : OffLineZeroClaim, i.e. DH has a zero s with |Re s − 0.808517| < 10⁻⁵ and
|Im s − 85.699348| < 10⁻⁵ (Spira 1994), proved by Newton–Kantorovich (newton_zero) around
z_A = 0.8085171825 + 85.6993484854 i:
* ‖DH z_A‖ ≤ 3.7·10⁻⁶ (ball for Gem 60 at z_A from ZeroCertA plus the Euler–Maclaurin error);
* the difference quotient over k = 10⁻⁵ i is within 10⁻⁴ of d = (c_B − c_A)/k, |d| ≥ 1.2323
  (ZeroCertB, plus a Cauchy estimate for the Euler–Maclaurin error);
* |DH| ≤ 157.6 on the disc of radius 0.74 (crude bounds) gives |DH''| ≤ 1151 near z_A.
Hence a zero within 5·10⁻⁶ of z_A.
-/

open Complex Set

namespace Fayed.DH.Zero
open Fayed.DH
open Fayed.FirstZero (nps gc norm_gc_le norm_le_of_sq_le)
open Fayed.FirstZero.Min (rpow_le_q)

local notation "zA" => ((((323406873 / 400000000 : ℝ) : ℝ) : ℝ) : ℂ) + ((((428496742427 / 5000000000 : ℝ) : ℝ) : ℝ) : ℂ) * I
local notation "zB" => ((((323406873 / 400000000 : ℝ) : ℝ) : ℝ) : ℂ) + ((((428496792427 / 5000000000 : ℝ) : ℝ) : ℝ) : ℂ) * I
local notation "kk" => (((1 / 100000 : ℝ)) : ℂ) * I

lemma ball_facts {z w : ℂ} {R : ℝ} (h : ‖z - w‖ ≤ R) :
    w.re - R ≤ z.re ∧ w.im - R ≤ z.im ∧ ‖z‖ ≤ ‖w‖ + R := by
  have h1 := Complex.abs_re_le_norm (z - w)
  have h2 := Complex.abs_im_le_norm (z - w)
  simp only [sub_re, sub_im] at h1 h2
  refine ⟨by linarith [neg_abs_le (z.re - w.re)], by linarith [neg_abs_le (z.im - w.im)], ?_⟩
  calc ‖z‖ = ‖(z - w) + w‖ := by ring_nf
    _ ≤ ‖z - w‖ + ‖w‖ := norm_add_le _ _
    _ ≤ ‖w‖ + R := by linarith

lemma zA_re : (zA).re = (323406873 / 400000000 : ℝ) := by simp
lemma zA_im : (zA).im = (428496742427 / 5000000000 : ℝ) := by simp
lemma norm_zA : ‖zA‖ ≤ (85703163 / 1000000 : ℝ) := norm_le_of_sq_le (by norm_num) (by simp; norm_num)
lemma norm_kk : ‖kk‖ = 1 / 100000 := by
  rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Real.norm_eq_abs]; norm_num
lemma hzB : zA + kk = zB := by apply Complex.ext <;> simp <;> norm_num

lemma Kup : κ ≤ (710197609601030740070729581 / 2500000000000000000000000000 : ℝ) + (1 / 5000000000000000000000000000 : ℝ) := by
  have := Kap; rw [abs_le] at this; linarith

lemma hG7 : ‖gc zA 7‖ ≤ ∏ j ∈ Finset.range 7, [(85718429 / 1000000 : ℝ), (85745357 / 1000000 : ℝ), (85783933 / 1000000 : ℝ), (2682317 / 31250 : ℝ), (2684249 / 31250 : ℝ), (4298469 / 50000 : ℝ), (86054351 / 1000000 : ℝ)].getD j 0 := by
  refine norm_gc_le _ _ 7 fun j hj => ?_
  interval_cases j <;> exact norm_le_of_sq_le (by norm_num) (by simp; norm_num)

/-- Euler–Maclaurin error at z_A. -/
theorem EbdA : Ebd 60 zA ≤ (903959 / 250000000000 : ℝ) :=
  Ebd_le (M := 60) (by norm_num) (σ1 := 4 / 5) (by rw [zA_re]; norm_num) (by norm_num)
    (rpow_le_q (n := 5) (by norm_num) 4 5 (by norm_num) (by norm_num) (v := (13797297 / 50000000 : ℝ)) (by norm_num)
      (by norm_num))
    (rpow_le_q (n := ((60 : ℕ) : ℝ)) (by norm_num) 39 5 (by norm_num) (by norm_num) (v := (27005377 / 2000000000000000000000 : ℝ))
      (by norm_num) (by norm_num))
    norm_zA hG7 Kup (by norm_num [Finset.prod_range_succ])

/-- Crude bound on the disc of radius 10⁻⁵ + 2·0.37. -/
theorem DH_crude : ∀ z, ‖z - zA‖ ≤ 1 / 100000 + 2 * (37 / 100) → ‖DH z‖ ≤ (787691 / 5000 : ℝ) := by
  intro z hz
  obtain ⟨hre, him, hn⟩ := ball_facts hz
  rw [zA_re] at hre
  rw [zA_im] at him
  have hS : ‖z‖ ≤ (86443173 / 1000000 : ℝ) := by have := norm_zA; linarith
  have hre0 : 0 < z.re := by norm_num at hre ⊢; linarith
  have hT : (424796692427 / 5000000000 : ℝ) ≤ z.im := by norm_num at him ⊢; linarith
  have hz1 : z ≠ 1 := by
    intro h; rw [h] at hT; norm_num at hT
  have hE := DH_sub_Gem_le (M := 60) (by norm_num) hre0 hz1
  have hEb : Ebd 60 z ≤ (11197281 / 20000000000 : ℝ) :=
    Ebd_le (M := 60) (by norm_num) (σ1 := 0) hre0.le le_rfl (v5 := 1) (by simp)
      (rpow_le_q (n := ((60 : ℕ) : ℝ)) (by norm_num) 7 1 (by norm_num) (by norm_num)
        (v := (1 / 2799360000000 : ℝ)) (by norm_num) (by norm_num))
      hS (norm_gc_le_S hS 7) Kup (by norm_num [Finset.prod_range_succ])
  have hQ : ∀ j : ℕ, 1 ≤ j → j ≤ 4 → ‖Qf z (((60 : ℕ) : ℝ) + j / 5)‖ ≤ (1342592652109 / 1000000000000 : ℝ) := by
    intro j h1 h4
    have hj : (j : ℝ) ≤ 4 := by exact_mod_cast h4
    have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
    refine (norm_Qf_le (xl := 60) (X := 61) (by norm_num) (by push_cast; linarith)
      (by push_cast; linarith) hS (by norm_num) hT (norm_gc_le_S hS 2) (norm_gc_le_S hS 4)
      (norm_gc_le_S hS 6)).trans ?_
    norm_num [Finset.prod_range_succ]
  have hG := norm_Gem_le (M := 60) (by norm_num) hre0.le hz1 hQ
  have hk := Kup
  calc ‖DH z‖ = ‖(DH z - Gem 60 z) + Gem 60 z‖ := by ring_nf
    _ ≤ ‖DH z - Gem 60 z‖ + ‖Gem 60 z‖ := norm_add_le _ _
    _ ≤ (11197281 / 20000000000 : ℝ) + ((60 : ℕ) + (1342592652109 / 1000000000000 : ℝ)) * (2 + 2 * ((710197609601030740070729581 / 2500000000000000000000000000 : ℝ) + (1 / 5000000000000000000000000000 : ℝ))) := by
        gcongr
        · exact hE.trans hEb
        · exact hG.trans (by gcongr)
    _ ≤ (787691 / 5000 : ℝ) := by norm_num

/-- Euler–Maclaurin error near z_A (for the difference quotient). -/
theorem E_small : ∀ z, ‖z - zA‖ ≤ ‖kk‖ + 3 / 10 → ‖DH z - Gem 60 z‖ ≤ (580261 / 20000000000 : ℝ) := by
  intro z hz
  rw [norm_kk] at hz
  obtain ⟨hre, him, hn⟩ := ball_facts hz
  rw [zA_re] at hre
  rw [zA_im] at him
  have hS : ‖z‖ ≤ (86003173 / 1000000 : ℝ) := by have := norm_zA; linarith
  have hre0 : 1 / 2 ≤ z.re := by norm_num at hre ⊢; linarith
  have hz1 : z ≠ 1 := by
    intro h; rw [h] at him; norm_num at him
  refine (DH_sub_Gem_le (M := 60) (by norm_num) (by linarith) hz1).trans ?_
  exact Ebd_le (M := 60) (by norm_num) (σ1 := 1 / 2) hre0 (by norm_num)
    (rpow_le_q (n := 5) (by norm_num) 1 2 (by norm_num) (by norm_num) (v := (4473 / 10000 : ℝ)) (by norm_num)
      (by norm_num))
    (rpow_le_q (n := ((60 : ℕ) : ℝ)) (by norm_num) 15 2 (by norm_num) (by norm_num) (v := (4611749 / 100000000000000000000 : ℝ))
      (by norm_num) (by norm_num))
    hS (norm_gc_le_S hS 7) Kup (by norm_num [Finset.prod_range_succ])

lemma E_diff : ∀ z, ‖z - zA‖ ≤ ‖kk‖ + 3 / 10 → DifferentiableAt ℂ (fun z => DH z - Gem 60 z) z := by
  intro z hz
  rw [norm_kk] at hz
  obtain ⟨-, him, -⟩ := ball_facts hz
  rw [zA_im] at him
  have hz1 : z ≠ 1 := by
    intro h; rw [h] at him; norm_num at him
  exact (differentiable_DH z).sub (differentiableAt_Gem 60 hz1)

lemma dq_bound {F G : ℂ → ℂ} {z0 k c0 c1 : ℂ} {r0 r1 b : ℝ} (hk : k ≠ 0)
    (h0 : ‖G z0 - c0‖ ≤ r0) (h1 : ‖G (z0 + k) - c1‖ ≤ r1)
    (hE : ‖(F (z0 + k) - G (z0 + k)) - (F z0 - G z0)‖ ≤ b) :
    ‖(F (z0 + k) - F z0) / k - (c1 - c0) / k‖ ≤ (r0 + r1 + b) / ‖k‖ := by
  have e : (F (z0 + k) - F z0) / k - (c1 - c0) / k =
      ((G (z0 + k) - c1) - (G z0 - c0) + ((F (z0 + k) - G (z0 + k)) - (F z0 - G z0))) / k := by
    field_simp; ring
  rw [e, norm_div]
  gcongr
  refine (norm_add_le _ _).trans ?_
  have := norm_sub_le (G (z0 + k) - c1) (G z0 - c0)
  linarith

/-- **The Davenport–Heilbronn off-line zero (Spira 1994), certified.** -/
theorem offLineZero_holds : OffLineZeroClaim := by
  have hk : kk ≠ 0 := by
    intro h; have := congrArg norm h; rw [norm_kk, norm_zero] at this; norm_num at this
  set d : ℂ := ((((1210082949655058269 / 500000000000000000000000 : ℝ) : ℝ) : ℂ) + (((770172634371022231 / 62500000000000000000000 : ℝ) : ℝ) : ℂ) * I - ((((4670756690002783 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) + (((-250271350221213 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) * I)) / kk
    with hd
  -- the difference quotient
  have hD : ‖(DH (zA + kk) - DH zA) / kk - d‖ ≤ (96710203 / 1000000000000 : ℝ) := by
    have hB := GemB
    rw [← hzB] at hB
    have hdq := diff_quot_le (E := fun z => DH z - Gem 60 z) (z0 := zA) (k := kk) (ρ := 3 / 10)
      (eE := (580261 / 20000000000 : ℝ)) (by norm_num) E_diff E_small
    have := dq_bound (F := DH) (G := Gem 60) hk GemA hB hdq
    rw [hd]
    refine this.trans ?_
    rw [norm_kk]
    norm_num
  have hdl : (1232301 / 1000000 : ℝ) ≤ ‖d‖ := by
    rw [hd, norm_div, norm_kk, le_div_iff₀ (by norm_num)]
    refine le_trans ?_ (Complex.abs_im_le_norm _)
    simp only [sub_im, add_im, ofReal_im, mul_im, ofReal_re, I_re, I_im, mul_zero, mul_one,
      zero_add, add_zero]
    norm_num [abs_of_pos]
  have hdpos : 0 < ‖d‖ := lt_of_lt_of_le (by norm_num) hdl
  -- value at z_A
  have hF0 : ‖DH zA‖ ≤ (3620759 / 1000000000000 : ℝ) := by
    have h1 := DH_sub_Gem_le (M := 60) (by norm_num) (s := zA) (by rw [zA_re]; norm_num) hzA1
    have h2 := GemA
    have h3 : ‖((((4670756690002783 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) + (((-250271350221213 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) * I)‖ ≤ |(4670756690002783 / 1000000000000000000000000 : ℝ)| + |(-250271350221213 / 1000000000000000000000000 : ℝ)| :=
      (norm_le_abs_re_add_abs_im _).trans (le_of_eq (by simp))
    calc ‖DH zA‖ = ‖(DH zA - Gem 60 zA) + (Gem 60 zA - ((((4670756690002783 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) + (((-250271350221213 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) * I)) +
          ((((4670756690002783 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) + (((-250271350221213 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) * I)‖ := by ring_nf
      _ ≤ ‖DH zA - Gem 60 zA‖ + ‖Gem 60 zA - ((((4670756690002783 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) + (((-250271350221213 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) * I)‖ +
          ‖((((4670756690002783 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) + (((-250271350221213 / 1000000000000000000000000 : ℝ) : ℝ) : ℂ) * I)‖ := norm_add₃_le
      _ ≤ (903959 / 250000000000 : ℝ) + (882403880393 / 5000000000000000000000000000 : ℝ) + (|(4670756690002783 / 1000000000000000000000000 : ℝ)| + |(-250271350221213 / 1000000000000000000000000 : ℝ)|) := by
          gcongr
          exact h1.trans EbdA
      _ ≤ (3620759 / 1000000000000 : ℝ) := by norm_num [abs_of_pos, abs_of_neg]
  -- second derivative
  have hmax : max (1 / 200000 : ℝ) ‖kk‖ = 1 / 100000 := by rw [norm_kk]; norm_num
  have hM2 : ∀ u, ‖u - zA‖ ≤ max (1 / 200000 : ℝ) ‖kk‖ → ‖deriv (deriv DH) u‖ ≤ (787691 / 5000 : ℝ) / (37 / 100) ^ 2 := by
    rw [hmax]; exact deriv2_le differentiable_DH (by norm_num) DH_crude
  have hq : (((1 / 200000 : ℝ) + ‖kk‖) * ((787691 / 5000 : ℝ) / (37 / 100) ^ 2) + (96710203 / 1000000000000 : ℝ)) / ‖d‖ ≤ (14087 / 1000000 : ℝ) := by
    rw [norm_kk]
    refine (div_le_div_of_nonneg_left (by norm_num) (by norm_num) hdl).trans ?_
    norm_num
  have hc : (3620759 / 1000000000000 : ℝ) / ‖d‖ ≤ (1 - (14087 / 1000000 : ℝ)) * (1 / 200000 : ℝ) := by
    refine (div_le_div_of_nonneg_left (by norm_num) (by norm_num) hdl).trans ?_
    norm_num
  obtain ⟨z, hz, hz0⟩ := newton_zero differentiable_DH (by norm_num) hk hM2 hF0 hD hdpos hq
    (by norm_num) hc
  refine ⟨z, hz0, ?_, ?_⟩
  · have h1 := Complex.abs_re_le_norm (z - zA)
    rw [sub_re, zA_re] at h1
    have h2 := h1.trans hz
    rw [abs_le] at h2
    rw [abs_lt]; constructor <;> norm_num at h2 ⊢ <;> linarith [h2.1, h2.2]
  · have h1 := Complex.abs_im_le_norm (z - zA)
    rw [sub_im, zA_im] at h1
    have h2 := h1.trans hz
    rw [abs_le] at h2
    rw [abs_lt]; constructor <;> norm_num at h2 ⊢ <;> linarith [h2.1, h2.2]

/-- Hence DH has a zero off the critical line, unconditionally. -/
theorem dh_offline_holds : ∃ s : ℂ, DH s = 0 ∧ s.re ≠ 1 / 2 :=
  Fayed.R81.dh_offline offLineZero_holds

end Fayed.DH.Zero

#print axioms Fayed.DH.Zero.offLineZero_holds
#print axioms Fayed.DH.Zero.dh_offline_holds
