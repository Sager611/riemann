import Riemann.Deep.WeilEF.Z23.FromPNTPlus.ZetaBounds

/-!
# A uniform bound for ζ on 1/2 ≤ Re s ≤ 3

‖ζ(σ + it)‖ ≤ 13 (1 + log|t|) |t|^{max(1-σ,0)} for 1/2 ≤ σ ≤ 3, |t| ≥ 1.

From the first-order Euler–Maclaurin formula `Zeta0EqZeta` (ported from PNT+) at N = ⌊|t|⌋:
ζ(s) = Σ_{n ≤ N} n^{-s} - N^{1-s}/(1-s) - N^{-s}/2 + s ∫_N^∞ (⌊x⌋+1/2-x) x^{-s-1} dx.
-/

open Complex Set MeasureTheory

namespace Fayed.InputsProof.PerronAux

theorem norm_zeta_le {σ t : ℝ} (h1 : 1 / 2 ≤ σ) (h2 : σ ≤ 3) (ht : 1 ≤ |t|) :
    ‖riemannZeta (σ + t * I)‖ ≤ 13 * (1 + Real.log |t|) * |t| ^ (max (1 - σ) 0) := by
  set N := ⌊|t|⌋₊ with hN
  have hNpos : 0 < N := Nat.floor_pos.mpr ht
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hNpos
  have hNt : (N : ℝ) ≤ |t| := Nat.floor_le (abs_nonneg _)
  have htN : |t| < N + 1 := Nat.lt_floor_add_one _
  have hσpos : 0 < σ := by linarith
  have hsre : ((σ : ℂ) + t * I).re = σ := by simp
  have hsim : ((σ : ℂ) + t * I).im = t := by simp
  have hs1 : (σ : ℂ) + t * I ≠ 1 := by
    intro h
    have := congrArg Complex.im h
    rw [hsim] at this
    simp only [one_im] at this
    rw [this] at ht
    simp at ht
    linarith
  set m := max (1 - σ) 0 with hm
  have hm0 : 0 ≤ m := le_max_right _ _
  have hNm : (1 : ℝ) ≤ (N : ℝ) ^ m := Real.one_le_rpow hN1 hm0
  have hpow : ∀ x : ℝ, 1 ≤ x → x ≤ N → x ^ (1 - σ) ≤ (N : ℝ) ^ m := fun x hx hxN =>
    (Real.rpow_le_rpow_of_exponent_le hx (le_max_left _ _)).trans
      (Real.rpow_le_rpow (by linarith) hxN hm0)
  have hlogN : 0 ≤ Real.log N := Real.log_nonneg hN1
  rw [← Zeta0EqZeta hNpos (by rw [hsre]; exact hσpos) hs1, riemannZeta0]
  -- (a) the finite sum
  have ha : ‖∑ n ∈ Finset.range (N + 1), 1 / (n : ℂ) ^ ((σ : ℂ) + t * I)‖ ≤
      (1 + Real.log N) * (N : ℝ) ^ m := by
    refine (norm_sum_le _ _).trans ?_
    have hb : ∀ n ∈ Finset.range (N + 1),
        ‖1 / (n : ℂ) ^ ((σ : ℂ) + t * I)‖ ≤ (n : ℝ)⁻¹ * (N : ℝ) ^ m := by
      intro n hn
      rcases Nat.eq_zero_or_pos n with h0 | hpos
      · subst h0
        have hne : (σ : ℂ) + t * I ≠ 0 := by
          intro h; have := congrArg Complex.re h; rw [hsre] at this; simp at this; linarith
        simp [Complex.zero_cpow hne]
      · have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hpos
        have hnN : (n : ℝ) ≤ N := by
          exact_mod_cast Nat.lt_succ_iff.mp (Finset.mem_range.mp hn)
        rw [norm_div, norm_one, norm_natCast_cpow_of_pos hpos, hsre]
        have e : (1 : ℝ) / (n : ℝ) ^ σ = (n : ℝ)⁻¹ * (n : ℝ) ^ (1 - σ) := by
          rw [Real.rpow_sub (by linarith), Real.rpow_one]
          field_simp
        rw [e]
        exact mul_le_mul_of_nonneg_left (hpow n hn1 hnN) (by positivity)
    refine (Finset.sum_le_sum hb).trans ?_
    rw [← Finset.sum_mul]
    refine mul_le_mul_of_nonneg_right ?_ (by positivity)
    have hr : Finset.range (N + 1) = Finset.Icc 0 N := by
      ext k; simp
    rw [hr, harmonic_eq_sum_Icc0]
    exact harmonic_le_one_add_log N
  -- (b)
  have hb : ‖(-(N : ℂ) ^ (1 - ((σ : ℂ) + t * I))) / (1 - ((σ : ℂ) + t * I))‖ ≤ (N : ℝ) ^ m := by
    rw [norm_div, norm_neg, norm_natCast_cpow_of_pos hNpos]
    have hre : (1 - ((σ : ℂ) + t * I)).re = 1 - σ := by simp
    rw [hre]
    have h1s : 1 ≤ ‖1 - ((σ : ℂ) + t * I)‖ := by
      refine le_trans ?_ (abs_im_le_norm _)
      simpa using ht
    calc (N : ℝ) ^ (1 - σ) / ‖1 - ((σ : ℂ) + t * I)‖ ≤ (N : ℝ) ^ (1 - σ) / 1 :=
          div_le_div_of_nonneg_left (by positivity) zero_lt_one h1s
      _ ≤ (N : ℝ) ^ m := by rw [div_one]; exact hpow N hN1 le_rfl
  -- (c)
  have hc : ‖(-(N : ℂ) ^ (-((σ : ℂ) + t * I))) / 2‖ ≤ (N : ℝ) ^ m := by
    rw [norm_div, norm_neg, norm_natCast_cpow_of_pos hNpos]
    have hre : (-((σ : ℂ) + t * I)).re = -σ := by simp
    rw [hre]
    have : (N : ℝ) ^ (-σ) ≤ 1 := Real.rpow_le_one_of_one_le_of_nonpos hN1 (by linarith)
    simp only [Complex.norm_ofNat]
    linarith
  -- (d)
  have hd : ‖((σ : ℂ) + t * I) *
      ∫ x in Ioi (N : ℝ), (⌊x⌋ + 1 / 2 - x) / (x : ℂ) ^ (((σ : ℂ) + t * I) + 1)‖ ≤
      10 * (N : ℝ) ^ m := by
    rw [norm_mul]
    have hI := ZetaBnd_aux1b N (by exact_mod_cast hNpos) (t := t) hσpos
    have hsn : ‖(σ : ℂ) + t * I‖ ≤ 5 * N := by
      refine (norm_add_le _ _).trans ?_
      simp only [Complex.norm_real, Real.norm_eq_abs, norm_mul, Complex.norm_I, mul_one,
        abs_of_pos hσpos]
      linarith
    have hNs : (N : ℝ) ^ (-σ) / σ ≤ 2 * (N : ℝ) ^ (-σ) := by
      rw [div_le_iff₀ hσpos]
      have : 0 ≤ (N : ℝ) ^ (-σ) := by positivity
      nlinarith
    have hNN : (N : ℝ) * (N : ℝ) ^ (-σ) = (N : ℝ) ^ (1 - σ) := by
      rw [Real.rpow_sub (by linarith), Real.rpow_one, Real.rpow_neg (by linarith)]
      field_simp
    calc ‖(σ : ℂ) + t * I‖ * ‖∫ x in Ioi (N : ℝ), (⌊x⌋ + 1 / 2 - x) /
          (x : ℂ) ^ (((σ : ℂ) + t * I) + 1)‖
        ≤ (5 * N) * (2 * (N : ℝ) ^ (-σ)) :=
          mul_le_mul hsn (hI.trans hNs) (norm_nonneg _) (by positivity)
      _ = 10 * (N : ℝ) ^ (1 - σ) := by rw [← hNN]; ring
      _ ≤ 10 * (N : ℝ) ^ m := mul_le_mul_of_nonneg_left (hpow N hN1 le_rfl) (by norm_num)
  have hsum := norm_add₄_le.trans (add_le_add (add_le_add (add_le_add ha hb) hc) hd)
  refine hsum.trans ?_
  have hlog : Real.log N ≤ Real.log |t| := Real.log_le_log (by linarith) hNt
  have hNmt : (N : ℝ) ^ m ≤ |t| ^ m := Real.rpow_le_rpow (by linarith) hNt hm0
  have hlog0 : 0 ≤ Real.log |t| := Real.log_nonneg ht
  calc (1 + Real.log N) * (N : ℝ) ^ m + (N : ℝ) ^ m + (N : ℝ) ^ m + 10 * (N : ℝ) ^ m
      ≤ 13 * (1 + Real.log N) * (N : ℝ) ^ m := by nlinarith
    _ ≤ 13 * (1 + Real.log |t|) * |t| ^ m := by
        gcongr

end Fayed.InputsProof.PerronAux
