import Riemann.Deep.Inputs.PerronZeta

/-!
# Contour shift for ζ(s+z) u^z / z

F(z) = ζ(s+z) u^z / z on the rectangle [1-σ-η, 2] × [-T, T], |t| < T. The only pole inside is
z = 1 - s (residue u^{1-s}/(1-s)); z = 0 lies to the left (1-σ-η > 0). So

(1/2π)∫_{-T}^{T} F(2+iy) dy - (1/2π)∫_{-T}^{T} F(1-σ-η+iy) dy - u^{1-s}/(1-s)
  = (1/2π)(i H₋ - i H₊),  H_± = ∫_{1-σ-η}^{2} F(x ± iT) dx,

and with T = u⁸ the horizontal sides are O((1 + log u) u^{E-8}), E = max(1-σ+7η, 2) < 4.
-/

open Complex Set MeasureTheory intervalIntegral Filter Topology
open scoped Real Interval

namespace Fayed.InputsProof.PerronAux

noncomputable section

/-- The contour integrand ζ(s+z) u^z / z. -/
def Fz (s : ℂ) (u : ℝ) (z : ℂ) : ℂ := riemannZeta (s + z) * (u : ℂ) ^ z / z

lemma Fz_differentiableAt (s : ℂ) {u : ℝ} (hu : 0 < u) {z : ℂ} (hz : z ≠ 0) (hz1 : s + z ≠ 1) :
    DifferentiableAt ℂ (Fz s u) z := by
  have hu' : (u : ℂ) ≠ 0 := by exact_mod_cast hu.ne'
  have h1 : DifferentiableAt ℂ (fun z => riemannZeta (s + z)) z :=
    (differentiableAt_riemannZeta hz1).comp z (differentiableAt_const _ |>.add differentiableAt_id)
  have h2 : DifferentiableAt ℂ (fun z : ℂ => (u : ℂ) ^ z) z :=
    DifferentiableAt.const_cpow differentiableAt_id (Or.inl hu')
  exact (h1.mul h2).div differentiableAt_id hz

theorem shift_identity (s : ℂ) {η u T : ℝ} (hσ0 : 0 < s.re) (hη0 : 0 < η)
    (hη2 : η < 1 - s.re) (hu : 0 < u) (hT : |s.im| < T) :
    ((1 / (2 * π) : ℝ) : ℂ) * (∫ y in (-T)..T, Fz s u ((2 : ℝ) + y * I)) -
      ((1 / (2 * π) : ℝ) : ℂ) * (∫ y in (-T)..T, Fz s u (((1 - s.re - η : ℝ) : ℂ) + y * I)) -
      (u : ℂ) ^ (1 - s) / (1 - s) =
    ((1 / (2 * π) : ℝ) : ℂ) * (I * (∫ x in (1 - s.re - η)..2, Fz s u (x + (-T : ℝ) * I)) -
      I * (∫ x in (1 - s.re - η)..2, Fz s u (x + (T : ℝ) * I))) := by
  set a := 1 - s.re - η with ha
  have ha0 : 0 < a := by rw [ha]; linarith
  have ha2 : a < 2 := by rw [ha]; linarith
  have hT0 : 0 < T := lt_of_le_of_lt (abs_nonneg _) hT
  have hu' : (u : ℂ) ≠ 0 := by exact_mod_cast hu.ne'
  set p : ℂ := 1 - s with hp
  have hp0 : p ≠ 0 := by
    intro h; have := congrArg Complex.re h; rw [hp] at this; simp at this; linarith
  have hpre : p.re = 1 - s.re := by rw [hp]; simp
  have hpim : p.im = -s.im := by rw [hp]; simp
  set A : ℂ := (u : ℂ) ^ (1 - s) / (1 - s) with hA
  have hres : RectangleIntegral' (Fz s u) (⟨a, -T⟩ : ℂ) ⟨2, T⟩ = A := by
    refine ResidueTheoremOnRectangleWithSimplePole' (p := p) (A := A)
      (by simp only; linarith) (by simp only; linarith) ?_ ?_ ?_
    · rw [rectangle_mem_nhds_iff]
      simp only [mem_reProdIm]
      rw [uIoo_of_le (by linarith), uIoo_of_le (by linarith), hpre, hpim]
      have h1 := neg_abs_le s.im
      have h2 := le_abs_self s.im
      exact ⟨⟨by rw [ha]; linarith, by linarith⟩, ⟨by linarith, by linarith⟩⟩
    · intro z hz
      have hzre : z.re ∈ [[a, 2]] := hz.1.1
      rw [uIcc_of_le ha2.le] at hzre
      have hz0 : z ≠ 0 := by
        intro h; rw [h] at hzre; simp at hzre; linarith
      have hz1 : s + z ≠ 1 := by
        intro h; apply hz.2; show z = p; rw [hp]; linear_combination h
      exact (Fz_differentiableAt s hu hz0 hz1).differentiableWithinAt
    · -- principal part at p
      set U : Set ℂ := {z | 0 < z.re}
      have hU : U ∈ 𝓝 p :=
        (isOpen_lt continuous_const Complex.continuous_re).mem_nhds (by
          show 0 < p.re; rw [hpre]; linarith)
      have hhol : HolomorphicOn (Fz s u) (U \ {p}) := by
        intro z hz
        have hz0 : z ≠ 0 := by intro h; have := hz.1; rw [h] at this; simp [U] at this
        have hz1 : s + z ≠ 1 := by
          intro h; apply hz.2; show z = p; rw [hp]; linear_combination h
        exact (Fz_differentiableAt s hu hz0 hz1).differentiableWithinAt
      have hshift : Tendsto (fun z => s + z) (𝓝[≠] p) (𝓝[≠] 1) := by
        rw [tendsto_nhdsWithin_iff]
        refine ⟨?_, ?_⟩
        · have hc : Continuous (fun z : ℂ => s + z) := continuous_const.add continuous_id
          have := hc.tendsto p
          rw [show s + p = 1 by rw [hp]; ring] at this
          exact this.mono_left nhdsWithin_le_nhds
        · filter_upwards [self_mem_nhdsWithin] with z hz
          intro h
          apply hz
          show z = p
          rw [hp]
          have : s + z = 1 := h
          linear_combination this
      have hζ := riemannZeta_residue_one.comp hshift
      have hcont : ContinuousAt (fun z : ℂ => (u : ℂ) ^ z / z) p :=
        ((DifferentiableAt.const_cpow differentiableAt_id (Or.inl hu')).div
          differentiableAt_id hp0).continuousAt
      have hlim : Tendsto (fun z => (z - p) * Fz s u z) (𝓝[≠] p)
          (𝓝 (1 * ((u : ℂ) ^ p / p))) := by
        refine (hζ.mul (hcont.tendsto.mono_left nhdsWithin_le_nhds)).congr (fun z => ?_)
        simp only [Function.comp, Fz, hp]
        ring
      have hA' : 1 * ((u : ℂ) ^ p / p) = A := by rw [one_mul, hA, hp]
      rw [hA'] at hlim
      obtain ⟨V, hV, hb⟩ := ResidueOfTendsTo hU hhol hlim
      have hO := BddAbove_to_IsBigO hV hb
      have e : (Fz s u - fun z => A / (z - p)) = (Fz s u - fun z => A * (z - p)⁻¹) := by
        funext z; simp [div_eq_mul_inv]
      rw [e]
      exact hO
  simp only [RectangleIntegral', RectangleIntegral, HIntegral, VIntegral, smul_eq_mul] at hres
  set Hb := ∫ x in a..2, Fz s u (x + (-T : ℝ) * I)
  set Ht := ∫ x in a..2, Fz s u (x + (T : ℝ) * I)
  set V2 := ∫ y in (-T)..T, Fz s u ((2 : ℝ) + y * I)
  set Va := ∫ y in (-T)..T, Fz s u ((a : ℂ) + y * I)
  have hπ' : (π : ℂ) ≠ 0 := by exact_mod_cast Real.pi_pos.ne'
  have h' : Hb - Ht + I * V2 - I * Va = 2 * π * I * A := by
    field_simp at hres
    linear_combination hres
  push_cast
  field_simp
  linear_combination (-I) * h' + (V2 - Va) * I_sq + (-2 * (π : ℂ) * A) * I_sq

/-- 7η < 3 + σ when η < σ and η < 1 - σ. -/
lemma seven_eta {σ η : ℝ} (h1 : η < σ) (h2 : η < 1 - σ) : 1 - σ + 7 * η < 4 := by linarith

theorem horiz_bound (s : ℂ) {η u : ℝ} (hη0 : 0 < η)
    (hη1 : η < s.re) (hη2 : η < 1 - s.re) (hu : 1 ≤ u) (hT : |s.im| + 1 ≤ u ^ 8)
    {y : ℝ} (hy : |y| = u ^ 8) :
    ‖∫ x in (1 - s.re - η)..2, Fz s u (x + y * I)‖ ≤
      52 * (1 + Real.log (2 * u ^ 8)) * u ^ (max (1 - s.re + 7 * η) 2) / u ^ 8 := by
  set a := 1 - s.re - η with ha
  set T := u ^ 8 with hTdef
  set E := max (1 - s.re + 7 * η) 2 with hE
  have hu0 : 0 < u := by linarith
  have hT0 : 0 < T := by positivity
  have ha0 : 0 < a := by rw [ha]; linarith
  have ha2 : a < 2 := by rw [ha]; linarith
  have hlog2T : 0 ≤ Real.log (2 * T) := Real.log_nonneg (by linarith [abs_nonneg s.im])
  set K := 26 * (1 + Real.log (2 * T)) * u ^ E / T with hK
  have hpt : ∀ x ∈ Ι a 2, ‖Fz s u (x + y * I)‖ ≤ K := by
    intro x hx
    rw [uIoc_of_le ha2.le] at hx
    obtain ⟨hxa, hx2⟩ := hx
    have hw : s + ((x : ℂ) + y * I) = ((s.re + x : ℝ) : ℂ) + ((s.im + y : ℝ) : ℂ) * I := by
      apply Complex.ext <;> simp
    have hre1 : 1 / 2 ≤ s.re + x := by rw [ha] at hxa; linarith
    have hre2 : s.re + x ≤ 3 := by linarith
    have hty1 : 1 ≤ |s.im + y| := by
      have := abs_sub_abs_le_abs_add y s.im
      rw [hy] at this
      rw [add_comm]
      linarith
    have hty2 : |s.im + y| ≤ 2 * T := by
      have := abs_add_le s.im y
      rw [hy] at this
      linarith [le_abs_self s.im, abs_nonneg s.im]
    have hz := norm_zeta_le hre1 hre2 hty1
    rw [← hw] at hz
    set m := max (1 - (s.re + x)) 0 with hm
    have hm0 : 0 ≤ m := le_max_right _ _
    have hm1 : m ≤ 1 := max_le (by linarith) zero_le_one
    have hzeta : ‖riemannZeta (s + ((x : ℂ) + y * I))‖ ≤
        13 * (1 + Real.log (2 * T)) * (2 * u ^ (8 * m)) := by
      refine hz.trans ?_
      have hl : Real.log |s.im + y| ≤ Real.log (2 * T) :=
        Real.log_le_log (by linarith) hty2
      have hp : |s.im + y| ^ m ≤ (2 * T) ^ m := Real.rpow_le_rpow (by linarith) hty2 hm0
      have hp2 : (2 * T) ^ m ≤ 2 * u ^ (8 * m) := by
        rw [Real.mul_rpow (by norm_num) hT0.le]
        have h2m : (2 : ℝ) ^ m ≤ 2 := by
          calc (2 : ℝ) ^ m ≤ 2 ^ (1 : ℝ) :=
                Real.rpow_le_rpow_of_exponent_le (by norm_num) hm1
            _ = 2 := Real.rpow_one 2
        have hTm : T ^ m = u ^ (8 * m) := by
          rw [hTdef, ← Real.rpow_natCast, ← Real.rpow_mul hu0.le]
          norm_num
        rw [hTm]
        have : 0 ≤ u ^ (8 * m) := by positivity
        nlinarith
      have hlog0 : 0 ≤ Real.log |s.im + y| := Real.log_nonneg hty1
      calc 13 * (1 + Real.log |s.im + y|) * |s.im + y| ^ m
          ≤ 13 * (1 + Real.log (2 * T)) * (2 * T) ^ m := by gcongr
        _ ≤ 13 * (1 + Real.log (2 * T)) * (2 * u ^ (8 * m)) := by gcongr
    have hcpow : ‖(u : ℂ) ^ ((x : ℂ) + y * I)‖ = u ^ x := by
      rw [Complex.norm_cpow_eq_rpow_re_of_pos hu0]; simp
    have hden : T ≤ ‖(x : ℂ) + y * I‖ := by
      rw [← hy]; simpa using abs_im_le_norm ((x : ℂ) + y * I)
    have hexp : x + 8 * m ≤ E := by
      rcases le_total (1 - (s.re + x)) 0 with h | h
      · rw [hm, max_eq_right h]; rw [hE]; linarith [le_max_right (1 - s.re + 7 * η) 2]
      · rw [hm, max_eq_left h]; rw [hE]
        rw [ha] at hxa
        linarith [le_max_left (1 - s.re + 7 * η) 2]
    have hux : u ^ x * u ^ (8 * m) ≤ u ^ E := by
      rw [← Real.rpow_add hu0]
      exact Real.rpow_le_rpow_of_exponent_le hu hexp
    rw [Fz, norm_div, norm_mul, hcpow]
    rw [div_le_iff₀ (lt_of_lt_of_le hT0 hden)]
    have hKT : K * T = 26 * (1 + Real.log (2 * T)) * u ^ E := by
      rw [hK]; field_simp
    calc ‖riemannZeta (s + ((x : ℂ) + y * I))‖ * u ^ x
        ≤ 13 * (1 + Real.log (2 * T)) * (2 * u ^ (8 * m)) * u ^ x := by gcongr
      _ = 26 * (1 + Real.log (2 * T)) * (u ^ x * u ^ (8 * m)) := by ring
      _ ≤ 26 * (1 + Real.log (2 * T)) * u ^ E := by gcongr
      _ = K * T := hKT.symm
      _ ≤ K * ‖(x : ℂ) + y * I‖ := by
          have : 0 ≤ K := by rw [hK]; positivity
          exact mul_le_mul_of_nonneg_left hden this
  refine (intervalIntegral.norm_integral_le_of_norm_le_const hpt).trans ?_
  rw [abs_of_pos (by linarith : (0 : ℝ) < 2 - a)]
  have hK0 : 0 ≤ K := by rw [hK]; positivity
  calc K * (2 - a) ≤ K * 2 := by gcongr; linarith
    _ = 52 * (1 + Real.log (2 * T)) * u ^ E / T := by rw [hK]; ring

end

end Fayed.InputsProof.PerronAux
