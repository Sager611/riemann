/- UNCOMPILED DRAFT (nblplean, 2026-11-12). Not checked by Lean: every compile of it exceeds the 290 MB
   per-process cap (imports alone: Pow.Real 334 MB, LpSeminorm.Defs 343 MB, LSeries.RiemannZeta 342 MB).
   See ../../LEAN.md. Do not import from checked files. -/
import Riemann.Deep.NBLp.Lean.draft.Sec7

/-!
# NBLp in Lean, part 5: §7 A4 (B(ε) ≤ K_p ε)

From A3 on the line σ' = 1/p − ε/2 (HM σ' Φ_ε Ψ_ε), Lemma H, Lemma Γ(b), the toolkit's
Plancherel identity and A4(b) (the exact value of J, JEval), PROOF.md l. 282–310.
-/

open Complex MeasureTheory Set Filter Topology
open scoped ENNReal

namespace NBLp

/-! ## integrals on the real line -/

lemma sq_norm_line (σ t : ℝ) : ‖(σ : ℂ) + t * I‖ ^ 2 = σ ^ 2 + t ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_add_mul_I]

lemma inv_sq_add_eq {σ : ℝ} (hσ : 0 < σ) (t : ℝ) :
    (σ ^ 2 + t ^ 2)⁻¹ = (σ ^ 2)⁻¹ * (1 + (t / σ) ^ 2)⁻¹ := by
  have h1 : σ ^ 2 + t ^ 2 ≠ 0 := by positivity
  have h2 : σ ≠ 0 := hσ.ne'
  field_simp

lemma integrable_inv_sq_add {σ : ℝ} (hσ : 0 < σ) :
    Integrable (fun t : ℝ => (σ ^ 2 + t ^ 2)⁻¹) := by
  have h := (integrable_inv_one_add_sq.comp_div hσ.ne').const_mul (σ ^ 2)⁻¹
  exact h.congr (Eventually.of_forall fun t => (inv_sq_add_eq hσ t).symm)

lemma integral_inv_sq_add {σ : ℝ} (hσ : 0 < σ) :
    ∫ t : ℝ, (σ ^ 2 + t ^ 2)⁻¹ = Real.pi / σ := by
  simp_rw [inv_sq_add_eq hσ]
  have h3 := Measure.integral_comp_div (fun u : ℝ => (1 + u ^ 2)⁻¹) σ
  simp only at h3
  rw [integral_const_mul, h3, integral_univ_inv_one_add_sq, abs_of_pos hσ, smul_eq_mul]
  field_simp

/-- (3 + |t|)^a (1/4 + t²)^{−1} is integrable for 0 ≤ a < 1. -/
lemma integrable_poly_div {a : ℝ} (ha0 : 0 ≤ a) (ha : a < 1) :
    Integrable (fun t : ℝ => (3 + |t|) ^ a * (1 / 4 + t ^ 2)⁻¹) := by
  have hmaj : Integrable (fun t : ℝ => 16 * (1 + ‖t‖ ^ 2) ^ (-(2 - a) / 2)) :=
    (integrable_rpow_neg_one_add_norm_sq (by simp; linarith)).const_mul 16
  refine hmaj.mono' ?_ (Eventually.of_forall fun t => ?_)
  · refine Continuous.aestronglyMeasurable ?_
    exact (Continuous.rpow_const (by fun_prop) (fun _ => Or.inr ha0)).mul
      (Continuous.inv₀ (by fun_prop) (fun t => by positivity))
  · have hw : 1 ≤ 1 + t ^ 2 := by nlinarith [sq_nonneg t]
    have hsq : |t| ≤ Real.sqrt (1 + t ^ 2) := by
      rw [← Real.sqrt_sq_eq_abs]; exact Real.sqrt_le_sqrt (by linarith)
    have hs1 : 1 ≤ Real.sqrt (1 + t ^ 2) := by
      rw [show (1:ℝ) = Real.sqrt 1 by simp]; exact Real.sqrt_le_sqrt hw
    have h3 : 3 + |t| ≤ 4 * Real.sqrt (1 + t ^ 2) := by linarith
    have hpos : 0 < 3 + |t| := by positivity
    have h4 : (3 + |t|) ^ a ≤ 4 * (1 + t ^ 2) ^ (a / 2) := by
      calc (3 + |t|) ^ a ≤ (4 * Real.sqrt (1 + t ^ 2)) ^ a := Real.rpow_le_rpow hpos.le h3 ha0
        _ = 4 ^ a * (1 + t ^ 2) ^ (a / 2) := by
            rw [Real.mul_rpow (by norm_num) (Real.sqrt_nonneg _), Real.sqrt_eq_rpow,
              ← Real.rpow_mul (by linarith)]
            ring_nf
        _ ≤ 4 * (1 + t ^ 2) ^ (a / 2) := by
            have : (4:ℝ) ^ a ≤ 4 := by
              calc (4:ℝ) ^ a ≤ 4 ^ (1:ℝ) := Real.rpow_le_rpow_of_exponent_le (by norm_num) ha.le
                _ = 4 := Real.rpow_one 4
            exact mul_le_mul_of_nonneg_right this (Real.rpow_nonneg (by linarith) _)
    have h5 : (1 / 4 + t ^ 2)⁻¹ ≤ ((1 + t ^ 2) / 4)⁻¹ :=
      inv_anti₀ (by positivity) (by nlinarith [sq_nonneg t])
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity), Real.norm_eq_abs, sq_abs]
    calc (3 + |t|) ^ a * (1 / 4 + t ^ 2)⁻¹ ≤ (4 * (1 + t ^ 2) ^ (a / 2)) * ((1 + t ^ 2) / 4)⁻¹ :=
          mul_le_mul h4 h5 (by positivity) (by positivity)
      _ = 16 * (1 + t ^ 2) ^ (-(2 - a) / 2) := by
          have e : -(2 - a) / 2 = a / 2 - 1 := by ring
          rw [e, Real.rpow_sub (by linarith), Real.rpow_one, inv_div]
          field_simp
          ring

lemma gCG_le (t : ℝ) : gCG t ≤ 12 * ((3 + |t|) ^ (3 / 4 : ℝ) * (1 / 4 + t ^ 2)⁻¹) := by
  have hy : 1 ≤ 2 + |t| := by linarith [abs_nonneg t]
  have hl : Real.log (2 + |t|) ≤ 4 * (2 + |t|) ^ (1 / 4 : ℝ) := by
    have := Real.log_le_rpow_div (by linarith : (0:ℝ) ≤ 2 + |t|) (by norm_num : (0:ℝ) < 1 / 4)
    linarith [this, show (2 + |t|) ^ (1 / 4 : ℝ) / (1 / 4) = 4 * (2 + |t|) ^ (1 / 4 : ℝ) by ring]
  have h1 : 1 ≤ (2 + |t|) ^ (1 / 4 : ℝ) := Real.one_le_rpow hy (by norm_num)
  have h8 : 8 + Real.log (2 + |t|) ≤ 12 * (2 + |t|) ^ (1 / 4 : ℝ) := by linarith
  have hm : (2 + |t|) ^ (1 / 4 : ℝ) * (2 + |t|) ^ (1 / 2 : ℝ) ≤ (3 + |t|) ^ (3 / 4 : ℝ) := by
    rw [← Real.rpow_add (by linarith)]
    norm_num
    exact Real.rpow_le_rpow (by linarith) (by linarith) (by norm_num)
  have hi : 0 ≤ (1 / 4 + t ^ 2)⁻¹ := inv_nonneg.mpr (by positivity)
  have hh : 0 ≤ (2 + |t|) ^ (1 / 2 : ℝ) := Real.rpow_nonneg (by linarith) _
  unfold gCG
  calc (8 + Real.log (2 + |t|)) * (2 + |t|) ^ (1 / 2 : ℝ) * (1 / 4 + t ^ 2)⁻¹
      ≤ (12 * (2 + |t|) ^ (1 / 4 : ℝ)) * (2 + |t|) ^ (1 / 2 : ℝ) * (1 / 4 + t ^ 2)⁻¹ := by
        gcongr
    _ = 12 * ((2 + |t|) ^ (1 / 4 : ℝ) * (2 + |t|) ^ (1 / 2 : ℝ)) * (1 / 4 + t ^ 2)⁻¹ := by ring
    _ ≤ 12 * (3 + |t|) ^ (3 / 4 : ℝ) * (1 / 4 + t ^ 2)⁻¹ := by gcongr
    _ = 12 * ((3 + |t|) ^ (3 / 4 : ℝ) * (1 / 4 + t ^ 2)⁻¹) := by ring

lemma continuous_gCG : Continuous gCG := by
  unfold gCG
  refine ((continuous_const.add (Continuous.log (by fun_prop) (fun t => by positivity))).mul
    (Continuous.rpow_const (by fun_prop) (fun _ => Or.inr (by norm_num)))).mul
    (Continuous.inv₀ (by fun_prop) (fun t => by positivity))

lemma integrable_gCG : Integrable gCG := by
  refine ((integrable_poly_div (a := 3 / 4) (by norm_num) (by norm_num)).const_mul 12).mono'
    continuous_gCG.aestronglyMeasurable (Eventually.of_forall fun t => ?_)
  rw [Real.norm_eq_abs, abs_of_nonneg (gCG_nonneg t)]
  exact gCG_le t

/-- e^{8ε}y^ε ≤ 1 + ε e⁴ (8 + log y) y^{1/2} for 0 ≤ ε ≤ 1/2, y ≥ 1 (A4(c)). -/
lemma exp_bound {ε y : ℝ} (hε : 0 ≤ ε) (hε2 : ε ≤ 1 / 2) (hy : 1 ≤ y) :
    Real.exp (8 * ε) * y ^ ε ≤ 1 + ε * (Real.exp 4 * ((8 + Real.log y) * y ^ (1 / 2 : ℝ))) := by
  have hly := Real.log_nonneg hy
  set u := 8 * ε + ε * Real.log y with hu
  have e : Real.exp (8 * ε) * y ^ ε = Real.exp u := by
    rw [hu, Real.exp_add, Real.rpow_def_of_pos (by linarith), mul_comm (Real.log y) ε]
  have hpos := Real.exp_pos u
  have h1 : Real.exp u ≤ 1 + u * Real.exp u := by
    have h0 := Real.add_one_le_exp (-u)
    rw [Real.exp_neg] at h0
    have : (1 - u) * Real.exp u ≤ 1 := by
      calc (1 - u) * Real.exp u = (-u + 1) * Real.exp u := by ring
        _ ≤ (Real.exp u)⁻¹ * Real.exp u := mul_le_mul_of_nonneg_right h0 hpos.le
        _ = 1 := inv_mul_cancel₀ hpos.ne'
    linarith
  have h3 : u * Real.exp u ≤ ε * (Real.exp 4 * ((8 + Real.log y) * y ^ (1 / 2 : ℝ))) := by
    rw [← e]
    have hl : 0 ≤ 8 + Real.log y := by linarith
    have hexp : Real.exp (8 * ε) ≤ Real.exp 4 := Real.exp_le_exp.mpr (by linarith)
    have hyε : y ^ ε ≤ y ^ (1 / 2 : ℝ) := Real.rpow_le_rpow_of_exponent_le hy hε2
    have hu' : u = ε * (8 + Real.log y) := by rw [hu]; ring
    rw [hu']
    calc ε * (8 + Real.log y) * (Real.exp (8 * ε) * y ^ ε)
        ≤ ε * (8 + Real.log y) * (Real.exp 4 * y ^ (1 / 2 : ℝ)) :=
          mul_le_mul_of_nonneg_left (mul_le_mul hexp hyε (Real.rpow_nonneg (by linarith) _)
            (Real.exp_pos 4).le) (mul_nonneg hε hl)
      _ = ε * (Real.exp 4 * ((8 + Real.log y) * y ^ (1 / 2 : ℝ))) := by ring
  rw [e]; linarith

/-! ## real values at real points -/

lemma RE_real {ε a : ℝ} (hε : 0 < ε) (ha : 1 < a) :
    RE ε (a : ℂ) = (((∑' n : ℕ, 1 / (n : ℝ) ^ a) / ∑' n : ℕ, 1 / (n : ℝ) ^ (a + ε) : ℝ) : ℂ) := by
  have h1 : (a : ℂ) + ε ≠ 1 := by
    intro h0; have := congrArg Complex.re h0; simp at this; linarith
  have e : (a : ℂ) + ε = ((a + ε : ℝ) : ℂ) := by push_cast; ring
  rw [RE, invZeta, if_neg h1, e, zeta_real ha, zeta_real (by linarith), Complex.ofReal_div,
    div_eq_mul_inv]

lemma RE_ge_one {ε a : ℝ} (hε : 0 < ε) (ha : 1 < a) : 1 ≤ (RE ε (a : ℂ)).re := by
  rw [RE_real hε ha, Complex.ofReal_re]
  have hpos : 0 < ∑' n : ℕ, 1 / (n : ℝ) ^ (a + ε) :=
    lt_of_lt_of_le (one_div_pos.mpr (by linarith)) (tz_ge (by linarith))
  rw [one_le_div hpos]
  exact tz_anti ha (by linarith)

lemma inv_zeta_real {ε : ℝ} (hε : 0 < ε) :
    (riemannZeta (1 + ε))⁻¹ = (((∑' n : ℕ, 1 / (n : ℝ) ^ (1 + ε))⁻¹ : ℝ) : ℂ) := by
  have e : (1 : ℂ) + ε = ((1 + ε : ℝ) : ℂ) := by push_cast; ring
  rw [e, zeta_real (by linarith), Complex.ofReal_inv]

lemma inv_tz_le {ε : ℝ} (hε : 0 < ε) : (∑' n : ℕ, 1 / (n : ℝ) ^ (1 + ε))⁻¹ ≤ ε := by
  have h := tz_ge (a := 1 + ε) (by linarith)
  rw [show 1 + ε - 1 = ε by ring] at h
  have hpos : 0 < ∑' n : ℕ, 1 / (n : ℝ) ^ (1 + ε) := lt_of_lt_of_le (one_div_pos.mpr hε) h
  rw [inv_le_comm₀ hpos hε]; simpa using h

/-! ## A4 -/

/-- §7 A4(b): J := ∫ R_ε(σ'+it) |s|^{−2} dt converges absolutely and equals
2π[R(2σ')/(2σ') − 1/((2σ'−1)ζ(1+ε))] (PROOF.md l. 289–303). -/
def JEval (p ε : ℝ) : Prop :=
  Integrable (fun t : ℝ => RE ε (((1 / p - ε / 2 : ℝ) : ℂ) + t * I) /
      ((‖((1 / p - ε / 2 : ℝ) : ℂ) + t * I‖ ^ 2 : ℝ) : ℂ)) ∧
  ∫ t : ℝ, RE ε (((1 / p - ε / 2 : ℝ) : ℂ) + t * I) /
      ((‖((1 / p - ε / 2 : ℝ) : ℂ) + t * I‖ ^ 2 : ℝ) : ℂ) =
    2 * Real.pi * (RE ε ((2 * (1 / p - ε / 2) : ℝ) : ℂ) / ((2 * (1 / p - ε / 2) : ℝ) : ℂ) -
      1 / (((2 * (1 / p - ε / 2) - 1 : ℝ) : ℂ) * riemannZeta (1 + ε)))

/-- §7 A4: B(ε) = ‖Φ_ε‖²_{L²_σ'} ≤ K_p ε. -/
theorem A4 {HM : ℝ → (ℝ → ℂ) → (ℝ → ℂ) → Prop} (hTK : MellinToolkit HM) (hH : LemmaH)
    (hGb : LemmaGb) {p ε : ℝ} (h : EpsOK p ε) (hi : ZeroFree (1 / p))
    (hA3 : HM (1 / p - ε / 2) (PhiE ε) (fun t => PsiE ε (((1 / p - ε / 2 : ℝ) : ℂ) + t * I)))
    (hJ : JEval p ε) :
    MemL2w (1 / p - ε / 2) (PhiE ε) ∧ L2w (1 / p - ε / 2) (PhiE ε) ≤ Kp p * ε := by
  refine ⟨hTK.mem hA3, ?_⟩
  have hσ2 : 1 / 2 ≤ 1 / p - ε / 2 := h.s1_ge
  have hσ0 : 0 < 1 / p - ε / 2 := by linarith
  have hσe : 1 / p - ε / 2 + ε < 1 := by linarith [h.s1e_lt]
  obtain ⟨hPint, hPl⟩ := hTK.plancherel hA3
  obtain ⟨hJint, hJv⟩ := hJ
  set σ : ℝ := 1 / p - ε / 2 with hσ
  have hsre : ∀ t : ℝ, ((σ : ℂ) + t * I).re = σ := by intro t; simp
  have hsim : ∀ t : ℝ, ((σ : ℂ) + t * I).im = t := by intro t; simp
  have hzeros : ∀ ρ : ℂ, riemannZeta ρ = 0 → 0 < ρ.re → ρ.re ≤ 1 / p :=
    fun ρ hρ _ => not_lt.mp fun h' => hi ρ h' hρ
  -- Lemma H and Lemma Γ(b) on the line
  have hR : ∀ t : ℝ, ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 ≤ Real.exp (8 * ε) * (2 + |t|) ^ ε := by
    intro t
    have hne1 : (σ : ℂ) + t * I ≠ 1 := by
      intro h1; have := congrArg Complex.re h1; rw [hsre] at this; simp at this
      linarith [h.ip1]
    have hse1 : (σ : ℂ) + t * I + ε ≠ 1 := by
      intro h1; have := congrArg Complex.re h1; simp at this; linarith
    have hHt := (hH (1 / p) ε ((σ : ℂ) + t * I) hzeros h.pos (by rw [hsre, hσ]) hne1
      (by rw [hsre]; exact hσ0)).2 hse1
    have hz0 : riemannZeta ((σ : ℂ) + t * I + ε) ≠ 0 := hi _ (by simp; linarith [h.pos])
    have hRle : ‖RE ε ((σ : ℂ) + t * I)‖ ≤ ‖Gr ε ((σ : ℂ) + t * I)‖ := by
      rw [RE, invZeta, if_neg hse1, norm_mul, norm_inv,
        mul_inv_le_iff₀ (norm_pos_iff.mpr hz0)]
      exact hHt.trans (le_of_eq (mul_comm _ _))
    have hG := (hGb ε ((σ : ℂ) + t * I) h.pos (by linarith [h.half]) (by rw [hsre]; exact hσ2)
      (fun u hu => Or.inl (by rw [hsre]; linarith [hu.2]))).2 (by rw [hsre]; exact hσe)
    rw [hsim] at hG
    calc ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 ≤ ‖Gr ε ((σ : ℂ) + t * I)‖ ^ 2 := by gcongr
      _ ≤ _ := hG
  -- pointwise expansion |Ψ|² = (1 − 2 Re R + |R|²)/|s|²
  have hpt : ∀ t : ℝ, ‖PsiE ε ((σ : ℂ) + t * I)‖ ^ 2 = (σ ^ 2 + t ^ 2)⁻¹ -
      2 * ((RE ε ((σ : ℂ) + t * I)).re * (σ ^ 2 + t ^ 2)⁻¹) +
      ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 * (σ ^ 2 + t ^ 2)⁻¹ := by
    intro t
    have hsq : ‖1 - RE ε ((σ : ℂ) + t * I)‖ ^ 2 =
        1 - 2 * (RE ε ((σ : ℂ) + t * I)).re + ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 := by
      rw [Complex.sq_norm, Complex.sq_norm, Complex.normSq_apply, Complex.normSq_apply]
      simp only [sub_re, one_re, sub_im, one_im]
      ring
    rw [PsiE, norm_div, div_pow, sq_norm_line, hsq]
    ring
  have I1 := integrable_inv_sq_add hσ0
  have I2 : Integrable (fun t : ℝ => (RE ε ((σ : ℂ) + t * I)).re * (σ ^ 2 + t ^ 2)⁻¹) := by
    refine hJint.re.congr (Eventually.of_forall fun t => ?_)
    simp only [RCLike.re_to_complex]
    rw [Complex.div_ofReal_re, sq_norm_line, div_eq_mul_inv]
  have I3 : Integrable (fun t : ℝ => ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 * (σ ^ 2 + t ^ 2)⁻¹) := by
    refine ((hPint.sub I1).add (I2.const_mul 2)).congr (Eventually.of_forall fun t => ?_)
    simp only [Pi.add_apply, Pi.sub_apply]
    rw [hpt t]; ring
  have I4 : Integrable (fun t : ℝ => (σ ^ 2 + t ^ 2)⁻¹ + ε * (Real.exp 4 * gCG t)) :=
    I1.add ((integrable_gCG.const_mul (Real.exp 4)).const_mul ε)
  have hb3 : ∀ t : ℝ, ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 * (σ ^ 2 + t ^ 2)⁻¹ ≤
      (σ ^ 2 + t ^ 2)⁻¹ + ε * (Real.exp 4 * gCG t) := by
    intro t
    have hy : 1 ≤ 2 + |t| := by linarith [abs_nonneg t]
    have hE := exp_bound h.pos.le h.half hy
    have hσt : (σ ^ 2 + t ^ 2)⁻¹ ≤ (1 / 4 + t ^ 2)⁻¹ :=
      inv_anti₀ (by positivity) (by nlinarith [hσ2])
    have hX : 0 ≤ Real.exp 4 * ((8 + Real.log (2 + |t|)) * (2 + |t|) ^ (1 / 2 : ℝ)) := by
      have := Real.log_nonneg hy
      have := Real.rpow_nonneg (by linarith : (0:ℝ) ≤ 2 + |t|) (1 / 2 : ℝ)
      positivity
    have hinv0 : 0 ≤ (σ ^ 2 + t ^ 2)⁻¹ := inv_nonneg.mpr (by positivity)
    calc ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 * (σ ^ 2 + t ^ 2)⁻¹
        ≤ (1 + ε * (Real.exp 4 * ((8 + Real.log (2 + |t|)) * (2 + |t|) ^ (1 / 2 : ℝ)))) *
            (σ ^ 2 + t ^ 2)⁻¹ := mul_le_mul_of_nonneg_right ((hR t).trans hE) hinv0
      _ = (σ ^ 2 + t ^ 2)⁻¹ + ε * (Real.exp 4 * ((8 + Real.log (2 + |t|)) *
            (2 + |t|) ^ (1 / 2 : ℝ))) * (σ ^ 2 + t ^ 2)⁻¹ := by ring
      _ ≤ (σ ^ 2 + t ^ 2)⁻¹ + ε * (Real.exp 4 * ((8 + Real.log (2 + |t|)) *
            (2 + |t|) ^ (1 / 2 : ℝ))) * (1 / 4 + t ^ 2)⁻¹ := by
          gcongr
      _ = (σ ^ 2 + t ^ 2)⁻¹ + ε * (Real.exp 4 * gCG t) := by unfold gCG; ring
  have hint3 : ∫ t : ℝ, ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 * (σ ^ 2 + t ^ 2)⁻¹ ≤
      Real.pi / σ + CG * ε := by
    calc ∫ t : ℝ, ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 * (σ ^ 2 + t ^ 2)⁻¹
        ≤ ∫ t : ℝ, ((σ ^ 2 + t ^ 2)⁻¹ + ε * (Real.exp 4 * gCG t)) := integral_mono I3 I4 hb3
      _ = Real.pi / σ + CG * ε := by
          rw [integral_add I1 ((integrable_gCG.const_mul (Real.exp 4)).const_mul ε),
            integral_inv_sq_add hσ0, integral_const_mul, integral_const_mul, CG]
          ring
  -- the value of ∫ Re R/|s|² (A4(b))
  have h2σ : 1 < 2 * σ := by have := h.two_s1; have := h.two_s1_pos; linarith
  have hReJ : ∫ t : ℝ, (RE ε ((σ : ℂ) + t * I)).re * (σ ^ 2 + t ^ 2)⁻¹ =
      2 * Real.pi * ((RE ε ((2 * σ : ℝ) : ℂ)).re / (2 * σ) -
        (∑' n : ℕ, 1 / (n : ℝ) ^ (1 + ε))⁻¹ / (2 * σ - 1)) := by
    have e1 : ∫ t : ℝ, (RE ε ((σ : ℂ) + t * I)).re * (σ ^ 2 + t ^ 2)⁻¹ =
        (∫ t : ℝ, RE ε ((σ : ℂ) + t * I) / ((‖(σ : ℂ) + t * I‖ ^ 2 : ℝ) : ℂ)).re := by
      rw [← RCLike.re_to_complex, ← integral_re hJint]
      congr 1; ext t
      rw [RCLike.re_to_complex, Complex.div_ofReal_re, sq_norm_line, div_eq_mul_inv]
    rw [e1, hJv, RE_real h.pos h2σ, ← Complex.ofReal_ofNat, mul_one_div, one_div_mul_eq_div]
    have hden : ((2 * σ - 1 : ℝ) : ℂ) * riemannZeta (1 + ε) =
        (((2 * σ - 1) * ∑' n : ℕ, 1 / (n : ℝ) ^ (1 + ε) : ℝ) : ℂ) := by
      have e : (1 : ℂ) + ε = ((1 + ε : ℝ) : ℂ) := by push_cast; ring
      rw [e, zeta_real (by linarith [h.pos]), Complex.ofReal_mul]
    rw [hden]
    simp only [← Complex.ofReal_div, ← Complex.ofReal_sub, ← Complex.ofReal_mul,
      ← Complex.ofReal_ofNat, Complex.ofReal_re, ← Complex.ofReal_one]
    rw [RE_real h.pos h2σ, Complex.ofReal_re]
    field_simp
  -- assemble
  have hRr : 1 ≤ (RE ε ((2 * σ : ℝ) : ℂ)).re := RE_ge_one h.pos h2σ
  have hZ : (∑' n : ℕ, 1 / (n : ℝ) ^ (1 + ε))⁻¹ / (2 * σ - 1) ≤ ε * (2 * p / (2 - p)) := by
    have hd : (2 - p) / (2 * p) ≤ 2 * σ - 1 := h.two_s1
    have hd0 := h.two_s1_pos
    rw [div_le_iff₀ (by linarith)]
    have h1 := inv_tz_le h.pos
    have h1' : 0 ≤ (∑' n : ℕ, 1 / (n : ℝ) ^ (1 + ε))⁻¹ := inv_nonneg.mpr
      (tsum_nonneg fun n => by positivity)
    calc (∑' n : ℕ, 1 / (n : ℝ) ^ (1 + ε))⁻¹ ≤ ε := h1
      _ = ε * (2 * p / (2 - p)) * ((2 - p) / (2 * p)) := by
          have : p ≠ 0 := h.p0.ne'
          have : 2 - p ≠ 0 := by linarith [h.p2]
          field_simp
      _ ≤ ε * (2 * p / (2 - p)) * (2 * σ - 1) := by
          apply mul_le_mul_of_nonneg_left hd
          exact mul_nonneg h.pos.le (div_nonneg (by linarith [h.p0]) (by linarith [h.p2]))
  have htot : ∫ t : ℝ, ‖PsiE ε ((σ : ℂ) + t * I)‖ ^ 2 =
      ∫ t : ℝ, (σ ^ 2 + t ^ 2)⁻¹ - 2 * ∫ t : ℝ, (RE ε ((σ : ℂ) + t * I)).re * (σ ^ 2 + t ^ 2)⁻¹ +
        ∫ t : ℝ, ‖RE ε ((σ : ℂ) + t * I)‖ ^ 2 * (σ ^ 2 + t ^ 2)⁻¹ := by
    simp_rw [hpt]
    rw [integral_add (I1.sub (I2.const_mul 2)) I3, integral_sub I1 (I2.const_mul 2),
      integral_const_mul]
  have hbig : ∫ t : ℝ, ‖PsiE ε ((σ : ℂ) + t * I)‖ ^ 2 ≤ 2 * Real.pi * (Kp p * ε) := by
    rw [htot, integral_inv_sq_add hσ0, hReJ]
    have hπ := Real.pi_pos
    have hK : 2 * Real.pi * (Kp p * ε) = 2 * Real.pi * (ε * (4 * p / (2 - p))) + CG * ε := by
      unfold Kp; field_simp; ring
    rw [hK]
    have hσinv : 0 < Real.pi / σ := div_pos hπ hσ0
    have e2 : 2 * (2 * Real.pi * ((RE ε ((2 * σ : ℝ) : ℂ)).re / (2 * σ))) =
        2 * (Real.pi / σ) * (RE ε ((2 * σ : ℝ) : ℂ)).re := by
      field_simp
    have e3 : 2 * Real.pi * (ε * (4 * p / (2 - p))) = 4 * Real.pi * (ε * (2 * p / (2 - p))) := by
      ring
    nlinarith [hint3, hRr, hZ, hσinv, e2, e3]
  have hPl' : L2w σ (PhiE ε) = (2 * Real.pi)⁻¹ * ∫ t : ℝ, ‖PsiE ε ((σ : ℂ) + t * I)‖ ^ 2 := hPl
  rw [hPl']
  have h2π : 0 < 2 * Real.pi := by positivity
  rw [inv_mul_le_iff₀ h2π]
  exact hbig

end NBLp

