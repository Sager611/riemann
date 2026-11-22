import Riemann.Deep.LeanNBRate.LemmaE
import Riemann.Deep.LeanNBRate.Counting

/-!
# U2 assembly tools: low/high split of the OU form, final combination inequality
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem part_eval_union (P : ExpSum) (T₁ T₂ : Set ℝ) (hdisj : Disjoint T₁ T₂)
    (hcov : ∀ j, P.node j ∈ T₁ ∪ T₂) (u : ℝ) :
    P.eval u = (P.part T₁ 0).eval u + (P.part T₂ 0).eval u := by
  classical
  rw [P.part_eval, P.part_eval, ← Finset.sum_add_distrib]
  unfold ExpSum.eval
  refine Finset.sum_congr rfl fun j _ => ?_
  rcases hcov j with h | h
  · have h' : P.node j ∉ T₂ := fun h2 => Set.disjoint_left.mp hdisj h h2
    rw [if_pos h, if_neg h', add_zero, sub_zero]
  · have h' : P.node j ∉ T₁ := fun h1 => Set.disjoint_left.mp hdisj h1 h
    rw [if_neg h', if_pos h, zero_add, sub_zero]

theorem part_dmod_eval_union (P : ExpSum) (T₁ T₂ : Set ℝ) (hdisj : Disjoint T₁ T₂)
    (hcov : ∀ j, P.node j ∈ T₁ ∪ T₂) (u : ℝ) :
    (P.dmod 0).eval u = ((P.part T₁ 0).dmod 0).eval u + ((P.part T₂ 0).dmod 0).eval u := by
  have hF : P.eval = fun v => (P.part T₁ 0).eval v + (P.part T₂ 0).eval v := by
    funext v; exact part_eval_union P T₁ T₂ hdisj hcov v
  have h1 := ((P.part T₁ 0).hasDerivAt_eval u).add ((P.part T₂ 0).hasDerivAt_eval u)
  have h2 := P.hasDerivAt_eval u
  rw [hF] at h2
  exact h2.unique h1

theorem oubil_swap (ℓ : ℝ) (f f' g g' : ℝ → ℂ) :
    oubil ℓ g g' f f' = conj (oubil ℓ f f' g g') := by
  unfold oubil
  simp only [intervalIntegral, map_sub, ← integral_conj]
  congr 2 <;> (ext u; simp only [map_add, map_mul, Complex.conj_conj, map_div₀, map_one, map_ofNat]; ring)

theorem oubil_add_re (ℓ : ℝ) (f₁ f₁' f₂ f₂' : ℝ → ℂ) (h₁ : Continuous f₁) (h₁' : Continuous f₁')
    (h₂ : Continuous f₂) (h₂' : Continuous f₂') :
    (oubil ℓ (fun u => f₁ u + f₂ u) (fun u => f₁' u + f₂' u) (fun u => f₁ u + f₂ u)
      (fun u => f₁' u + f₂' u)).re =
      (oubil ℓ f₁ f₁' f₁ f₁').re + (oubil ℓ f₂ f₂' f₂ f₂').re + 2 * (oubil ℓ f₁ f₁' f₂ f₂').re := by
  have h := oubil_sum2 (Finset.univ : Finset (Fin 2)) (Finset.univ : Finset (Fin 2))
    (fun i => ![f₁, f₂] i) (fun i => ![f₁', f₂'] i) (fun i => ![f₁, f₂] i) (fun i => ![f₁', f₂'] i)
    (fun i => by fin_cases i <;> assumption) (fun i => by fin_cases i <;> assumption)
    (fun i => by fin_cases i <;> assumption) (fun i => by fin_cases i <;> assumption) ℓ
  simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons] at h
  rw [h, oubil_swap ℓ f₁ f₁' f₂ f₂']
  simp only [Complex.add_re, Complex.conj_re]
  ring

theorem oubil_self_re (ℓ : ℝ) (hℓ : 0 ≤ ℓ) (f f' : ℝ → ℂ) (hf : Continuous f) (hf' : Continuous f') :
    (oubil ℓ f f' f f').re = ∫ u in (0 : ℝ)..ℓ, (‖f' u‖ ^ 2 + 1 / 4 * ‖f u‖ ^ 2) := by
  unfold oubil
  have e : ∀ u, f' u * conj (f' u) + (1 / 4 : ℂ) * (f u * conj (f u)) =
      ((‖f' u‖ ^ 2 + 1 / 4 * ‖f u‖ ^ 2 : ℝ) : ℂ) := by
    intro u
    rw [Complex.mul_conj, Complex.mul_conj, Complex.normSq_eq_norm_sq, Complex.normSq_eq_norm_sq]
    push_cast; ring
  simp_rw [e]
  rw [intervalIntegral.integral_ofReal, Complex.ofReal_re]

/-- The windowed derivative mass is at most the OU form. -/
theorem psi_mass_le_oubil {ℓ : ℝ} (hℓ : 0 < ℓ) (f f' : ℝ → ℂ) (hf : Continuous f) (hf' : Continuous f') :
    ∫ u : ℝ, ‖f' u * psiC ℓ u‖ ^ 2 ≤ (oubil ℓ f f' f f').re := by
  rw [oubil_self_re ℓ hℓ.le f f' hf hf']
  have hpsi : ∀ u, ‖f' u * psiC ℓ u‖ ^ 2 = ‖f' u‖ ^ 2 * twin ℓ u ^ 2 := by
    intro u
    rw [norm_mul, mul_pow]
    unfold psiC
    rw [Complex.norm_real, Real.norm_eq_abs, sq_abs]
  have hcont : Continuous fun u => ‖f' u * psiC ℓ u‖ ^ 2 :=
    ((hf'.mul (psiC_isTest hℓ).continuous).norm.pow 2)
  have hzero : ∀ u ∉ Icc 0 ℓ, ‖f' u * psiC ℓ u‖ ^ 2 = 0 := by
    intro u hu
    rw [hpsi u]
    have : twin ℓ u = 0 := by
      apply twin_eq_zero hℓ
      simp only [mem_Icc, not_and_or, not_le] at hu
      rcases hu with h | h
      · left; exact h.le
      · right; exact h.le
    rw [this]; ring
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero hzero, integral_Icc_eq_integral_Ioc,
    ← intervalIntegral.integral_of_le hℓ.le]
  refine intervalIntegral.integral_mono_on hℓ.le (hcont.intervalIntegrable _ _)
    (((hf'.norm.pow 2).add (continuous_const.mul (hf.norm.pow 2))).intervalIntegrable _ _) ?_
  intro u _
  rw [hpsi u]
  have h1 : twin ℓ u ^ 2 ≤ 1 := by
    have := twin_le_one ℓ u
    have := twin_nonneg ℓ u
    nlinarith
  have h2 : 0 ≤ ‖f' u‖ ^ 2 := sq_nonneg _
  have h3 : 0 ≤ 1 / 4 * ‖f u‖ ^ 2 := by positivity
  nlinarith

/-- **The combination inequality** (pure algebra). -/
theorem combine_bound {Q Qlo Qhi B Dlo Y a b θ Slo Shi p0 plo phi ℓ : ℝ}
    (hℓ : 0 < ℓ) (hDlo : 0 ≤ Dlo) (hY : 0 ≤ Y) (hSlo : 0 ≤ Slo) (hShi : 0 ≤ Shi)
    (hab : a + b + θ < 1) (hθ : 0 ≤ θ) (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hQ : Q = Qlo + Qhi + 2 * B) (hlo : (1 - a) * Dlo ≤ Qlo) (hhi : (1 - b) * Y ≤ Qhi)
    (hB : |B| ≤ θ * Real.sqrt Dlo * Real.sqrt Y)
    (hp : p0 ≤ plo + phi) (hp0 : 0 ≤ p0) (hplo0 : 0 ≤ plo) (hphi0 : 0 ≤ phi)
    (hplo : plo ^ 2 ≤ Slo * (Dlo / ℓ)) (hphi : phi ^ 2 ≤ Shi / ℓ * Y) :
    p0 ^ 2 ≤ (Slo + Shi) / ℓ * Q / (1 - a - b - θ) := by
  have hden : 0 < 1 - a - b - θ := by linarith
  have hsq : 2 * Real.sqrt Dlo * Real.sqrt Y ≤ Dlo + Y := by
    have := sq_nonneg (Real.sqrt Dlo - Real.sqrt Y)
    have e1 := Real.sq_sqrt hDlo
    have e2 := Real.sq_sqrt hY
    nlinarith
  have hQ' : (1 - a - b - θ) * (Dlo + Y) ≤ Q := by
    have hB' : -(θ * (Dlo + Y)) ≤ 2 * B := by
      have := (abs_le.mp hB).1
      have h3 : θ * (2 * Real.sqrt Dlo * Real.sqrt Y) ≤ θ * (Dlo + Y) := mul_le_mul_of_nonneg_left hsq hθ
      nlinarith
    have hb0 : 0 ≤ b * Dlo := mul_nonneg hb hDlo
    have ha0 : 0 ≤ a * Y := mul_nonneg ha hY
    nlinarith
  have hp2 : p0 ^ 2 ≤ (Slo + Shi) / ℓ * (Dlo + Y) := by
    have hplo' : plo ≤ Real.sqrt (Slo / ℓ) * Real.sqrt Dlo := by
      rw [← Real.sqrt_mul (by positivity)]
      exact Real.le_sqrt_of_sq_le (by rw [div_mul_eq_mul_div, mul_div_assoc] at *; linarith)
    have hphi' : phi ≤ Real.sqrt (Shi / ℓ) * Real.sqrt Y := by
      rw [← Real.sqrt_mul (by positivity)]
      exact Real.le_sqrt_of_sq_le hphi
    have hcs := Finset.sum_mul_sq_le_sq_mul_sq (Finset.univ : Finset (Fin 2))
      ![Real.sqrt (Slo / ℓ), Real.sqrt (Shi / ℓ)] ![Real.sqrt Dlo, Real.sqrt Y]
    simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons] at hcs
    rw [Real.sq_sqrt (by positivity), Real.sq_sqrt (by positivity), Real.sq_sqrt hDlo,
      Real.sq_sqrt hY] at hcs
    have hsum : p0 ≤ Real.sqrt (Slo / ℓ) * Real.sqrt Dlo + Real.sqrt (Shi / ℓ) * Real.sqrt Y := by
      linarith
    have := pow_le_pow_left₀ hp0 hsum 2
    have e : (Slo + Shi) / ℓ = Slo / ℓ + Shi / ℓ := by ring
    rw [e]; linarith
  rw [le_div_iff₀ hden]
  calc p0 ^ 2 * (1 - a - b - θ) ≤ (Slo + Shi) / ℓ * (Dlo + Y) * (1 - a - b - θ) :=
        mul_le_mul_of_nonneg_right hp2 hden.le
    _ = (Slo + Shi) / ℓ * ((1 - a - b - θ) * (Dlo + Y)) := by ring
    _ ≤ (Slo + Shi) / ℓ * Q := mul_le_mul_of_nonneg_left hQ' (by positivity)

end Riemann.Deep.LeanNBRate

