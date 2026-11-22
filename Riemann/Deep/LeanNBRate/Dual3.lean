import Riemann.Deep.LeanNBRate.Dual2
import Riemann.Deep.LeanNBRate.U2final

/-!
# Dual problem: the representer bound (from U2) and the Cauchy–Schwarz bound for E_low
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate
open scoped InnerProductSpace

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem oubil_self_eq (ℓ : ℝ) (f f' : ℝ → ℂ) :
    oubil ℓ f f' f f' = ((∫ u in (0 : ℝ)..ℓ, (‖f' u‖ ^ 2 + 1 / 4 * ‖f u‖ ^ 2) : ℝ) : ℂ) := by
  unfold oubil
  have e : ∀ u, f' u * conj (f' u) + (1 / 4 : ℂ) * (f u * conj (f u)) =
      ((‖f' u‖ ^ 2 + 1 / 4 * ‖f u‖ ^ 2 : ℝ) : ℂ) := by
    intro u
    rw [Complex.mul_conj, Complex.mul_conj, Complex.normSq_eq_norm_sq, Complex.normSq_eq_norm_sq]
    push_cast; ring
  simp_rw [e]
  rw [intervalIntegral.integral_ofReal]

/-- N(f, f) is real: ∫(|f′|²+¼|f|²) + ½(|f(0)|²+|f(ℓ)|²). -/
def Nself (ℓ : ℝ) (f f' : ℝ → ℂ) : ℝ :=
  (∫ u in (0 : ℝ)..ℓ, (‖f' u‖ ^ 2 + 1 / 4 * ‖f u‖ ^ 2)) + 1 / 2 * (‖f 0‖ ^ 2 + ‖f ℓ‖ ^ 2)

theorem NForm_self (ℓ : ℝ) (f f' : ℝ → ℂ) : NForm ℓ f f' f f' = (Nself ℓ f f' : ℂ) := by
  unfold NForm Nself
  rw [oubil_self_eq, Complex.mul_conj, Complex.mul_conj, Complex.normSq_eq_norm_sq,
    Complex.normSq_eq_norm_sq]
  push_cast; ring

theorem Nself_nonneg {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (f f' : ℝ → ℂ) : 0 ≤ Nself ℓ f f' := by
  unfold Nself
  have : 0 ≤ ∫ u in (0 : ℝ)..ℓ, (‖f' u‖ ^ 2 + 1 / 4 * ‖f u‖ ^ 2) :=
    intervalIntegral.integral_nonneg hℓ fun u _ => by positivity
  positivity

theorem oubil_re_le_Nself (ℓ : ℝ) (f f' : ℝ → ℂ) : (oubil ℓ f f' f f').re ≤ Nself ℓ f f' := by
  rw [oubil_self_eq, Complex.ofReal_re]
  unfold Nself
  have : 0 ≤ 1 / 2 * (‖f 0‖ ^ 2 + ‖f ℓ‖ ^ 2) := by positivity
  linarith

theorem ZeroConfig.withCoef {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) (c : Fin P.J → ℕ → ℂ) :
    ZeroConfig (P.withCoef c) T :=
  ⟨hP.inj, hP.zero, hP.height, hP.mult_pos, hP.mult_le⟩

/-- **The representer bound** from U2: N(p*) ≤ (S + ε)/ℓ. -/
theorem representer_bound {P : ExpSum} {T ℓ ε : ℝ} (hℓ : 0 < ℓ) (hP : ZeroConfig P T)
    (hU2 : ∀ (Q : ExpSum), ZeroConfig Q T →
      ‖Q.eval 0‖ ^ 2 ≤ (zeroSum + ε) / ℓ * (oubil ℓ Q.eval (Q.dmod 0).eval Q.eval (Q.dmod 0).eval).re)
    (hε : 0 ≤ zeroSum + ε) {c : Fin P.J → ℕ → ℂ}
    (hc : ∀ d : Fin P.J → ℕ → ℂ, NForm ℓ (P.withCoef d).eval ((P.withCoef d).dmod 0).eval
      (P.withCoef c).eval ((P.withCoef c).dmod 0).eval = (P.withCoef d).eval 0) :
    Nself ℓ (P.withCoef c).eval ((P.withCoef c).dmod 0).eval ≤ (zeroSum + ε) / ℓ := by
  set p := P.withCoef c
  set N := Nself ℓ p.eval (p.dmod 0).eval
  have hN0 : 0 ≤ N := Nself_nonneg hℓ.le _ _
  have h1 : p.eval 0 = (N : ℂ) := by rw [← hc c, NForm_self]
  have h2 := hU2 p (hP.withCoef c)
  rw [h1, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hN0] at h2
  have h3 : N ^ 2 ≤ (zeroSum + ε) / ℓ * N :=
    h2.trans (mul_le_mul_of_nonneg_left (oubil_re_le_Nself _ _ _) (by positivity))
  rcases eq_or_lt_of_le hN0 with h | h
  · rw [← h]; positivity
  · have : N * N ≤ (zeroSum + ε) / ℓ * N := by nlinarith
    exact le_of_mul_le_mul_right this h

/-- e_s(u) = e^{-(s-1/2)u} and its derivative. -/
def eS (s : ℂ) (u : ℝ) : ℂ := cexp (-(s - 1 / 2) * u)
def eS' (s : ℂ) (u : ℝ) : ℂ := -(s - 1 / 2) * cexp (-(s - 1 / 2) * u)

theorem continuous_eS (s : ℂ) : Continuous (eS s) := by unfold eS; fun_prop
theorem continuous_eS' (s : ℂ) : Continuous (eS' s) := by unfold eS'; fun_prop

/-- E_low(s) = N(e_s, p). -/
def Elow (ℓ : ℝ) (p p' : ℝ → ℂ) (s : ℂ) : ℂ := NForm ℓ (eS s) (eS' s) p p'

theorem norm_NForm_le {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (F F' p p' : ℝ → ℂ) (hF : Continuous F) (hF' : Continuous F')
    (hp : Continuous p) (hp' : Continuous p') :
    ‖NForm ℓ F F' p p'‖ ≤ Real.sqrt (Nself ℓ F F') * Real.sqrt (Nself ℓ p p') := by
  rw [← inner_Temb hℓ F F' p p' hF hF' hp hp']
  refine (norm_inner_le_norm _ _).trans ?_
  have h1 := norm_sq_Temb hℓ F F' hF hF'
  have h2 := norm_sq_Temb hℓ p p' hp hp'
  rw [NForm_self, Complex.ofReal_re] at h1 h2
  rw [← h1, ← h2, Real.sqrt_sq (norm_nonneg _), Real.sqrt_sq (norm_nonneg _), mul_comm]

theorem norm_Elow_le {ℓ : ℝ} (hℓ : 0 ≤ ℓ) (p p' : ℝ → ℂ) (hp : Continuous p) (hp' : Continuous p') (s : ℂ) :
    ‖Elow ℓ p p' s‖ ≤ Real.sqrt (Nself ℓ (eS s) (eS' s)) * Real.sqrt (Nself ℓ p p') :=
  norm_NForm_le hℓ _ _ _ _ (continuous_eS s) (continuous_eS' s) hp hp'

end Riemann.Deep.LeanNBRate

