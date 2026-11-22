import Riemann.Deep.LeanNBRate.Bernstein
import Riemann.Deep.LeanNBRate.FourierDeriv

/-!
# [N] Lemma Λ: frequency decay of a windowed group component

For `r` with nodes in `[−w, w]` and the time window `ψ = twin ℓ`,
`F(ζ) = ∫ r(u) ψ(u) e^{−iζu} du = Λ_r(Θ_ζ)`, `Θ_ζ(θ) = ∫ ψ(u) e^{−i(ζ−θ)u} du`; the multiplier lemma
with `G = Θ_ζ · cutoff(−w, w, η)` gives decay `(K/(ℓ d))^N` in `d = dist(ζ, [−w, w]) ≥ 2η`.
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The time window as a complex function. -/
def psiC (ℓ : ℝ) (u : ℝ) : ℂ := (twin ℓ u : ℂ)

theorem twin_hasCompactSupport {ℓ : ℝ} (hℓ : 0 < ℓ) : HasCompactSupport (twin ℓ) :=
  cutoff_hasCompactSupport (by positivity)

theorem psiC_isTest {ℓ : ℝ} (hℓ : 0 < ℓ) : IsTest (psiC ℓ) :=
  ⟨ofRealCLM.contDiff.comp (twin_gevB hℓ).smooth,
    (twin_hasCompactSupport hℓ).comp_left ofReal_zero⟩

/-- The modulated window `u ↦ ψ(u) e^{−iζu}`. -/
def modPsi (ℓ ζ : ℝ) (u : ℝ) : ℂ := psiC ℓ u * cexp (-(I * ζ * u))

theorem modPsi_isTest {ℓ : ℝ} (hℓ : 0 < ℓ) (ζ : ℝ) : IsTest (modPsi ℓ ζ) := by
  refine ⟨(psiC_isTest hℓ).smooth.mul ?_, (psiC_isTest hℓ).supp.mul_right⟩
  exact (contDiff_const.mul ofRealCLM.contDiff).neg.cexp

/-- `Λ_r(G)` at `x = 0` is `Σ a_{jk} G^{(k)}(δ_j)`. -/
theorem ExpSum.Lam_zero (r : ExpSum) (G : ℝ → ℂ) :
    r.Lam 0 G = ∑ j, ∑ k ∈ Finset.range (r.mult j), r.coef j k * iteratedDeriv k G (r.node j) := by
  unfold ExpSum.Lam
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
  rw [Finset.sum_range_succ]
  have hz : ∀ l ∈ Finset.range k, (Nat.choose k l : ℂ) * (I * ((0 : ℝ) : ℂ)) ^ (k - l) *
      iteratedDeriv l G (r.node j) = 0 := by
    intro l hl
    rw [Finset.mem_range] at hl
    rw [ofReal_zero, mul_zero, zero_pow (by omega), mul_zero, zero_mul]
  rw [Finset.sum_eq_zero hz]
  simp

/-- `F(ζ) = ∫ r ψ e^{−iζu} = Λ_r(Θ_ζ)`, `Θ_ζ = fint (modPsi ℓ ζ)`. -/
theorem ExpSum.integral_mul_modPsi (r : ExpSum) {ℓ : ℝ} (hℓ : 0 < ℓ) (ζ : ℝ) :
    ∫ u : ℝ, r.eval u * modPsi ℓ ζ u = r.Lam 0 (fint (modPsi ℓ ζ)) := by
  have hT := modPsi_isTest hℓ ζ
  rw [r.Lam_zero]
  simp_rw [hT.iteratedDeriv_fint]
  unfold fint pmul ExpSum.eval
  have hint : ∀ (j : Fin r.J) (k : ℕ), Integrable (fun u : ℝ =>
      r.coef j k * ((I * u) ^ k * modPsi ℓ ζ u * cexp (I * r.node j * u))) := by
    intro j k
    refine Integrable.const_mul ?_ _
    have hb := (hT.pmul_test k).integrable
    refine (hb.mul_bdd (g := fun u : ℝ => cexp (I * r.node j * u)) (c := 1) ?_ ?_).congr
      (Eventually.of_forall fun u => ?_)
    · exact (by fun_prop : Continuous fun u : ℝ => cexp (I * r.node j * u)).aestronglyMeasurable
    · refine Eventually.of_forall fun u => ?_
      rw [show (I * r.node j * u : ℂ) = ((r.node j * u : ℝ) : ℂ) * I by push_cast; ring,
        Complex.norm_exp_ofReal_mul_I]
    · simp only [pmul]
  have hexp : ∀ u : ℝ, (∑ j, ∑ k ∈ Finset.range (r.mult j),
      r.coef j k * (I * u) ^ k * cexp (I * r.node j * u)) * modPsi ℓ ζ u =
      ∑ j, ∑ k ∈ Finset.range (r.mult j),
        r.coef j k * ((I * u) ^ k * modPsi ℓ ζ u * cexp (I * r.node j * u)) := by
    intro u
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun k _ => ?_
    ring
  simp_rw [hexp]
  rw [integral_finset_sum _ fun j _ => integrable_finset_sum _ fun k _ => hint j k]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [integral_finset_sum _ fun k _ => hint j k]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [integral_const_mul]

/-- `Θ_ζ^{(a)}(θ) = 2π · kern ((iu)^a ψ) (ζ − θ)`. -/
theorem iteratedDeriv_fint_modPsi {ℓ : ℝ} (hℓ : 0 < ℓ) (ζ : ℝ) (a : ℕ) (θ : ℝ) :
    iteratedDeriv a (fint (modPsi ℓ ζ)) θ =
      (2 * Real.pi : ℂ) * kern (pmul a (psiC ℓ)) (ζ - θ) := by
  rw [(modPsi_isTest hℓ ζ).iteratedDeriv_fint, kern_eq]
  have hπ : (2 * Real.pi : ℂ) ≠ 0 := by exact_mod_cast (show (0:ℝ) < 2 * Real.pi by positivity).ne'
  rw [← mul_assoc, mul_inv_cancel₀ hπ, one_mul]
  unfold fint pmul modPsi
  refine integral_congr_ae (Eventually.of_forall fun u => ?_)
  simp only
  have : cexp (-(I * ζ * u)) * cexp (I * θ * u) = cexp (-(I * ((ζ - θ : ℝ) : ℂ) * u)) := by
    rw [← Complex.exp_add]; push_cast; ring_nf
  calc (I * u) ^ a * (psiC ℓ u * cexp (-(I * ζ * u))) * cexp (I * θ * u)
      = (I * u) ^ a * psiC ℓ u * (cexp (-(I * ζ * u)) * cexp (I * θ * u)) := by ring
    _ = cexp (-(I * ((ζ - θ : ℝ) : ℂ) * u)) * ((I * u) ^ a * psiC ℓ u) := by rw [this]; ring
    _ = cexp (-(I * u * ((ζ - θ : ℝ) : ℂ))) * ((I * u) ^ a * psiC ℓ u) := by ring_nf

end Riemann.Deep.LeanNBRate
