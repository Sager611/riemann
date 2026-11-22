import Riemann.Deep.LeanNBRate.Kernel

/-!
# The multiplier lemma (replaces [N] Lemma A)

For an exponential sum `r` (nodes `δ_j`, jets) and a test function `G`,
`∫ kern G v · r(x+v) dv = Σ_{j,k} a_{jk} e^{iδ_j x} Σ_{l ≤ k} C(k,l) (ix)^{k−l} G^{(l)}(δ_j)`,
and, by Turán–Nazarov outside `[0,ℓ]`, this is at most
`c^{M−1}·sup_{[0,ℓ]}|r| · ∫ |kern G v| (1 + |v|/ℓ)^{M−1} dv`.
-/

open MeasureTheory Filter Topology Set Complex FourierTransform

noncomputable section

namespace Riemann.Deep.LeanNBRate

namespace IsTest

variable {G : ℝ → ℂ} (hG : IsTest G)
include hG

theorem integrable_kern : Integrable (kern G) := by
  unfold kern
  have hπ : (2 * Real.pi : ℝ) ≠ 0 := by positivity
  exact (hG.integrable_fourier.comp_div hπ).const_mul _

theorem integrable_pow_mul_kern (n : ℕ) : Integrable (fun v : ℝ => (I * v) ^ n * kern G v) := by
  have := (hG.iterDeriv n).integrable_kern
  refine this.congr (Eventually.of_forall fun v => ?_)
  exact hG.kern_iteratedDeriv n v

theorem integrable_norm_pow_mul_kern (n : ℕ) : Integrable (fun v : ℝ => |v| ^ n * ‖kern G v‖) := by
  have := (hG.integrable_pow_mul_kern n).norm
  refine this.congr (Eventually.of_forall fun v => ?_)
  simp [norm_mul, norm_pow, Complex.norm_I, Complex.norm_real, Real.norm_eq_abs]

end IsTest

namespace ExpSum

variable (r : ExpSum)

/-- `Λ_r(e^{iθx} G)`, written through the jets of `G` at the nodes. -/
def Lam (x : ℝ) (G : ℝ → ℂ) : ℂ :=
  ∑ j, ∑ k ∈ Finset.range (r.mult j), r.coef j k * cexp (I * r.node j * x) *
    ∑ l ∈ Finset.range (k + 1), (Nat.choose k l : ℂ) * (I * x) ^ (k - l) *
      iteratedDeriv l G (r.node j)

theorem integral_kern_mul_eval {G : ℝ → ℂ} (hG : IsTest G) (x : ℝ) :
    ∫ v : ℝ, kern G v * r.eval (x + v) = r.Lam x G := by
  have hterm : ∀ (j : Fin r.J) (k l : ℕ), Integrable (fun v : ℝ =>
      r.coef j k * cexp (I * r.node j * x) * ((Nat.choose k l : ℂ) * (I * x) ^ (k - l) *
        (kern G v * (I * v) ^ l * cexp (I * r.node j * v)))) := by
    intro j k l
    refine Integrable.const_mul ?_ _
    refine Integrable.const_mul ?_ _
    have hb : Integrable (fun v : ℝ => (I * v) ^ l * kern G v) := hG.integrable_pow_mul_kern l
    refine (hb.mul_bdd (g := fun v : ℝ => cexp (I * r.node j * v)) (c := 1) ?_ ?_).congr
      (Eventually.of_forall fun v => ?_)
    · exact (by fun_prop : Continuous fun v : ℝ => cexp (I * r.node j * v)).aestronglyMeasurable
    · refine Eventually.of_forall fun v => ?_
      rw [show (I * r.node j * v : ℂ) = ((r.node j * v : ℝ) : ℂ) * I by push_cast; ring,
        Complex.norm_exp_ofReal_mul_I]
    · simp only; ring
  have hexp : ∀ v : ℝ, kern G v * r.eval (x + v) =
      ∑ j, ∑ k ∈ Finset.range (r.mult j), ∑ l ∈ Finset.range (k + 1),
        r.coef j k * cexp (I * r.node j * x) * ((Nat.choose k l : ℂ) * (I * x) ^ (k - l) *
          (kern G v * (I * v) ^ l * cexp (I * r.node j * v))) := by
    intro v
    unfold eval
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    have h1 : (I * ((x + v : ℝ) : ℂ)) ^ k =
        ∑ l ∈ Finset.range (k + 1), (I * v) ^ l * (I * x) ^ (k - l) * (Nat.choose k l : ℂ) := by
      rw [show (I * ((x + v : ℝ) : ℂ)) = I * v + I * x by push_cast; ring, add_pow]
    have h2 : cexp (I * r.node j * ((x + v : ℝ) : ℂ)) =
        cexp (I * r.node j * x) * cexp (I * r.node j * v) := by
      rw [← Complex.exp_add]; push_cast; ring_nf
    rw [h1, h2]
    simp only [Finset.mul_sum, Finset.sum_mul]
    refine Finset.sum_congr rfl fun l _ => ?_
    ring
  simp_rw [hexp]
  rw [integral_finset_sum _ fun j _ => integrable_finset_sum _ fun k _ =>
    integrable_finset_sum _ fun l _ => hterm j k l]
  unfold Lam
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [integral_finset_sum _ fun k _ => integrable_finset_sum _ fun l _ => hterm j k l]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [integral_finset_sum _ fun l _ => hterm j k l, Finset.mul_sum]
  refine Finset.sum_congr rfl fun l _ => ?_
  rw [integral_const_mul, integral_const_mul, hG.integral_kern_pow_mul_exp l (r.node j)]

end ExpSum

/-- Turán–Nazarov outside `[0, ℓ]`: `|r(x+v)| ≤ (c(1 + |v|/ℓ))^{M−1} K` for `x ∈ [0,ℓ]`. -/
theorem turan_shift (hT : Nazarov1993_Turan) : ∃ c : ℝ, 0 < c ∧
    ∀ (r : ExpSum) (ℓ K x v : ℝ), 0 < ℓ → 1 ≤ r.M → x ∈ Icc 0 ℓ →
      (∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ K) →
      ‖r.eval (x + v)‖ ≤ (c * (1 + |v| / ℓ)) ^ (r.M - 1) * K := by
  obtain ⟨c, hc, hT⟩ := turan_expSum hT
  refine ⟨c, hc, fun r ℓ K x v hℓ hM hx hK => ?_⟩
  have hK0 : 0 ≤ K := (norm_nonneg _).trans (hK 0 ⟨le_rfl, hℓ.le⟩)
  set y := x + v
  have hab : min 0 y < max ℓ y := lt_of_le_of_lt (min_le_left _ _) (lt_of_lt_of_le hℓ (le_max_left _ _))
  have hsub : Icc 0 ℓ ⊆ Icc (min 0 y) (max ℓ y) :=
    Icc_subset_Icc (min_le_left _ _) (le_max_left _ _)
  have hvol : (volume (Icc (0:ℝ) ℓ)).toReal = ℓ := by
    rw [Real.volume_Icc, ENNReal.toReal_ofReal (by linarith)]; ring
  have hpos : 0 < volume (Icc (0:ℝ) ℓ) := by
    rw [Real.volume_Icc]; exact ENNReal.ofReal_pos.mpr (by linarith)
  have hy : y ∈ Icc (min 0 y) (max ℓ y) := ⟨min_le_right _ _, le_max_right _ _⟩
  have h := hT r (min 0 y) (max ℓ y) (Icc 0 ℓ) K hM hab measurableSet_Icc hsub hpos hK y hy
  rw [hvol] at h
  refine h.trans (mul_le_mul_of_nonneg_right (pow_le_pow_left₀ ?_ ?_ _) hK0)
  · apply div_nonneg (by nlinarith [hab]) hℓ.le
  · rw [div_le_iff₀ hℓ, mul_assoc, add_mul, one_mul, div_mul_cancel₀ _ hℓ.ne']
    gcongr
    rcases le_total 0 y with h0 | h0 <;> rcases le_total ℓ y with h1 | h1
    · rw [min_eq_left h0, max_eq_right h1]; linarith [hx.1, hx.2, le_abs_self v]
    · rw [min_eq_left h0, max_eq_left h1]; linarith [abs_nonneg v]
    · linarith
    · rw [min_eq_right h0, max_eq_left h1]; linarith [hx.1, hx.2, neg_abs_le v]

end Riemann.Deep.LeanNBRate
