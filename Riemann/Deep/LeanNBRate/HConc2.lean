import Riemann.Deep.LeanNBRate.HConc1
import Riemann.Deep.LeanNBRate.U2aux
import Riemann.Deep.LeanNBRate.Contour8
import Riemann.Deep.LeanNBRate.Window

/-!
# The concrete zero configuration up to height T and the removability of F = H/ζ at the low zeros
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem zeta_order_ne_top {ρ : ℂ} (hρ : ρ ≠ 1) : analyticOrderAt riemannZeta ρ ≠ ⊤ := by
  have h2 : riemannZeta 2 ≠ 0 := riemannZeta_ne_zero_of_one_lt_re (by norm_num)
  have h2' : analyticOrderAt riemannZeta 2 ≠ ⊤ := by
    rw [analyticOrderAt_eq_zero.mpr (Or.inr h2)]; exact ENat.zero_ne_top
  exact analyticOn_riemannZeta.analyticOrderAt_ne_top_of_isPreconnected
    (isConnected_compl_singleton_of_one_lt_rank (by simp) 1).isPreconnected (x := 2)
    (by simp) (by simpa using hρ) h2'

theorem mult_pos_of_NTZ {ρ : ℂ} (hρ : ρ ∈ NTZ) : 1 ≤ mult ρ := by
  have hρ1 := NTZ_ne_one hρ
  have hfin := zeta_order_ne_top hρ1
  by_contra h
  push Not at h
  have h0 : mult ρ = 0 := by omega
  have hz : analyticOrderAt riemannZeta ρ = 0 := by
    rw [← Nat.cast_analyticOrderNatAt hfin]
    unfold mult at h0
    rw [h0]; rfl
  rcases analyticOrderAt_eq_zero.mp hz with h1 | h1
  · exact h1 (analyticOn_riemannZeta ρ hρ1)
  · exact h1 hρ.1

/-- The zero configuration: ordinates γ with 1/2 + iγ a nontrivial zero, |γ| ≤ T. -/
def cfg (T : ℝ) : ExpSum where
  J := (lowOrd T).card
  node := fun j => (((lowOrd T).equivFin.symm j : lowOrd T) : ℝ)
  mult := fun j => mult ((1 / 2 : ℂ) + ((((lowOrd T).equivFin.symm j : lowOrd T) : ℝ) : ℂ) * I)
  coef := fun _ _ => 0

theorem cfg_zeroConfig (T : ℝ) : ZeroConfig (cfg T) T := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro i j h
    exact (lowOrd T).equivFin.symm.injective (Subtype.ext h)
  · intro j
    exact (mem_lowOrd.mp ((lowOrd T).equivFin.symm j).2).1
  · intro j
    exact (mem_lowOrd.mp ((lowOrd T).equivFin.symm j).2).2
  · intro j
    exact mult_pos_of_NTZ (mem_lowOrd.mp ((lowOrd T).equivFin.symm j).2).1
  · intro j
    exact le_rfl

theorem cfg_node_of_zero (hRH : RiemannHypothesis) {T : ℝ} {ρ : ℂ} (hρ : ρ ∈ NTZ) (hT : |ρ.im| ≤ T) :
    ∃ j₀ : Fin (cfg T).J, (cfg T).node j₀ = -ρ.im ∧ mult ρ ≤ (cfg T).mult j₀ := by
  have hre := NTZ_re_half hRH hρ
  have hconj : (1 / 2 : ℂ) + ((-ρ.im : ℝ) : ℂ) * I = conj ρ := by
    apply Complex.ext <;> simp [hre]
  have hmem : -ρ.im ∈ lowOrd T := by
    rw [mem_lowOrd, hconj]
    exact ⟨conj_mem_NTZ hρ, by rw [abs_neg]; exact hT⟩
  refine ⟨(lowOrd T).equivFin ⟨-ρ.im, hmem⟩, ?_, ?_⟩
  · show ((((lowOrd T).equivFin.symm ((lowOrd T).equivFin ⟨-ρ.im, hmem⟩)) : lowOrd T) : ℝ) = -ρ.im
    rw [Equiv.symm_apply_apply]
  · show mult ρ ≤ mult ((1 / 2 : ℂ) +
      ((((lowOrd T).equivFin.symm ((lowOrd T).equivFin ⟨-ρ.im, hmem⟩)) : lowOrd T) : ℝ) * I)
    rw [Equiv.symm_apply_apply]
    show mult ρ ≤ mult ((1 / 2 : ℂ) + ((-ρ.im : ℝ) : ℂ) * I)
    rw [hconj, mult_conj ρ (NTZ_ne_one hρ)]

/-- The concrete H = (1 - E_low) Φ. -/
def Hc (ℓ a : ℝ) (n : ℕ) (r : ExpSum) (s : ℂ) : ℂ :=
  (1 - Elow ℓ r.eval (r.dmod 0).eval s) * Phi a n s

/-- **HZeros for the concrete H**: F = H/ζ is removable at every zero below T₂. -/
theorem hZeros_Hc (hRH : RiemannHypothesis) {T₂ ℓ a : ℝ} {n : ℕ} (hT₂ : 0 < T₂) (hℓ : 0 ≤ ℓ)
    (hgood : ∀ σ ∈ Icc (-1 : ℝ) 2, riemannZeta (σ + T₂ * I) ≠ 0) {c : Fin (cfg T₂).J → ℕ → ℂ}
    (hc : ∀ d : Fin (cfg T₂).J → ℕ → ℂ, NForm ℓ ((cfg T₂).withCoef d).eval (((cfg T₂).withCoef d).dmod 0).eval
      ((cfg T₂).withCoef c).eval (((cfg T₂).withCoef c).dmod 0).eval = ((cfg T₂).withCoef d).eval 0) :
    HZeros (Hc ℓ a n ((cfg T₂).withCoef c)) T₂ := by
  refine ⟨hT₂, hgood, ?_⟩
  intro ρ hρ0 hre0 hre1 him
  have hρ1 : ρ ≠ 1 := fun h => by rw [h] at hre1; simp at hre1
  have hnt : ¬∃ k : ℕ, ρ = -2 * (k + 1) := by
    rintro ⟨k, hk⟩
    have h := congrArg Complex.re hk
    simp at h
    have : (0 : ℝ) ≤ k := k.cast_nonneg
    linarith
  have hhalf : ρ.re = 1 / 2 := hRH ρ hρ0 hnt hρ1
  have hNTZ : ρ ∈ NTZ := ⟨hρ0, by rw [hhalf]; norm_num, by rw [hhalf]; norm_num⟩
  obtain ⟨j₀, hnode, hmult⟩ := cfg_node_of_zero hRH hNTZ (T := T₂) him.le
  have hζa : AnalyticAt ℂ riemannZeta ρ := analyticOn_riemannZeta ρ hρ1
  obtain ⟨g₁, hg₁a, hg₁0, hζeq⟩ := hζa.analyticOrderAt_ne_top.mp (zeta_order_ne_top hρ1)
  obtain ⟨r, hr⟩ : ∃ r : ExpSum, r = (cfg T₂).withCoef c := ⟨_, rfl⟩
  have hEd : Differentiable ℂ (Elow ℓ r.eval (r.dmod 0).eval) :=
    differentiable_Elow hℓ r.eval (r.dmod 0).eval ((r.dmod 0).dmod 0).eval
      (fun x => (r.dmod 0).hasDerivAt_eval x) r.continuous_eval (r.dmod 0).continuous_eval
      ((r.dmod 0).dmod 0).continuous_eval
  have hfa : AnalyticAt ℂ (fun s => 1 - Elow ℓ r.eval (r.dmod 0).eval s) ρ :=
    ((differentiable_const _).sub hEd).analyticAt ρ
  have hρf : ρ - 1 / 2 = -(((cfg T₂).node j₀ : ℝ) : ℂ) * I := by
    rw [hnode]
    apply Complex.ext <;> simp [hhalf]
  have hord : ((mult ρ : ℕ) : ℕ∞) ≤ analyticOrderAt (fun s => 1 - Elow ℓ r.eval (r.dmod 0).eval s) ρ := by
    rw [natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero hfa]
    intro i hi
    rw [hr]
    exact iteratedDeriv_one_sub_Elow hℓ hc (j₀ := j₀) (by omega) hρf
  obtain ⟨g₂, hg₂a, hfeq⟩ := (natCast_le_analyticOrderAt hfa).mp hord
  refine ⟨fun z => g₂ z * Phi a n z / g₁ z, ?_, ?_⟩
  · exact (hg₂a.continuousAt.mul (differentiable_Phi a n).continuous.continuousAt).div hg₁a.continuousAt hg₁0
  · have h1 : ∀ᶠ z in 𝓝[≠] ρ, z ≠ ρ := self_mem_nhdsWithin
    filter_upwards [nhdsWithin_le_nhds hζeq, nhdsWithin_le_nhds hfeq, h1] with z hz1 hz2 hz3
    simp only [smul_eq_mul] at hz1 hz2
    have hne : (z - ρ) ^ (mult ρ) ≠ 0 := pow_ne_zero _ (sub_ne_zero.mpr hz3)
    rw [← hr]
    unfold Fq Hc
    rw [hz1, hz2]
    have e : (z - ρ) ^ mult ρ * g₂ z * Phi a n z / ((z - ρ) ^ analyticOrderNatAt riemannZeta ρ * g₁ z) =
        (z - ρ) ^ mult ρ * (g₂ z * Phi a n z) / ((z - ρ) ^ mult ρ * g₁ z) := by
      unfold mult; ring
    rw [e, mul_div_mul_left _ _ hne]

end Riemann.Deep.LeanNBRate

