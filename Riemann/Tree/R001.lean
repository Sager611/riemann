import Riemann.Fayed.Theorem1

/-!
# Round 1 — frontier R1: repair step (55) of Fayed

| id | direction | status |
|----|-----------|--------|
| R1.1 | gap-aware L'Hôpital / Stolz along b_M | DEAD: valid versions see F/G → 0 |
| R1.2 | σ₀-sensitivity audit of (34),(35) | DEAD: profile realizable ∀σ₀ (fayed_profile_realizable) |
| R1.3 | other comparisons (s̄, 1-s, χ, ∂σ, many p) | DEAD: restate zero symmetry; same profile |
| R1.4 | Nyman–Beurling via B̄ₖ | CIRCULAR: density ≡ RH |
| R1.5 | function-field transfer | DEAD: no Frobenius / Riemann–Roch positivity analog |
| R1.6 | Weil positivity / Li λₙ | CIRCULAR: ≡ RH |
| R1.7 | Hilbert–Pólya operator from B̄ₖ flow | CIRCULAR: open |
| R1.8 | de Bruijn–Newman heat flow | CIRCULAR: Λ ≤ 0 ≡ RH |
| R1.9 | Bohr / period means in u | CIRCULAR: B̄ₖ mean zero (bernoulli_mean_zero); positivity ≡ RH |
| R1.10 | automated identity search over {s,1-s,s̄,1-s̄} | CIRCULAR: every orbit value vanishes; restates FE |

Meta-lesson: Fayed's local data (Euler–Maclaurin + Perron + reflection) hold for every σ₀;
any contradiction needs a GLOBAL σ-sensitive inequality.
Filter for all later rounds: Davenport–Heilbronn (FE, real coefficients, off-line zeros) —
an argument blind to the Euler product would also "prove" RH for it, so it is dead.
-/

open Filter Topology Set Real

namespace Fayed.R1

/-- The paper's own profile at a hypothetical zero σ₀: F ≈ (1-2σ₀)·(first Fourier mode of B̄₂)/u. -/
noncomputable def Fσ (σ₀ : ℝ) (u : ℝ) : ℝ := (1 - 2 * σ₀) * Fc u

lemma deriv_Fσ (σ₀ : ℝ) {u : ℝ} (hu : 0 < u) : deriv (Fσ σ₀) u = (1 - 2 * σ₀) * deriv Fc u :=
  deriv_const_mul _ (hasDerivAt_Fc hu).differentiableAt

/-- **R1.2 no-go.** For every σ₀, a zeta-free pair (F,G) has exactly the limit profile the paper
derives for its own (F,G) at an off-line zero: (39),(42),(43),(50) with the paper's signs, (54).
So no argument using only these facts can force σ₀ = 1/2. -/
theorem fayed_profile_realizable (σ₀ : ℝ) {γ : ℝ} (hγ : 0 < γ) (hγ' : γ < 1 / 4) (M₀ : ℕ) :
    (∀ u, 0 < u → DifferentiableAt ℝ (Fσ σ₀) u) ∧
    Tendsto (Fσ σ₀) (along (Dγ γ M₀)) (𝓝 0) ∧
    Tendsto Gc (along (Dγ γ M₀)) (𝓝 0) ∧
    Tendsto (fun u => Fσ σ₀ u / Gc u) (along (Dγ γ M₀)) (𝓝 0) ∧
    (σ₀ < 1 / 2 → Tendsto (fun u => deriv (Fσ σ₀) u / deriv Gc u) (along (Dγ γ M₀)) atBot) ∧
    (1 / 2 < σ₀ → Tendsto (fun u => deriv (Fσ σ₀) u / deriv Gc u) (along (Dγ γ M₀)) atTop) ∧
    (∀ M, ∀ u ∈ Ioo (aγ γ M₀ M) (bγ γ M₀ M),
      |Fσ σ₀ u / Gc u - Fσ σ₀ (bγ γ M₀ M) / Gc (bγ γ M₀ M)| ≤
        |1 - 2 * σ₀| * (2 * (√(aγ γ M₀ M))⁻¹)) := by
  obtain ⟨hF, hG, -, -, -, -, hF0, hG0, hR, hosc, -, hFG, -⟩ := counterexample hγ hγ' M₀
  have hpos : ∀ᶠ u in along (Dγ γ M₀), 0 < u :=
    (eventually_gt_atTop 0).filter_mono (along_le_atTop _)
  have hratio : ∀ᶠ u in along (Dγ γ M₀),
      (1 - 2 * σ₀) * (deriv Fc u / deriv Gc u) = deriv (Fσ σ₀) u / deriv Gc u := by
    filter_upwards [hpos] with u hu
    rw [deriv_Fσ σ₀ hu, mul_div_assoc]
  refine ⟨fun u hu => (hF u hu).const_mul _, ?_, hG0, ?_, ?_, ?_, ?_⟩
  · have := hF0.const_mul (1 - 2 * σ₀)
    rw [mul_zero] at this
    exact this
  · have := hFG.const_mul (1 - 2 * σ₀)
    rw [mul_zero] at this
    refine this.congr (fun u => ?_)
    simp only [Fσ, mul_div_assoc]
  · intro h
    exact ((hR.const_mul_atBot (by linarith))).congr' hratio
  · intro h
    exact ((hR.const_mul_atBot_of_neg (by linarith))).congr' hratio
  · intro M u hu
    have := hosc M u hu
    calc |Fσ σ₀ u / Gc u - Fσ σ₀ (bγ γ M₀ M) / Gc (bγ γ M₀ M)|
        = |1 - 2 * σ₀| * |Fc u / Gc u - Fc (bγ γ M₀ M) / Gc (bγ γ M₀ M)| := by
          rw [← abs_mul]; congr 1; simp only [Fσ]; ring
      _ ≤ |1 - 2 * σ₀| * (2 * (√(aγ γ M₀ M))⁻¹) :=
          mul_le_mul_of_nonneg_left this (abs_nonneg _)

/-- **R1.9 no-go.** Every B̄ₖ (k ≥ 1) has mean zero over a period, so period averages in u
annihilate Fayed's (1-2σ₀)·B̄ₖ signal. -/
theorem bernoulli_mean_zero {k : ℕ} (hk : k ≠ 0) : ∫ x in (0 : ℝ)..1, bernoulliFun k x = 0 :=
  integral_bernoulliFun_eq_zero hk

end Fayed.R1
