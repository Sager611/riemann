module

public import Riemann.Deep.LeanLi.Thm2M
public import Riemann.Deep.LeanLi.Part13

/-!
# LeanLi: the main theorem — the Li-coefficient remainder laws (Claim of `discoveries/li-remainder.md`)
-/

@[expose] public section

open Real Finset Filter

namespace LeanLi

lemma Tf_id' {N : ℕ} (hN : 3 ≤ N) : (N : ℝ) * Tf N * Real.log N = N * √(N : ℝ) * √(Real.log N) := by
  have hl := log_nat_ge_one hN
  have hq : 0 < √(Real.log N) := Real.sqrt_pos.2 (by linarith)
  obtain ⟨q, hq0, hℓq⟩ : ∃ q : ℝ, 0 < q ∧ Real.log N = q ^ 2 :=
    ⟨√(Real.log N), hq, (Real.sq_sqrt (by linarith)).symm⟩
  have hsq : √(Real.log N) = q := by rw [hℓq, Real.sqrt_sq hq0.le]
  unfold Tf
  rw [Real.sqrt_div' _ (by linarith), hsq, hℓq]
  field_simp

/-- **Claim part 1, third clause** (Theorem 2 with Montgomery's pair-correlation theorem). -/
theorem part1_mont (L : LiteratureInputs) (hRH : RiemannHypothesis) :
    ∃ C : ℝ, ∀ᶠ N : ℕ in atTop, ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤
      C * ((N : ℝ) ^ (3 / 2 : ℝ) * Real.log N ^ (1 / 2 : ℝ) * Real.log (Real.log N)) := by
  obtain ⟨C, hC⟩ := thm2_mont L hRH
  refine ⟨C, ?_⟩
  filter_upwards [hC, eventually_ge_atTop 3] with N h hN
  have hN1 : (1 : ℝ) ≤ N := by
    have h3 : (3 : ℝ) ≤ N := by exact_mod_cast hN
    linarith
  rw [rpow_three_halves _ (by linarith), ← Real.sqrt_eq_rpow]
  have e := Tf_id' hN
  calc ∑ n ∈ Finset.Ioc N (2 * N), R n ^ 2 ≤ C * ((N : ℝ) * Tf N * Real.log N * Real.log (Real.log N)) := h
    _ = C * ((N : ℝ) * √(N : ℝ) * √(Real.log N) * Real.log (Real.log N)) := by rw [e]

/-- **Claim part 1.** -/
theorem part1 (L : LiteratureInputs) : Part1 :=
  ⟨c₁_num, fun hRH => ⟨part1_lower L hRH, part1_upper L hRH, part1_mont L hRH⟩⟩

/-- Claim 1's lower bound from RH and `LitCore` alone: no pair correlation (`MontgomeryPC` is not a field
of `LitCore`), no simplicity (zeros counted with multiplicity). -/
theorem part1_lower_rh_only (L : LitCore) (hRH : RiemannHypothesis) :
    ∀ ε > (0 : ℝ), ∀ᶠ N : ℕ in atTop,
      (c₁ - ε) * ((N : ℝ) ^ (3 / 2 : ℝ) * Real.log N ^ (1 / 2 : ℝ)) ≤ ∑ n ∈ Finset.Icc 1 N, R n ^ 2 :=
  part1_lower L hRH

/-- **Main theorem: the Li-coefficient remainder laws** (parts 0–5 of the Claim). -/
theorem li_remainder_laws (L : LiteratureInputs) : LiRemainderLaws :=
  ⟨part0 L, part1 L, part2 L, part3 L, part4 L, part5 L⟩

end LeanLi
