import Mathlib

/-!
# Round 50 — fresh directions after 49 rounds (8 lenses)

| lens | direction | odds |
|------|-----------|------|
| analytic number theory | self-adjoint operator with det = ξ (Hilbert–Pólya) | very low |
| operator theory | unitary adelic scattering operator, Cayley transform | very low |
| mathematical physics | Hilbert–Pólya from the explicit-formula kernel | very low |
| logic | cut-elimination for an explicit-formula theory forbidding off-line witnesses | very low |
| additive combinatorics | off-line zero ⇒ excess additive energy of prime phases; inverse theorem | very low |
| probability | random prime-power signs; L² drift contradiction | very low |
| arithmetic geometry | Arakelov H¹ with Deligne-style purity | very low |
| complexity | least counterexample: incompressibility vs analytic compression | very low |

Three of eight lenses converge on Hilbert–Pólya; its finite core (self-adjoint ⇒ real spectrum) is below,
together with the least-counterexample principle used by the complexity lens.
-/

namespace Fayed.R50

/-- Hilbert–Pólya core: eigenvalues of a symmetric operator are real. -/
theorem symmetric_eigenvalue_real {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℂ E]
    {T : E →ₗ[ℂ] E} (hT : T.IsSymmetric) {μ : ℂ} (hμ : Module.End.HasEigenvalue T μ) :
    μ.im = 0 :=
  Complex.conj_eq_iff_im.mp (hT.conj_eigenvalue_eq_self hμ)

/-- Least counterexample principle. -/
theorem exists_least {P : ℕ → Prop} [DecidablePred P] (h : ∃ n, P n) :
    ∃ n, P n ∧ ∀ m, m < n → ¬ P m :=
  ⟨Nat.find h, Nat.find_spec h, fun _ hm => Nat.find_min h hm⟩

end Fayed.R50
