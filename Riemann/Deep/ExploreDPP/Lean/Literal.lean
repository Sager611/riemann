import Riemann.Deep.ExploreDPP.Lean.WeakDuality

/-!
# Lemma 1 with the brief's literal hypothesis list

The brief's Lemma 1 also assumes k0 ≥ 1 and concludes J_θ(f) ≤ f^(0). weak_duality drops
k0 ≥ 1 (not needed) and writes f^(0) = 1; this wrapper restores both.
-/

open Real MeasureTheory FourierTransform
open scoped Real

namespace ExploreDPP

theorem weak_duality_literal {h lam θ : ℝ} {k0 : ℕ} {η : ℕ → ℝ}
    (hh : 0 < h) (_hk0 : 1 ≤ k0) (hlam : (k0 : ℝ) * h = lam) (hη : Summable η)
    (hcore : ∀ k : ℕ, 1 ≤ k → k < k0 → η k = -1) (hge : ∀ k : ℕ, k0 ≤ k → -1 ≤ η k)
    (hθ0 : 0 ≤ θ) (hθ : θ ≤ 1 / (2 * h))
    (hP : ∀ α, 0 ≤ α → α ≤ 1 / (2 * h) → (if α ≤ θ then α else 0) ≤ P h η α)
    {f : ℝ → ℝ} (hf : Admissible lam f) : J θ f ≤ (𝓕 (fun x => (f x : ℂ)) 0).re := by
  rw [hf.fourier_zero, Complex.one_re]
  exact weak_duality hh hlam hη hcore hge hθ0 hθ hP hf

end ExploreDPP

