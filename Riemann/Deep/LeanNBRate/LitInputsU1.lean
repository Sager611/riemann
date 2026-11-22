import Riemann.Deep.LeanNBRate.LitInputs

/-!
# Literature inputs used by the upper bound (U1 and the closure)
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- **Plancherel's theorem on L¹ ∩ L²** (W. Rudin, *Real and Complex Analysis*, 3rd ed., McGraw-Hill
1987, Theorem 9.13, with the Fourier transform of 9.1: "$\hat f(t) = \int_{-\infty}^{\infty} f(x)
e^{-ixt}\,dm(x)$", "$dm(x) = (2\pi)^{-1/2}\,dx$"): "One can associate to each $f \in L^2$ a function
$\hat f \in L^2$ so that the following properties hold: (a) If $f \in L^1 \cap L^2$, then $\hat f$ is
the previously defined Fourier transform of $f$. (b) For every $f \in L^2$, $\|\hat f\|_2 = \|f\|_2$."
Lean form: for f ∈ L¹ ∩ L² (Lebesgue measure), (a) and (b) give t ↦ ∫ f(x)e^{−ixt}dx ∈ L² and
∫|∫ f(x)e^{−ixt}dx|²dt = 2π ∫|f|² (Rudin's normalisation dm = dx/√(2π) multiplied out). -/
def Rudin1987_Thm9_13 : Prop :=
  ∀ f : ℝ → ℂ, Integrable f → MemLp f 2 →
    MemLp (fun t : ℝ => ∫ x : ℝ, f x * Complex.exp (-(x * t * I))) 2 ∧
    ∫ t : ℝ, ‖∫ x : ℝ, f x * Complex.exp (-(x * t * I))‖ ^ 2 = 2 * Real.pi * ∫ x : ℝ, ‖f x‖ ^ 2

end Riemann.Deep.LeanNBRate

