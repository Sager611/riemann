/-
Riemann/Deep/Yoshida/Windows.lean — Weil positivity on every window `L ≤ 1/3`, sorry-free.

Route (`Hybrid.lean`): below `log 2`, `Re Q(f) = (1/2π) ∫ |F(r)|² (2π μ(r) + J(r)) dr`, where the
pole term became the weight `J(r) = ∫_{−2a}^{2a} 2cosh(x/2) cos(rx) dx` (`PoleKernel.lean`).
Layer cake with `|F|² ≤ 2aN`, `∫|F|² = 2πN`; the weight is bounded below cell by cell
(digamma series, Taylor bounds for cos/sin, Lipschitz bound for `J`) by the generated
certificate `Third.lean` (139 cells, `.scratch/yoshida2/gen.py`).
-/
import Riemann.Deep.Yoshida.Third

namespace Fayed.Yoshida

open Riemann.Deep.Weil

/-- **Weil positivity on the window `1/3`.** -/
theorem weilPos_third : WeilPos (1 / 3) := by
  have h := Third.weilPos
  norm_num at h
  exact h

/-- **Yoshida-type window**: Weil positivity for every window `L ≤ 1/3`. -/
theorem weilPos_le_third (L : ℝ) (hL : L ≤ 1 / 3) : WeilPos L :=
  WeilPos.mono hL weilPos_third

end Fayed.Yoshida
