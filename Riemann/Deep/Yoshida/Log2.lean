/-
Riemann/Deep/Yoshida/Log2.lean — Weil positivity on every window `L < log 2` (Yoshida's range), sorry-free.

Route: `Affine.lean` (affine cells, chord bounds, wider support) + `AReflect.lean` (rational checkers) +
the generated certificate `WLog2.lean` at `a₁ = 1733/5000 > log 2 / 2`. For `f` supported in `[−a, a]`,
`2a < log 2`, the prime sum vanishes and the pole kernel may be cut at `2a₁`, so one certificate at `a₁`
covers every smaller window.
-/
import Riemann.Deep.Yoshida.WLog2

namespace Fayed.Yoshida

open Riemann.Deep.Weil

/-- **Yoshida's window**: Weil positivity for every `L < log 2`. -/
theorem weilPos_lt_log2 (L : ℝ) (hL : L < Real.log 2) : WeilPos L :=
  WLog2.weilPos hL (by
    have := Real.log_two_lt_d9
    simp only [WLog2.P]
    push_cast
    norm_num at this ⊢
    linarith)

/-- The round-4 targets, as corollaries. -/
theorem weilPos_le_65 (L : ℝ) (hL : L ≤ 13 / 20) : WeilPos L :=
  weilPos_lt_log2 L (by have := Real.log_two_gt_d9; norm_num at this; linarith)

theorem weilPos_le_60 (L : ℝ) (hL : L ≤ 3 / 5) : WeilPos L :=
  weilPos_le_65 L (by linarith)

end Fayed.Yoshida
