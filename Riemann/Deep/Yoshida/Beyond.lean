/-
Riemann/Deep/Yoshida/Beyond.lean — Weil positivity beyond `log 2`: every window `L ≤ 3/4`, sorry-free.

Route: `BeyondCore` (Fourier form with the prime 2 and free cosines `cos(r jℓ/8)`, `j ≥ 9`, which integrate to
zero against `|F|²`), `BeyondBall`/`BeyondReflect` (`ℤ` fixed-point ball checker), and the generated
certificate `BeyondW75` at `a = 3/8` (3440 cells, Taylor order 40, margin 0.0016 in 2π units).
-/
import Riemann.Deep.Yoshida.BeyondW75

namespace Fayed.Yoshida

open Riemann.Deep.Weil

/-- **Weil positivity on every window `L ≤ 3/4`** (beyond Yoshida's `log 2 ≈ 0.693`; the prime `2` enters). -/
theorem weilPos_le_75 (L : ℝ) (hL : L ≤ 3 / 4) : WeilPos L :=
  Beyond.BeyondW75.weilPos (by
    simp only [Beyond.BeyondW75.C]
    push_cast
    linarith)

end Fayed.Yoshida
