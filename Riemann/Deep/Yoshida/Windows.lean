/-
Riemann/Deep/Yoshida/Windows.lean — Weil positivity on every window `L ≤ 11/20`, sorry-free.

Route (`Hybrid.lean`): below `log 2`, `Re Q(f) = (1/2π) ∫ |F(r)|² (2π μ(r) + J(r)) dr`.
Round 2 (`Third.lean`, window 1/3) bounded `∫ |F|² t` by the cap `|F|² ≤ 2aN`. Round 3
(`Moment.lean`, `Reflect.lean`, `Half.lean`) bounds it through the Taylor moments of `F`, a Bessel
inequality and a rational PSD certificate; all numerics are checked by `decide +kernel`.
-/
import Riemann.Deep.Yoshida.Half
import Riemann.Deep.Yoshida.W55

namespace Fayed.Yoshida

open Riemann.Deep.Weil

/-- **Weil positivity on the window `1/2`.** -/
theorem weilPos_half : WeilPos (1 / 2) := by
  have h := Half.weilPos
  have e : (2 * ((Half.P.a : ℚ) : ℝ)) = 1 / 2 := by
    simp only [Half.P]
    norm_num
  rw [e] at h
  exact h

/-- **Yoshida-type window**: Weil positivity for every window `L ≤ 1/2`. -/
theorem weilPos_le_half (L : ℝ) (hL : L ≤ 1 / 2) : WeilPos L :=
  WeilPos.mono hL weilPos_half

/-- **Weil positivity on the window `11/20`.** -/
theorem weilPos_55 : WeilPos (11 / 20) := by
  have h := W55.weilPos
  have e : (2 * ((W55.P.a : ℚ) : ℝ)) = 11 / 20 := by
    simp only [W55.P]
    norm_num
  rw [e] at h
  exact h

/-- **Largest window so far**: Weil positivity for every `L ≤ 11/20 = 0.55`. -/
theorem weilPos_le_55 (L : ℝ) (hL : L ≤ 11 / 20) : WeilPos L :=
  WeilPos.mono hL weilPos_55

/-- Round-2 statement, now a corollary. -/
theorem weilPos_le_third (L : ℝ) (hL : L ≤ 1 / 3) : WeilPos L :=
  weilPos_le_half L (by linarith)

end Fayed.Yoshida
