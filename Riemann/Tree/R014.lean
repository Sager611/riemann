import Mathlib

/-!
# Round 14 — composites of individually dead ingredients

| id | combination | status |
|----|-------------|--------|
| R14.1 | Weil windows ≤ 1.6 + RH to 3·10¹² + symbol growth | DEAD: no coercivity margin for low/high cross terms |
| R14.2 | de Bruijn–Newman (t ≥ 0.2) + GORZ Jensen hyperbolicity | DEAD: continuity cannot cross t = 0.2 |
| R14.3 | Levinson + Levinson–Montgomery left zeros of ζ' | DEAD: needs κ < 0.2915 left-zero density (critical_line_from_left_zeros) |
| R14.4 | Fayed's Perron kernel as a Weil pairing | DEAD: J is not a ξ'/ξ pairing |
| R14.5 | Báez-Duarte d_N + Korobov–Vinogradov Mertens bound | DEAD: near-1 regions only |
| R14.6 | Bagchi + Voronin + zero density | CIRCULAR |
| R14.7 | Hardy + symmetry + Hadamard stability | DEAD: line membership not perturbation-stable |
| R14.8 | Selberg CLT + Balazard–Saias–Yor | DEAD: CLT is distributional, BSY needs deterministic control |
-/

namespace Fayed.R14

/-- R14.3 arithmetic: if off-line zeros are at most twice the left zeros of ζ' (Levinson–Montgomery,
up to lower order) and those are at most κ N, then the critical-line proportion is ≥ 1 - 2κ. -/
theorem critical_line_from_left_zeros {N N₀ D κ : ℝ} (hN : 0 < N) (h1 : N - N₀ ≤ 2 * D)
    (h2 : D ≤ κ * N) : 1 - 2 * κ ≤ N₀ / N := by
  rw [le_div_iff₀ hN]
  nlinarith

end Fayed.R14
