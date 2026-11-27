import Riemann.Deep.LeanGale.Main
import Riemann.Deep.LeanGale.BallGeo

/-!
# LeanGale T6 (3/3): `SimplexPartitionBound` is a theorem

`simplexPartitionBound`: Voronov's literature input `SimplexPartitionBound` (Main.lean; β_{n+1}(Bⁿ) ≤ σₙ for
n ≥ 2, B the closed ball of diameter 1) holds, by the regular-simplex partition (BallGeo
`beta_closedBall_le_sigma`, BallAlg `psi_nonneg`). `gale_ball_fails'`: the ball form of Gale's conjecture fails in
every dimension n ≥ 3 with `Hinrichs2002` as the only hypothesis.
-/

namespace LeanGale

/-- **The regular-simplex partition bound** (Voronov arXiv 2610.28828, main_en.tex l. 96–107; Gale 1953), proved:
β_{n+1}(B) ≤ σₙ for every n ≥ 2, B the closed ball of diameter 1 in ℝⁿ. -/
theorem simplexPartitionBound : SimplexPartitionBound := fun _ hn =>
  beta_closedBall_le_sigma (by omega)

/-- **The "ball is extremal" form of Gale's conjecture fails in every dimension n ≥ 3**, given only
`Hinrichs2002` (for n > 1000). -/
theorem gale_ball_fails' (hHin : Hinrichs2002) : ∀ n : ℕ, 3 ≤ n → ¬ GaleConjectureBall n :=
  gale_ball_fails hHin simplexPartitionBound

end LeanGale
