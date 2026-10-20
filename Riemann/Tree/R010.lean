import Mathlib

/-!
# Round 10 — radical reformulations: change the object that carries the zeros

| id | reformulation | status |
|----|---------------|--------|
| R10.1 | sign changes of Hardy's Z vs N(T) | CIRCULAR: sign changes miss even-order zeros; N₀ = N is RH |
| R10.2 | Biane–Pitman–Yor: 2ξ as a Mellin moment of a positive law | DEAD: positive laws allow off-line Mellin zeros |
| R10.3 | Remmen amplitude 1/Ξ, unitarity | CIRCULAR: partial-wave positivity ≡ RH |
| R10.4 | Redheffer matrices, det R_n = M(n) | DEAD: spectral product is M(n) itself (formalized: Riemann/Deep/Redheffer) |
| R10.5 | arithmetic topology, torsion of Spec ℤ | DEAD: no arithmetic 3-manifold / flow |
| R10.6 | Busy Beaver / Matiyasevich encoding | DEAD: any halting bound decides RH |
| R10.7 | compensated partial sums P_N → ζ (Fayed (2),(15)), Hurwitz | DEAD: Hurwitz transfers ζ-zeros to P_N, not conversely |
| R10.8 | ξ as regularized determinant (Connes, Deninger) | CIRCULAR: self-adjointness is the missing RH input |
| R10.9 | Newton dynamics of ξ | DEAD: basins give no localization (newton_fixed_of_zero) |
| R10.10 | archimedean Iwasawa module | CIRCULAR: no canonical module |
-/

namespace Fayed.R10

/-- R10.7: Fayed's compensated partial sums P_N(s) = Σ_{n≤N} n^{-s} - N^{1-s}/(1-s). -/
noncomputable def PN (N : ℕ) (s : ℂ) : ℂ :=
  ∑ n ∈ Finset.Icc 1 N, (n : ℂ) ^ (-s) - (N : ℂ) ^ (1 - s) / (1 - s)

/-- R10.9: Newton map of f. -/
noncomputable def newtonMap (f : ℂ → ℂ) (z : ℂ) : ℂ := z - f z / deriv f z

/-- R10.9: zeros are fixed points of the Newton map — and that is all the dynamics gives for free. -/
theorem newton_fixed_of_zero {f : ℂ → ℂ} {z : ℂ} (hz : f z = 0) : newtonMap f z = z := by
  simp [newtonMap, hz]

end Fayed.R10
