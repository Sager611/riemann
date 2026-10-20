import Riemann.Fayed.Theorem1

/-!
# Round 6 — "Fayed + X": the root identity combined with one global input

| id | X | status |
|----|---|--------|
| R6.1 | −ζ'/ζ in the Perron step, Λ ≥ 0 | DEAD: coefficient positivity is lost for t ≠ 0; residues restate zeros |
| R6.2 | 1/ζ in the Perron step (Möbius) | DEAD: residues restate ρ ↦ 1-ρ̄ symmetry |
| R6.3 | mean square in t (Montgomery–Vaughan) | DEAD: zeros are measure-zero samples |
| R6.4 | Hadamard product / Guinand | CIRCULAR: sign-indefinite zero sums; positivity ≡ RH |
| R6.5 | Jensen / Carathéodory–Borel | DEAD: counts zeros, cannot locate real parts |
| R6.6 | approximate functional equation, σ-sensitive χ | DEAD: at a zero S + χT + E = 0 absorbs χ |
| R6.7 | complex p, p-families | DEAD: I(u,s₀) ≡ 0 for every p (Ifun_zero) |
| R6.8 | Hurwitz ζ(s,a), a-averages | CIRCULAR: a-average vanishes; nothing singles out a = 1 |
| R6.9 | small u → 0⁺ | DEAD: D_γ is bounded away from 0 (no_small_u_in_Dγ) |
| R6.10 | log ζ with 3-4-1 inside J | DEAD: Euler product only for Re s > 1 |
-/

open Filter Topology Set
open scoped ComplexConjugate

namespace Fayed.R6

/-- R6.7: at a zero s₀ (so ζ(1 - s̄₀) = 0 too) Fayed's I vanishes for EVERY complex p and u,
so no family of p's carries information. -/
theorem Ifun_zero {s : ℂ} (hs : riemannZeta s = 0) (hdual : riemannZeta (1 - conj s) = 0)
    (p u : ℂ) : u ^ (s + p) * riemannZeta s - u ^ (1 - conj s + p) * riemannZeta (1 - conj s) = 0 := by
  rw [hs, hdual]; ring

/-- R6.6: at a zero the approximate functional equation is a bare cancellation. -/
theorem afe_at_zero {S T χ E ζs : ℂ} (hafe : ζs = S + χ * T + E) (hz : ζs = 0) :
    S + χ * T + E = 0 := hafe ▸ hz

/-- R6.9: D_γ stays above 1/2, so no sequence in D_γ tends to 0. -/
theorem no_small_u_in_Dγ {γ : ℝ} (hγ : 0 < γ) (M₀ : ℕ) :
    ¬ ∃ u : ℕ → ℝ, Tendsto u atTop (𝓝 0) ∧ ∀ n, u n ∈ Dγ γ M₀ := by
  rintro ⟨u, hu, hD⟩
  have hev : ∀ᶠ n in atTop, u n < 1 / 2 := hu (Iio_mem_nhds (by norm_num))
  obtain ⟨n, hn⟩ := hev.exists
  obtain ⟨M, hM⟩ := mem_iUnion.1 (hD n)
  have h1 := hM.1
  simp only [aγ] at h1
  have : (0 : ℝ) ≤ ((M₀ + M : ℕ) : ℝ) := Nat.cast_nonneg _
  linarith

end Fayed.R6
