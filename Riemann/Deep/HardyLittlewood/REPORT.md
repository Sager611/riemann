# Hardy–Littlewood criterion, sufficiency — `Suff.lean`

- `Fayed.HL.RH_of_hlBound`: (∀ ε>0, ∃ C, ∀ x≥1, |hlF x| ≤ C x^{-1/4+ε}) → RiemannHypothesis. Sorry-free; axioms: propext, Classical.choice, Quot.sound. `lake build Riemann.Deep.HardyLittlewood.Suff` green.
- **Definition change (needed):** the briefed `hlF x = ∑' (μ n/n)·exp(-x/n²)` is not absolutely summable (terms ≈ μ(n)/n), so Lean's `tsum` makes it ≡ 0 and the theorem would just restate RH. We use `hlF x = ∑' (μ n/n)(exp(-x/n²) − 1)`. Without extra assumptions this equals Σ_{k≥1} (−x)^k/(k! ζ(2k+1)). It also equals the classical conditionally convergent Möbius series, because Σ μ(n)/n = 0. Neither identity is formalized; both are stated in a comment.
- Proof route:
  1. hlF' = −g, where g = Σ μ(n) n⁻³ e^{−x/n²}.
  2. `hasSum_mellin` gives M[g](z) = Γ(z) L(μ, 3−2z) for 0 < Re z < 1.
  3. Integration by parts (`integral_Ioi_mul_deriv_eq_deriv_mul`) gives M[hlF](s) = M[g](s+1)/s = Γ(s)L(μ,1−2s) for −1 < Re s < −1/2. This step uses only bounds that need no hypothesis: |hlF x| ≤ S_b x^b for 0 < b ≤ 1, and |g x| ≤ ζ(2) x^{−1/2}.
  4. The hypothesis makes M[hlF] holomorphic on −1 < Re s < 1/4.
  5. The identity theorem, applied on two convex half-strips that avoid s = 0, gives ζ(1−2s)M(s) = Γ(s) ≠ 0 for 0 < Re s < 1/4. So ζ(w) ≠ 0 for 1/2 < Re w < 1. Finish with `R7.RH_iff_QuasiRH_half`.
- Differences from Riesz: Riesz's series identity holds directly on its strip. Here the Möbius series only converges for Re s < 0, and on that region the term-by-term Mellin integrals diverge, so the integration-by-parts step (3) is needed.
- Residual: the power-series identity is not formalized, and the necessity direction is not done.
