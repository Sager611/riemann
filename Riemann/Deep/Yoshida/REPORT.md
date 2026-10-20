# Yoshida deep dive: Weil positivity on small windows (Lean)

**Verdict: partial.** Sorry-free, axioms = propext, Classical.choice, Quot.sound:
- `Fayed.Yoshida.weilPos_small : ∀ L ≤ 1/16, WeilPos L` (f supported in [−1/32, 1/32]). Window 1/16 = 0.0625, not log 2 ≈ 0.693.
- `Fayed.Yoshida.weilPos_iff_primeFree (hL : L < log 2)`: WeilPos L ⟺ 0 ≤ 2 Re(F(i/2)·conj F(−i/2)) + ∫|F(r)|² μ(r) dr for every test f (F = paperFT f, μ = (Re ψ(1/4+ir/2) − log π)/2π). This is the exact target for the full theorem.
- File: `Riemann/Deep/Yoshida/Small.lean` (472 lines). Builds with `lake build Riemann.Deep.Yoshida.Small`. Not imported by Riemann.lean.

## Proof at 1/16 (layer cake; N = ∫|f|²)
1. Pole term ≥ −‖F(i/2) − F(−i/2)‖² ≥ −(2a∫|f|)² ≥ −8a³N.
2. Plancherel ∫|F|² = 2πN, from Fourier inversion of the autocorrelation at 0 (zeta23 `paper_inversion`).
3. |F(r)|² ≤ (∫|f|)² ≤ 2aN (Cauchy–Schwarz written as AM–GM), so ∫_{|r|<R}|F|² ≤ 4aRN.
4. μ is even, increasing in |r|, μ ≥ μ(0) > −1, and μ(R) ≥ (log(R/2π) − 20/R²)/2π (vertical Stirling). All of this was already proved in the ported zeta23 `GammaFacts`.
5. Four steps R = 4, 21/2, 69/4, 113/4 give Q/N ≥ 0.166 > 0 at a = 1/32.

## Why this does not reach log 2 (numerics in .scratch/yoshida/)
- Every bound of the form "|F(r)|² ≤ b(r), total mass 2πN" (bathtub) stops short of log 2 (`lp2.py`): |F|² ≤ 4aN: 0.088; Cauchy–Schwarz: 0.175; plus slope and pole: 0.278; plus exact projection onto constants: 0.304. With crude Stirling and 20 steps: 0.096 (`steps.py`).
- The true λ_min at log 2 is only 1.33e-3·‖f‖² (Weil/REPORT.md). The worst f is close to the first prolate function, by a prolate estimate (c ≈ 2.2) about 90% of its |f̂|² mass sits where μ < 0, and only the pole term saves it. So any proof has to use the rank-2 pole term together with the prolate concentration of f̂ (Yoshida; Connes–Consani), not pointwise bounds.

## Next rung (open)
- Prove the strict window log 2 − ε. Route A: a Lean-checked Galerkin certificate for εI + C ⪰ 0 (cert.py method, no primes, R ≈ 7, small basis). Needs Legendre/Bessel tail bounds and ψ quadrature in Lean, which is heavy. Route B: Connes–Consani prolate/Sonin argument, which needs prolate spheroidal theory that Mathlib lacks.
- Cheap gains, of little value: more steps plus a series bound for μ(0) push 1/16 to about 0.1–0.15.
