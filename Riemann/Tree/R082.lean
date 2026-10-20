/-!
# Round 82 — which small-window positivity technique scales beyond 1/16?

| technique | window reached | Lean feasibility |
|-----------|----------------|------------------|
| Yoshida 1992 (Fourier inequalities + Γ-factor estimates) | L < log 2 | low |
| Connes–Consani 2021 prolate concentration + archimedean trace formula | potentially far beyond log 2 | low |
| Paley–Wiener / Bernstein from F(±i/2) | no gain (F can vanish at ±i/2) | low |
| Cauchy–Schwarz + pole lower bound | no gain (no uniform pole lower bound) | low |
| certified smallest eigenvalue of the window operator | up to log 2, certificate-dependent | low |
| positive-definite kernel decomposition | up to log 2 with spectral certification | low |
| Weil / Burnol archimedean positivity | L < log 2 | low |
| Taylor expansion in L | incremental beyond 1/16 | low |

Every scalable route is heavy; the Taylor route (being tried by the agent: Taylor.lean, MuBounds.lean,
PoleKernel.lean) is the realistic incremental one.
-/
