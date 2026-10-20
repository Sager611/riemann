# Nyman–Beurling sufficiency — REPORT

Status: DONE. `Riemann/Deep/NymanBeurling/Suff.lean` (namespace `Fayed.NB`), no sorry.
`lake build Riemann.Deep.NymanBeurling.Suff` green; `#print axioms Fayed.NB.RH_of_NBApprox` →
[propext, Classical.choice, Quot.sound]. Target statements (rhoNB, NBApprox, RH_of_NBApprox) verbatim.

Key lemmas:
- `F s := mellin Q (-s)`, Q = {t}·1_{t≥1}; `differentiableAt_F` (Re s > 0).
- `zeta_eq_of_one_lt`: ζ(s) = s/(s-1) - s F(s), Re s > 1 (Abel summation via `LSeries_eq_mul_integral'` with f = 1).
- `zeta_eq`: same on {Re s > 0} \ {1} (identity theorem, reuses `Fayed.Mertens.isPreconnected_U`).
- `M_eq`: M(s) := ∫_0^∞ {1/x} x^{s-1} = 1/(1-s) + F(s), 0 < Re s < 1 (`mellin_comp_inv`, split (0,1) ∪ [1,∞)).
- `M_eq_zero`: ζ(s) = 0 ⇒ M(s) = 0.
- `hasMellin_rho`: Mellin of ρ_θ = θ^{s}M(s) - θ M(s) (`mellin_comp_mul_left`); `integral_rho_eq_zero` (ρ_θ = 0 on x > 1).
- `no_zero_of_NBApprox`: ∫_0^1 (1-f) x^{s-1} = 1/s, AM–GM pointwise bound |h|y ≤ (λh² + y²/λ)/2 replaces
  Cauchy–Schwarz; ε = |1/s|²/(B+1), B = ∫_0^1 x^{2σ-2}.
- `RH_of_NBApprox` via `Fayed.R7.RH_iff_QuasiRH_half` + `riemannZeta_ne_zero_of_one_le_re`.

Remaining gaps: none. Necessity half (RH ⇒ NBApprox) not attempted (out of scope).
