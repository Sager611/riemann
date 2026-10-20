# Certified first zero of ζ — REPORT

**Result: done, no sorry.** `Fayed.FirstZero.first_zero : ∃ t : ℝ, 14.13 < t ∧ t < 14.14 ∧ riemannZeta (1/2 + t·I) = 0`
(Riemann/Deep/FirstZero/FirstZero.lean). `#print axioms` = [propext, Classical.choice, Quot.sound].
Built with `lake build Riemann.Deep.FirstZero.FirstZero` (green). Not yet imported from Riemann.lean (outside my edit scope; orchestrator: add `import Riemann.Deep.FirstZero.FirstZero`).

## Route
1. Λ(1/2+it) is real (Fayed.Hardy.Lam_conj) and continuous; Λ = ζ·Γℝ.
2. ζ at t₁ = 14.1301, t₂ = 14.1399 by order-8 Euler–Maclaurin, N = 6: |error| ≤ 4.6e-4 vs |ζ| ≈ 3.7e-3.
3. Γℝ(s₂)/Γℝ(s₁) within 0.3 of 1 (no Stirling needed): Euler limit GammaSeq, telescoped n^{iδ}, paired factors 1 + O(δ|z|/(j+1)²).
4. Re(ζ₂·conj ζ₁·ρ) ≤ −4.9e-6 < 0 ⇒ Λ(s₁)Λ(s₂) < 0 ⇒ IVT ⇒ Λ = 0 ⇒ ζ = 0.

## Files (Riemann/Deep/FirstZero/)
| file | content |
|---|---|
| Zeta0.lean | port of PNT+ first-order formula `PNT.Zeta0EqZeta` (ζ = ζ₀ N on Re s > 0, s ≠ 1), via zeta23 copy |
| EM.lean | `zeta_EM8`: order-8 Euler–Maclaurin with explicit remainder; `abs_BP8_le` via Bernstein coefficients |
| GammaRatio.lean | `Gamma_ratio_le`: ‖Γ(z+iδ)/Γ(z) − 1‖ ≤ 1/4 (Re z = 1/4, 7 ≤ Im z ≤ 7.1, δ ≤ 1/200) |
| Ball.lean | complex ball arithmetic (add/mul/real-scale), `cis_ball` (π-reduction + degree-19 Taylor), `inv_sqrt_ball` |
| Crit.lean | p^{-s} = p^{-1/2}e^{-it log p}, log 2/3/5 to 1e-9, N = 6 main-term identity `em_main_eq`, `zeta_ball` |
| Cert.lean | generated (.scratch/firstzero/gen_cert.py): `zeta_t1`, `zeta_t2` enclosures, all side conditions by norm_num |
| FirstZero.lean | `Gammaℝ_ratio`, `re_mul_conj_neg`, `first_zero` |

## Gaps
None. No numeric certificate Props needed (ladder step 2 unused). Slack: EM error bound 4.6e-4 vs |ζ| ≈ 3.7e-3; final sign margin ≈ 33 %.

## Prior art
- zeta23 (Anthropic): `rg '14\.13|first zero|zetaZero'` finds nothing; it has the PNT+ first-order Euler–Maclaurin (reused here) and assumes Gamma facts as hypotheses.
- Mathlib: `riemannZetaZeros` (discreteness) only; no explicit zero, no Stirling for complex Γ, digamma without bounds.
