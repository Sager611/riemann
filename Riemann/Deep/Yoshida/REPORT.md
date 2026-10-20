# Yoshida deep dive: Weil positivity on small windows (Lean)

**Verdict: strict improvement, not log 2.** Sorry-free, axioms = propext, Classical.choice, Quot.sound:
- `Fayed.Yoshida.weilPos_le_third : ∀ L ≤ 1/3, WeilPos L` (f supported in [−1/6, 1/6]). Window 1/3, up from 1/16; log 2 ≈ 0.693.
- Kept from round 1: `weilPos_small` (L ≤ 1/16, `Small.lean`), `weilPos_iff_primeFree` (exact prime-free target below log 2).
- New exact identity `weilQ_re_eq_W` (below log 2): `Re Q(f) = (1/2π) ∫ |F(r)|² (2π μ(r) + J(r)) dr`, with
  `J(r) = ∫_{−2a}^{2a} 2cosh(x/2) cos(rx) dx = 4((1/2) sinh a cos 2ar + r cosh a sin 2ar)/(1/4 + r²)`.
  The pole term is no longer bounded crudely; it becomes part of a single weight on the Fourier side.
- Build: `lake build Riemann.Deep.Yoshida.Windows` (≈ 3 min; peak memory ≈ 5.8 GB, about the same as importing Mathlib, because the certificate is elaborated sequentially). Not imported by Riemann.lean.

## Files
| file | content |
|---|---|
| `Taylor.lean` | alternating Taylor bounds for cos and sin on x ≥ 0, every degree, by induction (comparing derivatives) |
| `PoleKernel.lean` | `J` closed form (FTC), evenness, Lipschitz bound (constant 8a² cosh a), tail bound; `pole_re_eq_J`: Re pole = (1/2π)∫|F|²J (Fourier inversion of the autocorrelation + Fubini) |
| `MuBounds.lean` | `2πμ(r) ≥` partial sum of the digamma series; γ < 0.5792 (H₂₅₆ − log 256); log π < 1.1448; Taylor bounds for exp |
| `Hybrid.lean` | `weilQ_re_eq_W`; layer cake with a step function over a list of cells; `weilPos_of_cert` (certificate ⇒ WeilPos); per-cell and tail lemmas |
| `Third.lean` | GENERATED certificate, 139 cells, s = 9/20, R = 18, margin 0.082 (in W units) |
| `gen_cells.py` | generator (Python stdlib, exact rationals; checks every inequality before writing Lean) |
| `Windows.lean` | `weilPos_third`, `weilPos_le_third` |

## Proof (f supported in [−a, a], N = ∫|f|², φ = |F|²)
1. Pole term = ∫ k(x) 2cosh(x/2) dx, k = autocorrelation (lives in [−2a, 2a]); Fourier inversion + Fubini give (1/2π)∫ φ J.
2. Hence Re Q = (1/2π)∫ φ W, W = 2πμ + J. Plancherel ∫φ = 2πN; pointwise φ ≤ 2aN.
3. If W(r) ≥ s − t(r) with t ≥ 0 a step function, then ∫φW ≥ 2πN s − 2aN ∫t. Positive when 2a∫t ≤ 2πs.
4. Per cell lo ≤ |r| < hi: μ(r) ≥ μ(lo) (monotone), J(r) ≥ J(mid) − Lipschitz·half-width, J(mid) from Taylor bounds of cos/sin, μ(lo) from 50 digamma-series terms. For |r| ≥ 18: μ(18) minus a crude bound on J.
5. The cells and s come from a dynamic program in `gen_cells.py` and are checked by `norm_num` in Lean.

## Why it stops near 1/3 (numerics in .scratch/yoshida2/, not tracked)
- With this weight and the cap φ ≤ 2aN, the best possible bound (exact μ, infinitely fine cells) turns negative at about L ≈ 0.355 (`hyb.py`). At 1/3 the certified margin is small (0.082 out of an ideal 0.28) because cells lose some sharpness.
- **Next rung, about 3/8 (numerically checked, not formalised).** For any b ≥ a, the pole also equals (1/2π)∫φ J_b, because the autocorrelation vanishes beyond 2a. So W = 2πμ + Σ λ_k J_{b_k} with Σλ_k = 1 is allowed. With 3 values of b and about 210 cells, a certified-style margin of 0.043 is reached at L = 3/8 (`genL_probe.py`); L = 0.4 still fails (−0.12). Missing Lean pieces: upper bounds for J (some λ < 0), range reduction for sin and cos at arguments up to about 20, and a telescoping bound for the series tail.
- **Toward log 2.** Without the cap the dual problem is exact: find H equal to 2cosh(x/2) on [−2a, 2a], free outside, with μ + Ĥ/2π ≥ 0. A linear program finds such H at L = 0.5 (margin 0.0026), but H has about 300 nodes (`lp.py`). The true lowest eigenvalue at log 2 is 1.3e-3. A fixed-size cell certificate cannot close that gap; a proof at log 2 needs the prolate/Sonin structure (Connes–Consani) or a Galerkin certificate with tail bounds (Route A of round 1).
