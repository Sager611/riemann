# Weil positivity ⇒ RH in Lean (WeilEF deep dive)

**Status 2026-10-22:** the archTerm bug below is fixed in Riemann/Deep/Weil/Defs.lean (the integral is now parenthesised); Riemann.Deep.WeilEF.Suff is imported by Riemann.lean; Bug.lean (which records the old bug) is not imported.

**Verdict: top rung reached.** Weil's criterion is proved in Lean, both directions, sorry-free, axioms = propext, Classical.choice, Quot.sound only:
`Fayed.WeilSuff.weil_criterion : RiemannHypothesis ↔ ∀ L > 0, WeilPosFix L` (Suff.lean).
**Bug found in Riemann/Deep/Weil/Defs.lean:** `archTerm` reads `(1/2π) * ∫ r, h r * Re ψ - g 0 * log π`; the integral swallows `- g 0 * log π`, the integrand is then not integrable, and Lean's junk value makes the term 0. Proved: `Fayed.WeilSuff.Bug.not_weilPos : 0 < L → ¬ WeilPos L` (Bug.lean). So every window Prop in Defs.lean (Yoshida, Zhu, Liu, Mori, CertWindow) is false as written, `WeilCriterion` there is equivalent to ¬RH, and Desogus/Paper.lean's `WeilPosOdd` uses the same broken `weilQ`. The round-55 audit missed this.
**Fix (orchestrator, one line):** in Defs.lean write `(1 / (2 * Real.pi) : ℂ) * (∫ r : ℝ, ft g r * (...).re) - g 0 * ...`, i.e. parenthesise the integral. Then `archTerm = archTermFix` by `rfl`, and `RH_of_weilPos (fun _ => rfl) h` is the brief's theorem verbatim.

## Ladder
1. Port (done): the zeta23 explicit-formula closure (59 files, ~18.8k lines, adapted from anthropics/formal-math zeta23, Apache-2.0; LICENSE+NOTICE copied) lives in `Z23/` under module prefix `Riemann.Deep.WeilEF.Z23`. Six small fixes for Mathlib drift (prod_le_prod₀, prod_le_one₀, logDeriv_fun_prod, three `erw` for logDeriv_mul, one `module` header). Port script: .scratch/weilef/port.py; build loop: .scratch/weilef/buildloop.sh.
2. Explicit formula in our normalisation (Bridge.lean): `weilQFix f = literatureRHS (autocorr f) = Σ_ρ m_ρ F(γ_ρ)·conj F(conj γ_ρ)`, summable, for every `TestFn`; `F = ∫ f e^{izu}`, `γ_ρ = (ρ−1/2)/i`, zeros and multiplicities from Mathlib's ζ (`analyticOrderAt`). The `−g(0) log π` comes from Fourier inversion.
3. Sufficiency (Core.lean + Suff.lean), new argument, not in zeta23:
   - `D τ = Σ w_ρ e^{(ρ−1/2)τ}`, `w_ρ = m_ρ F(γ_ρ) conj F(conj γ_ρ)`. Testing `f + c·f(·−τ)` gives `Q = (1+|c|²)D 0 + c D τ + conj c·D(−τ)`; the zero reflection `ρ ↦ 1−conj ρ` gives `D(−τ) = conj D τ`; positivity for all `c` gives `‖D τ‖ ≤ Re D 0`.
   - Core lemma (no zeta): a bounded exponential sum has zero weight at any exponent with Re > 0 that has no exponent of larger real part within ordinate distance η. Proof: Laplace transform is holomorphic on Re s > 0 and ≤ M/Re s; equals Σ w/(s−p) on Re s > 1/2 (Fubini); identity theorem on a convex rectangle; then x·L(p₀+x) → 0 but → w₀.
   - Pick ρ₀ = largest real part in a finite ordinate window (finite_window), η from finitely many neighbours; test function = bump·e^{−i t₀ u} makes w_{ρ₀} = m·(∫bump e^{bu})(∫bump e^{−bu}) > 0. Contradiction.
   - Necessity: under RH each term is m|F(γ)|² ≥ 0.

## Files
Core.lean (329 lines), Bridge.lean (136), Suff.lean (376), Bug.lean (179), Z23/ (port). Top module: `Riemann.Deep.WeilEF.Bug` (imports all). Not yet imported by Riemann.lean; first build of the Z23 closure ≈ 15 min, then seconds.

## Significance
This is the textbook Weil criterion (Weil 1952, Bombieri 2000), not progress on RH: it moves RH to "positivity on all windows", which R3.1 already found to be doubly-exponentially hard per window. New: a sorry-free Lean proof of the criterion (zeta23 has the explicit formula but not the criterion), and the Defs.lean bug.
